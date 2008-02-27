! $Id: poisson.f90,v 1.42 2008-02-27 12:35:06 wlyra Exp $

!
!  This module solves the Poisson equation
!    (d^2/dx^2 + d^2/dy^2 + d^2/dz^2 - h) f = RHS(x,y,z)
!  [which for h/=0 could also be called inhomogenous nonuniform Helmholtz
!  equation] for the function f(x,y,z).
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

  real :: kmax=0.0, expand_factor=2.
  logical :: lrazor_thin=.false., lsemispectral=.false., lklimit_shear=.false.

  include 'poisson.h'

  namelist /poisson_init_pars/ &
      lsemispectral, kmax, lrazor_thin, lklimit_shear, expand_factor

  namelist /poisson_run_pars/ &
      lsemispectral, kmax, lrazor_thin, lklimit_shear, expand_factor

  contains

!***********************************************************************
    subroutine initialize_poisson()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  18-oct-07/anders: adapted
!
      if (lrazor_thin) then
        if (nzgrid/=1) then
          if (lroot) print*, 'inverse_laplacian_fft: razon-thin approximation only works with nzgrid==1'
          call fatal_error('inverse_laplacian_fft','')
        endif
      endif
!
!  Limit the wavenumber to the maximum circular region that is always available
!  in k-space. The radial wavenumber kx changes with time due to shear as
!
!    kx = kx0+qshear*Omega*t*ky
!
!  Considering the available (kx,ky) space, it turns slowly from a square to a
!  parallellogram (the hole for kx,ky<1 is ignored here):
!
!       - - - -                  - - - -
!      |       |               /       /
!      |       |    ->       /       /
!      |       |           /       /
!      |       |         /       /
!       - - - -          - - - -
!     
!  To make the gravity force isotropic at small scales one can limit k to
!  the largest circular region that is present in both the square and the
!  parallellogram. The circle has radius kmax=kNy/sqrt(2). See Gammie (2001).
!
      if (lklimit_shear) kmax = kx_ny/sqrt(2.0)
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
      if (lsemispectral) then
!
        call inverse_laplacian_semispectral(phi)
!
      else
        if (lcylindrical_coords) then
          if (lmpicomm) then
            call inverse_laplacian_fft_cyl_mpi(phi)
          else
            call inverse_laplacian_fft_cyl(phi)
          endif
        else
          call inverse_laplacian_fft(phi)
        endif
      endif
!
    endsubroutine inverse_laplacian
!***********************************************************************
    subroutine inverse_laplacian_fft(phi)
!
!  Solve the Poisson equation by Fourier transforming on a periodic grid.
!  This method works both with and without shear.
!
!  15-may-2006/anders+jeff: coded
!
      real, dimension (nx,ny,nz) :: phi
!
      real, dimension (nx,ny,nz) :: b1
      real :: k2
      integer :: ikx, iky, ikz
!
!  Identify version.
!
      if (lroot .and. ip<10) call cvs_id( &
        "$Id: poisson.f90,v 1.42 2008-02-27 12:35:06 wlyra Exp $")
!
!  The right-hand-side of the Poisson equation is purely real.
!
      b1 = 0.0
!
!  Forward transform (to k-space).
!
      if (lshear) then
        call fourier_transform_shear(phi,b1)
      else
        call fourier_transform(phi,b1)
      endif
!
!  Solve Poisson equation.
!
      do ikz=1,nz; do iky=1,ny; do ikx=1,nx
        if ((kx_fft(ikx)==0.0) .and. &
            (ky_fft(iky)==0.0) .and. (kz_fft(ikz)==0.0) ) then
          phi(ikx,iky,ikz) = 0.0
          b1(ikx,iky,ikz) = 0.0
        else
          if (.not. lrazor_thin) then
            if (lshear) then
!
!  Take into account that the Fourier transform has been done in shearing
!  coordinates, and that the kx of each Fourier mode is different in the normal
!  frame. The connection between kx0 (in the shearing frame) and the actual kx
!  is
!    kx = kx0 + qshear*Omega*t*ky
!  Writing deltay/2 = qshear*Omega*Lx/2*t, one gets the expression
!    kx = kx0 + deltay/Lx*ky
!  The parallel Fourier transform returns the array in a tranposed state, so
!  must be able to identify the x-direction in order to take shear into account.
!  (see the subroutine transform_fftpack_shear in Mpicomm for details).
!
              if (nzgrid/=1) then ! Order (kz,ky',kx)
                k2 = (kx_fft(ikz+ipz*nz)+deltay/Lx*ky_fft(iky+ipy*ny))**2 + &
                      ky_fft(iky+ipy*ny)**2 + kz_fft(ikx)**2
              else                ! Order (kx,ky',kz)
                k2 = (kx_fft(ikx)+deltay/Lx*ky_fft(iky+ipy*ny))**2 + &
                      ky_fft(iky+ipy*ny)**2 + kz_fft(ikz+ipz*nz)**2
              endif
!  The ordering of the array is not important here, because there is no shear!
            else
              k2 = kx_fft(ikx)**2 + ky_fft(iky+ipy*ny)**2 + kz_fft(ikz+ipz*nz)**2
            endif
!
!  Solution of Poisson equation.
!
            phi(ikx,iky,ikz) = -phi(ikx,iky,ikz) / k2
            b1(ikx,iky,ikz)  = - b1(ikx,iky,ikz) / k2
!
          else
!
!  Razor-thin approximation. Here we solve the equation
!    del2Phi=4*pi*G*Sigma(x,y)*delta(z)
!  The solution at scale k=(kx,ky) is
!    Phi(x,y,z)=-(2*pi*G/|k|)*Sigma(x,y)*exp[i*(kx*x+ky*y)-|k|*|z|]
!
            if (lshear) then
              k2 = (kx_fft(ikx)+deltay/Lx*ky_fft(iky+ipy*ny))**2+ky_fft(iky+ipy*ny)**2
            else
              k2 = kx_fft(ikx)**2+ky_fft(iky+ipy*ny)**2
            endif

            print*,phi(ikx,iky,ikz),b1(ikx,iky,ikz),ikx,iky
            
            phi(ikx,iky,ikz) = -.5*phi(ikx,iky,ikz) / sqrt(k2)
            b1(ikx,iky,ikz)  = -.5* b1(ikx,iky,ikz) / sqrt(k2)
          endif
        endif
!
!  Limit |k| < kmax
!
          if (kmax>0.0) then
            if (sqrt(k2)>=kmax) then
              phi(ikx,iky,ikz) = 0.
              b1(ikx,iky,ikz) = 0.
            endif
          endif
      enddo; enddo; enddo
!
!  Inverse transform (to real space).
!
      if (lshear) then
        call fourier_transform_shear(phi,b1,linv=.true.)
      else
        call fourier_transform(phi,b1,linv=.true.)
      endif
!
    endsubroutine inverse_laplacian_fft
!***********************************************************************
    subroutine inverse_laplacian_fft_cyl(phi)
!
!  Solve the 2D Poisson equation in cylindrical coordinates
!  by transforming to a periodic cartesian grid before 
!  Fourier transforming 
!
!  This method works only for serial runs
!
! 01-12-07/wlad: coded
!
      real, dimension (nx,ny,nz) :: phi
      real, dimension (2*nx,2*ny) :: nphi,nb1
!
      real, dimension(nx) :: rad
      real, dimension(ny) :: tht
!
      real, dimension(2*nx) :: xc,kkx_fft
      real, dimension(2*ny) :: yc
      real, dimension(2*nygrid) :: kky_fft
!
      real :: r0,rn,x0,xn,y0,yn,dxc,dyc,dr,dth,radius,rr
      real :: delr,delp,theta,fr,fp,p1,p2,p3,p4,interp_pot
      real :: distx,disty,fx,fy,delx,dely,xp,yp,k2,dxc1,dyc1
      real :: Lxn,Lyn,theta0,theta1,dr1,dth1
      integer :: nr,nth,ir,im,ix1,ix2,iy1,iy2,nnxgrid,nnygrid
      integer :: i,ir1,ir2,ip1,ip2,nnx,nny,j,ikx,iky
      integer :: nrgrid,nthgrid,nnghost
!
! transform the grid to cartesian
! keep the notation consistent.
!
      rad=x(l1:l2)     ; tht=y(m1:m2)
      nr =nx           ; nth=ny
      r0=rad(1)        ; theta0=xyz0(2)+.5*dth 
      rn=rad(nr)       ; theta1=xyz1(2)-.5*dth
      nrgrid=nxgrid    ; nthgrid=nygrid
      dr=dx            ; dth=dy
      dr1=1./dr        ; dth1=1./dth
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
              delp=theta-tht(ip1) + 2*pi
            else
              delp=theta-tht(ip1)
            endif
            if (ip2==nthgrid+1) ip2=1
!
! Bilinear interpolation
!
            p1=phi(ir2,ip1,nnghost);p2=phi(ir1,ip1,nnghost)
            p3=phi(ir2,ip2,nnghost);p4=phi(ir1,ip2,nnghost)
!
            fr=delr*dr1
            fp=delp*dth1
!
            nphi(i,m)=fr*fp*(p3+p2-p1-p4) + fr*(p1-p2) + fp*(p4-p2) + p2
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
!              if (ipy==0) print*,kkx_fft(ikx),kky_fft(iky+ipy*nny)
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
!  Convert back to cylindrical
!
      do ip=1,Nth
        do ir=1,Nr
!
          xp=rad(ir)*cos(tht(ip))
          yp=rad(ir)*sin(tht(ip))
!
          distx = xp - x0
          disty = yp - y0
!
          ix1 = floor(distx*dxc1)+1 ; ix2 = ix1+1
          iy1 = floor(disty*dyc1)+1 ; iy2 = iy1+1
!
          if (ix1 .lt.  1)      call fatal_error("","ix1 lt 1")
          if (iy1 .lt.  1)      call fatal_error("","iy1 lt 1")
          if (ix2 .gt. nnxgrid) call fatal_error("","ix2 gt nnxgrid")
          if (iy2 .gt. nnygrid) call fatal_error("","iy2 gt nnygrid")
!
          delx=xp-xc(ix1);fx=delx*dxc1
          dely=yp-yc(iy1);fy=dely*dyc1
!
! Bilinear interpolation
!
          p1=nphi(ix2,iy1);p2=nphi(ix1,iy1)
          p3=nphi(ix2,iy2);p4=nphi(ix1,iy2)
!
          interp_pot=fx*fy*(p3+p2-p1-p4) + fx*(p1-p2) + fy*(p4-p2) + p2
!
          do n=1,nz
            phi(ir,ip,n)=interp_pot
          enddo
!
        enddo
      enddo
!
    endsubroutine inverse_laplacian_fft_cyl
!***********************************************************************
    subroutine inverse_laplacian_fft_cyl_mpi(phi)
!
!  Solve the 2D Poisson equation in cylindrical coordinates
!  by transforming to a periodic cartesian grid before 
!  Fourier transforming. 
!
!  This is for parallel runs, but it still seems to involve an 
!  awful lots of communication. 
!
!  20-feb-08/wlad : coded
!
      use Mpicomm
!
      real, dimension (nx,ny,nz)  :: phi
      real, dimension (2*nx,2*ny) :: nphi,nb1
!
      real, dimension(nx) :: rad
      real, dimension(ny) :: tht
!
      real, dimension(nygrid)   :: theta_serial
!
      real, dimension(2*nx)     :: xc,kkx_fft
      real, dimension(2*ny)     :: yc
      real, dimension(2*nygrid) :: kky_fft,yserial
!
      real   , dimension (0:ncpus-1,nx*ny)     :: phi_send,phi_recv
      integer, dimension (0:ncpus-1,nx*ny)     :: sri_send,sri_recv
      integer, dimension (0:ncpus-1,0:ncpus-1) :: nmig
!
      real :: r0,rn,x0,xn,y0,yn,dxc,dyc,dr,dth,radius,rr
      real :: delr,delp,theta,fr,fp,p1,p2,p3,p4,interp_pot
      real :: distx,disty,fx,fy,delx,dely,xp,yp,k2
      real :: Lxn,Lyn,th0,thn,theta0,theta1
      real :: yps,dxc1,dyc1,dr1,dth1
!
      integer :: nr,nth,ir,im,ix1,ix2,iy1,iy2,nnxgrid,nnygrid
      integer :: i,j,ir1,ir2,ip1,ip2,nnx,nny
      integer :: iproc_ip1,iproc_ip2,iphi,ik
      integer :: iproc_iy1,iproc_iy2,ikx,iky
      integer :: tmp,tmp1,tmp2,nrgrid,nthgrid
!
      logical :: lprint
!
      lprint = lroot.and.(ip<=9)
!
! transform the grid to cartesian
! keep the notation consistent.
!
      rad=x(l1:l2)     ; tht=y(m1:m2)
      nr =nx           ; nth=ny
      r0=rad(1)        ; theta0=xyz0(2)+.5*dth 
      rn=rad(nr)       ; theta1=xyz1(2)-.5*dth
      nrgrid=nxgrid    ; nthgrid=nygrid
      dr=dx            ; dth=dy
      dr1=1./dr        ; dth1=1./dth
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
! Need the serial theta later in order to compute the 
! azimuthal displacement in parallel
!
      do i=1,nygrid
        theta_serial(i)=1.*(i-1)/(nthgrid-1)*(theta1-theta0)+theta0
      enddo
      yserial = xc
!
! Loop throught the cylindrical grid, check where they pertain, and send 
! a package (like a tarball) from processor to processor
! like a tarball
!
      nmig=0.
      do iphi=1,nth
        do ir=1,nr
!
          yp=rad(ir)*sin(tht(iphi))
          iy1 = floor((yp-y0)*dyc1)+1 ; iy2 = iy1+1
          if (iy1 .lt.  1)      call fatal_error("","iy1 lt 1")
          if (iy2 .gt. nnygrid) call fatal_error("","iy2 gt nnygrid")
!
! Which Cartesian processors to these cells belong to?          
!
          iproc_iy1 = (iy1-1)/nny 
          iproc_iy2 = (iy2-1)/nny
!
          if (iproc_iy1 .ge. ncpus) call fatal_error("","iproc_iy1 doesn't exist: up")
          if (iproc_iy2 .ge. ncpus) call fatal_error("","iproc_iy2 doesn't exist: up")
          if (iproc_iy1 .lt. 0    ) call fatal_error("","iproc_iy1 doesn't exist: do")
          if (iproc_iy2 .lt. 0    ) call fatal_error("","iproc_iy2 doesn't exist: do")
!
! Prepare the arrays to be sent
!
          if (iproc_iy1 /= iproc) then
            nmig(iproc,iproc_iy1) = nmig(iproc,iproc_iy1) + 1
            phi_send(iproc_iy1,nmig(iproc,iproc_iy1))=phi(ir,iphi,npoint-nghost)
            sri_send(iproc_iy1,nmig(iproc,iproc_iy1))=ir+(iphi+nth*iproc-1)*nr !serial index
          endif
          !if iproc_iy2/=iproc_iy1, another processor will receive
          if ((iproc_iy2 /= iproc).and.(iproc_iy2 /= iproc_iy1)) then
            nmig(iproc,iproc_iy2) = nmig(iproc,iproc_iy2) + 1
            phi_send(iproc_iy2,nmig(iproc,iproc_iy2))=phi(ir,iphi,npoint-nghost)
            sri_send(iproc_iy2,nmig(iproc,iproc_iy2))=ir+(iphi+nth*iproc-1)*nr
          endif
!
! ugly, but let's see... send all first and last rows of each proc to the neighbour processor
!          
          if (iproc_iy2/=ncpus-1) then !not last row of last processor
            tmp2=iy2-iproc_iy2*nny
            if ((tmp2==nny).and.(iproc_iy2+1/=iproc)) then
              nmig(iproc,iproc_iy2+1) = nmig(iproc,iproc_iy2+1) + 1
              phi_send(iproc_iy2+1,nmig(iproc,iproc_iy2+1))=&
                   phi(ir,iphi,npoint-nghost)
              sri_send(iproc_iy2+1,nmig(iproc,iproc_iy2+1))=&
                   ir+(iphi+nth*iproc-1)*nr
            endif
          endif
!
          if (iproc_iy1/=0) then !not first row of first processor
            tmp1=iy1-iproc_iy1*nny
            if ((tmp1==1).and.(iproc_iy1-1/=iproc)) then
              nmig(iproc,iproc_iy1-1) = nmig(iproc,iproc_iy1-1) + 1
              phi_send(iproc_iy1-1,nmig(iproc,iproc_iy1-1))=&
                   phi(ir,iphi,npoint-nghost)
              sri_send(iproc_iy1-1,nmig(iproc,iproc_iy1-1))=&
                   ir+(iphi+nth*iproc-1)*nr
            endif
          endif
!
        enddo
      enddo
!
!  Share information about size of the migrating chunks
!
      do i=0,ncpus-1
        if (iproc/=i) then
          call mpirecv_int(nmig(i,iproc), 1, i, 111)
        else
          do j=0,ncpus-1
            if (iproc/=j) call mpisend_int(nmig(iproc,j), 1, j, 111)
          enddo
        endif
      enddo
!
!  Do the communication between all processors - no need for broadcasting
!  Each processors knows what the others need.
!     
      do i=0,ncpus-1
        !directed send
        if (iproc==i) then
          !send to all processors
          do j=0,ncpus-1
            if (iproc/=j .and. nmig(iproc,j)/=0) then
              call mpisend_real(phi_send(j,1:nmig(iproc,j)),nmig(iproc,j),j,112)
              call mpisend_int (sri_send(j,1:nmig(iproc,j)),nmig(iproc,j),j,113)
            endif
          enddo
        else 
          !set to receive
          if (nmig(i,iproc)/=0) then
            call mpirecv_real(phi_recv(i,1:nmig(i,iproc)),nmig(i,iproc),i,112)
            call mpirecv_int (sri_recv(i,1:nmig(i,iproc)),nmig(i,iproc),i,113)
          endif
        endif
      enddo
!
! Transform to the Cartesian grid. All the information needed is stored in phi_recv
!
      do m=1,nny
        do i=1,nnx  
!
          radius=sqrt(xc(i)**2+yc(m)**2)
!
          if ((radius.ge.r0).and.(radius.le.rn)) then 
            ir1=floor((radius-r0)*dr1 ) +1;ir2=ir1+1
!
            delr=radius-rad(ir1)
            fr=delr*dr1
!
! this should never happen, but is here for warning
!
            if (ir1.lt.1 ) call fatal_error("","cyl2cart: ir1<1")
            if (ir2.gt.nr) call fatal_error("","cyl2cart: ir2>nr")
!
            theta=atan2(yc(m),xc(i))
            !for getting the crossed one, the index 
            !must be the serial index
            ip1=floor((theta -theta0)*dth1)+1;ip2=ip1+1
            if (ip1==0) then
              ip1=nthgrid
              delp=theta-theta_serial(ip1) + 2*pi
            else
              delp=theta-theta_serial(ip1)
            endif
            if (ip2==nthgrid+1) ip2=1
!
            !processor theta belongs to
            iproc_ip1 = (ip1-1)/nth
            iproc_ip2 = (ip2-1)/nth
!
            if (iproc_ip1 .ge. ncpus) call fatal_error("","iproc_ip1 doesn't exist: up")
            if (iproc_ip2 .ge. ncpus) call fatal_error("","iproc_ip2 doesn't exist: up")
            if (iproc_ip1 .lt. 0    ) call fatal_error("","iproc_ip1 doesn't exist: do")
            if (iproc_ip2 .lt. 0    ) call fatal_error("","iproc_ip2 doesn't exist: do")
!
            if (iproc_ip1 == iproc) then
              p1=phi(ir1,ip1-nth*iproc,npoint-nghost)
              p2=phi(ir2,ip1-nth*iproc,npoint-nghost)
            else
              !loop through the nmig to find ip1
              if (nmig(iproc_ip1,iproc) /=0) then 
                do j=1,nmig(iproc_ip1,iproc)
                  if (sri_recv(iproc_ip1,j) == ir1+(ip1-1)*nr) &
                       p1=phi_recv(iproc_ip1,j)
                  if (sri_recv(iproc_ip1,j) == ir2+(ip1-1)*nr) &
                       p2=phi_recv(iproc_ip1,j)
                enddo
              else
                print*,iproc_ip1,iproc
                call fatal_error("","nmig is empty."//&
                     " Something is wrong. Go check.")
              endif
            endif
!
            if (iproc_ip2 == iproc) then
              p3=phi(ir1,ip2-nth*iproc,npoint-nghost)
              p4=phi(ir2,ip2-nth*iproc,npoint-nghost)
            else
              if (nmig(iproc_ip2,iproc) /=0) then 
                do j=1,nmig(iproc_ip2,iproc)
                  if (sri_recv(iproc_ip2,j) == ir1+(ip2-1)*nr) &
                       p3=phi_recv(iproc_ip2,j)
                  if (sri_recv(iproc_ip2,j) == ir2+(ip2-1)*nr) &
                       p4=phi_recv(iproc_ip2,j)
                enddo
              else
                print*,iproc_ip2,iproc
                call fatal_error("","nmig is empty."//&
                     " Something is wrong. Go check.")
              endif
            endif
!
            !p1=phi(ir1,ip1),p2=phi(ir2,ip1)
            !p3=phi(ir1,ip2),p4=phi(ir2,ip2)
!
            fp=delp*dth1
            interp_pot=fr*fp*(p1-p2-p3+p4) + fr*(p2-p1) + fp*(p3-p1) + p1
!
            nphi(i,m)=interp_pot
          else 
            nphi(i,m)=0.
          endif
        enddo
      enddo
      if (lprint) print*,'finished transforming to Cartesian'
!      
!  Forward transform (to k-space).
!
      nb1=0.
!
      call fourier_transform_xy_xy_other(nphi,nb1)
!
!  Solve Poisson equation
!
      Lxn=2*xc(nnx)
      Lyn=Lxn
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
              call stop_it("3d case not implemented yet")
!
!  Razor-thin approximation. Here we solve the equation
!    del2Phi=4*pi*G*Sigma(x,y)*delta(z)
!  The solution at scale k=(kx,ky) is
!    Phi(x,y,z)=-(2*pi*G/|k|)*Sigma(x,y)*exp[i*(kx*x+ky*y)-|k|*|z|]
!
            else
            
              k2 = (kkx_fft(ikx)**2+kky_fft(iky+ipy*nny)**2)
!              if (ipy==0) print*,kkx_fft(ikx),kky_fft(iky+ipy*nny)
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
!  Transform back from Cartesian to cylindrical
!
      !
      !prepare the grid for directed send
      ! 
      
      nmig=0.
      do m=1,nny
        do i=1,nnx
          radius=sqrt(xc(i)**2+yc(m)**2)
          !if ((radius.ge.r0).and.(radius.le.rn)) then 
          if ((radius.ge.(r0-2*dr)).and.(radius.le.(rn+2*dr))) then !give some room for ir=1 and ir=nr
            theta=atan2(yc(m),xc(i))
            ip1=floor((theta -theta0)*dth1)+1;ip2=ip1+1
            if (ip1==0)         ip1=nthgrid
            if (ip2==nthgrid+1) ip2=1
!
! Which cylindrical processors to these cells belong to?          
!
            iproc_ip1 = (ip1-1)/nth 
            iproc_ip2 = (ip2-1)/nth
!
            if (iproc_ip1 .ge. ncpus) call fatal_error("","iproc_ip1 doesn't exist: up")
            if (iproc_ip2 .ge. ncpus) call fatal_error("","iproc_ip2 doesn't exist: up")
            if (iproc_ip1 .lt. 0    ) call fatal_error("","iproc_ip1 doesn't exist: do")
            if (iproc_ip2 .lt. 0    ) call fatal_error("","iproc_ip2 doesn't exist: do")
!
! Prepare the arrays to be sent
!
            if (iproc_ip1 /= iproc) then
              nmig(iproc,iproc_ip1) = nmig(iproc,iproc_ip1) + 1
              phi_send(iproc_ip1,nmig(iproc,iproc_ip1))=nphi(i,m)
              sri_send(iproc_ip1,nmig(iproc,iproc_ip1))=i+(m+nny*iproc-1)*nnx !serial index
            endif
          !if iproc_iy2/=iproc_iy1, another processor will receive
            if ((iproc_ip2 /= iproc).and.(iproc_ip2 /= iproc_ip1)) then
              nmig(iproc,iproc_ip2) = nmig(iproc,iproc_ip2) + 1
              phi_send(iproc_ip2,nmig(iproc,iproc_ip2))=nphi(i,m)
              sri_send(iproc_ip2,nmig(iproc,iproc_ip2))=i+(m+nny*iproc-1)*nnx
            endif
!
! ugly, but let's see... send all first and last rows of each proc to the neighbour processor
!          
            !last 3 rows of iproc_ip2
            tmp2=ip2-iproc_ip2*nth
            if (tmp2 >= nth-5) then
              if (iproc_ip2/=ncpus-1) then
                tmp=iproc_ip2+1
              else
                tmp=0 !root
              endif
              if (iproc/=tmp) then
                nmig(iproc,tmp) = nmig(iproc,tmp) + 1
                phi_send(tmp,nmig(iproc,tmp))=nphi(i,m)
                sri_send(tmp,nmig(iproc,tmp))=i+(m+nny*iproc-1)*nnx
              endif
            endif
!
            !first rows of iproc_ip1
            tmp1=ip1-iproc_ip1*nth
            if (tmp1<= 6) then
              if (iproc_ip1/=0) then
                tmp=iproc_ip1-1
              else
                tmp=ncpus-1 !last
              endif
              if (iproc/=tmp) then
                nmig(iproc,tmp) = nmig(iproc,tmp) + 1
                phi_send(tmp,nmig(iproc,tmp))=nphi(i,m)
                sri_send(tmp,nmig(iproc,tmp))=i+(m+nny*iproc-1)*nnx
              endif
            endif
!
          endif
        enddo
      enddo
!
!  Share information about size of the migrating chunks
!
      do i=0,ncpus-1
        if (iproc/=i) then
          call mpirecv_int(nmig(i,iproc), 1, i, 222)
        else
          do j=0,ncpus-1
            if (iproc/=j) call mpisend_int(nmig(iproc,j), 1, j, 222)
          enddo
        endif
      enddo
!
!  Do the communication between all processors - no need for broadcasting
!  Each processors knows what the others need.
!     
      do i=0,ncpus-1
        !directed send
        if (iproc==i) then
          !send to all processors
          do j=0,ncpus-1
            if (iproc/=j .and. nmig(iproc,j)/=0) then
              call mpisend_real(phi_send(j,1:nmig(iproc,j)),nmig(iproc,j),j,223)
              call mpisend_int (sri_send(j,1:nmig(iproc,j)),nmig(iproc,j),j,224)
            endif
          enddo
        else 
          !set to receive
          if (nmig(i,iproc)/=0) then
            call mpirecv_real(phi_recv(i,1:nmig(i,iproc)),nmig(i,iproc),i,223)
            call mpirecv_int (sri_recv(i,1:nmig(i,iproc)),nmig(i,iproc),i,224)
          endif
        endif
      enddo
!
!  Convert back to cylindrical. All the information needed is stored in phi_recv
!
      do ip=1,Nth
        do ir=1,Nr
!
          xp=rad(ir)*cos(tht(ip))
          yp=rad(ir)*sin(tht(ip))
!
          distx = xp - x0
          disty = yp - y0
!
          ix1 = floor(distx*dxc1)+1 ; ix2 = ix1+1
          iy1 = floor(disty*dyc1)+1 ; iy2 = iy1+1
!
          if (ix1 .lt.  1)      call fatal_error("","ix1 lt 1")
          if (iy1 .lt.  1)      call fatal_error("","iy1 lt 1")
          if (ix2 .gt. nnxgrid) call fatal_error("","ix2 gt nnxgrid")
          if (iy2 .gt. nnygrid) call fatal_error("","iy2 gt nnygrid")
!
! What processor does iy1 belong to?
! 
          iproc_iy1 = (iy1-1)/nny
          iproc_iy2 = (iy2-1)/nny
!
          if (iproc_iy1 .ge. ncpus) call fatal_error("","iproc_iy1 doesn't exist: up")
          if (iproc_iy2 .ge. ncpus) call fatal_error("","iproc_iy2 doesn't exist: up")
          if (iproc_iy1 .lt. 0    ) call fatal_error("","iproc_iy1 doesn't exist: do")
          if (iproc_iy2 .lt. 0    ) call fatal_error("","iproc_iy2 doesn't exist: do")
!
          if (iproc_iy1 == iproc) then
            p1=nphi(ix1,iy1-nny*iproc)
            p2=nphi(ix2,iy1-nny*iproc)
          else
            !loop through the nmig to find iy1
            if (nmig(iproc_iy1,iproc) /=0) then 
              do j=1,nmig(iproc_iy1,iproc)
                if (sri_recv(iproc_iy1,j) == ix1+(iy1-1)*nnx) &
                     p1=phi_recv(iproc_iy1,j)
                if (sri_recv(iproc_iy1,j) == ix2+(iy1-1)*nnx) &
                     p2=phi_recv(iproc_iy1,j)
              enddo
            else
              print*,iproc_iy1,iproc
              call fatal_error("","nmig iproc_iy1 is empty."//&
                   " Something is wrong. Go check.")
            endif
          endif
!
          if (iproc_iy2 == iproc) then
            p3=nphi(ix1,iy2-nny*iproc)
            p4=nphi(ix2,iy2-nny*iproc)
          else
            if (nmig(iproc_iy2,iproc) /=0) then 
              do j=1,nmig(iproc_iy2,iproc)
                if (sri_recv(iproc_iy2,j) == ix1+(iy2-1)*nnx) &
                     p3=phi_recv(iproc_iy2,j)
                if (sri_recv(iproc_iy2,j) == ix2+(iy2-1)*nnx) &
                     p4=phi_recv(iproc_iy2,j)
              enddo
            else
              print*,iproc_iy2,iproc
              call fatal_error("","nmig iproc_iy2 is empty."//&
                   " Something is wrong. Go check.")
            endif
          endif

          delx=xp-xc(ix1);fx=delx*dxc1
          dely=yp-yserial(iy1);fy=dely*dyc1
!
! Bilinear interpolation
!
          interp_pot=fx*fy*(p1-p2-p3+p4) + fx*(p2-p1) + fy*(p3-p1) + p1
!
          do n=1,nz
            phi(ir,ip,n)=interp_pot
          enddo
!
        enddo
      enddo
!
    endsubroutine inverse_laplacian_fft_cyl_mpi
!***********************************************************************
    subroutine inverse_laplacian_semispectral(phi)
!
!  Solve the Poisson equation by Fourier transforming in the xy-plane and
!  solving the discrete matrix equation in the z-direction.
!
!  19-dec-2006/anders: coded
!
      use General, only: tridag
      use Mpicomm, only: transp_xz, transp_zx
!
      real, dimension (nx,ny,nz) :: phi
!
      real, dimension (nx,ny,nz) :: b1
      real, dimension (nzgrid,nx/nprocz) :: rhst, b1t
      real, dimension (nzgrid) :: a_tri, b_tri, c_tri, r_tri, u_tri
      real :: k2
      integer :: ikx, iky, ikz
      logical :: err
!
!  identify version
!
      if (lroot .and. ip<10) call cvs_id( &
        "$Id: poisson.f90,v 1.42 2008-02-27 12:35:06 wlyra Exp $")
!
!  The right-hand-side of the Poisson equation is purely real.
!
      b1 = 0.0
!  Forward transform (to k-space).
      if (lshear) then
        call fourier_transform_shear_xy(phi,b1)
      else
        call fourier_transform_xy(phi,b1)
      endif
!
!  Solve for discrete z-direction with zero density above and below z-boundary.
!
      do iky=1,ny
        call transp_xz(phi(:,iky,:),rhst)
        a_tri(:)=1.0/dz**2
        c_tri(:)=1.0/dz**2
        do ikx=1,nxgrid/nprocz
          k2=kx_fft(ikx+nz*ipz)**2+ky_fft(iky)**2
          b_tri=-2.0/dz**2-k2
!
          if (k2==0.0) then
            b_tri(1)=-2.0/dz**2-2*dz/xyz0(3)
            c_tri(1)=1.0/dz**2+1.0
            b_tri(nzgrid)=-2.0/dz**2-2*dz/xyz1(3)
            a_tri(nzgrid)=1.0/dz**2+1.0
          else
            b_tri(1)=-2.0/dz**2-2*sqrt(k2)*dz
            c_tri(1)=1.0/dz**2+1.0
            b_tri(nzgrid)=-2.0/dz**2-2*sqrt(k2)*dz
            a_tri(nzgrid)=1.0/dz**2+1.0
          endif
!
          r_tri=rhst(:,ikx)
          call tridag(a_tri,b_tri,c_tri,r_tri,u_tri,err)
          rhst(:,ikx)=u_tri
        enddo
        call transp_zx(rhst,phi(:,iky,:))
      enddo
!
      do iky=1,ny
        call transp_xz(b1(:,iky,:),b1t)
        a_tri(:)=1.0/dz**2
        c_tri(:)=1.0/dz**2
        do ikx=1,nxgrid/nprocz
          k2=kx_fft(ikx+nz*ipz)**2+ky_fft(iky)**2
          b_tri=-2.0/dz**2-k2
!
          if (k2==0.0) then
            b_tri(1)=-2.0/dz**2-2*dz/xyz0(3)
            c_tri(1)=1.0/dz**2+1.0
            b_tri(nzgrid)=-2.0/dz**2-2*dz/xyz1(3)
            a_tri(nzgrid)=1.0/dz**2+1.0
          else
            b_tri(1)=-2.0/dz**2-2*sqrt(k2)*dz
            c_tri(1)=1.0/dz**2+1.0
            b_tri(nzgrid)=-2.0/dz**2-2*sqrt(k2)*dz
            a_tri(nzgrid)=1.0/dz**2+1.0
          endif
!
          r_tri=b1t(:,ikx)
          call tridag(a_tri,b_tri,c_tri,r_tri,u_tri,err)
          b1t(:,ikx)=u_tri
        enddo
        call transp_zx(b1t,b1(:,iky,:))
      enddo
!  Inverse transform (to real space).
      if (lshear) then
        call fourier_transform_shear_xy(phi,b1,linv=.true.)
      else
        call fourier_transform_xy(phi,b1,linv=.true.)
      endif
!
    endsubroutine inverse_laplacian_semispectral
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
