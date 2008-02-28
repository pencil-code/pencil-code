! $Id: poisson.f90,v 1.44 2008-02-28 11:03:53 wlyra Exp $

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

  real :: kmax=0.0
  logical :: lrazor_thin=.false., lsemispectral=.false., lklimit_shear=.false.

  include 'poisson.h'

  namelist /poisson_init_pars/ &
      lsemispectral, kmax, lrazor_thin, lklimit_shear

  namelist /poisson_run_pars/ &
      lsemispectral, kmax, lrazor_thin, lklimit_shear

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
          call inverse_laplacian_fft_cyl(phi)
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
        "$Id: poisson.f90,v 1.44 2008-02-28 11:03:53 wlyra Exp $")
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
!  Fourier transforming. 
!
!  This method works only for 2D. The frequent broadcast of 
!  a big array gave problems at the PIA cluster in Heidelberg
!  after some thousands of time-steps. The problem was probably
!  due to memory leaking. No problem occured at the UPPMAX
!  cluster in Uppsala. So beware that using this broadcasting
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
      real, dimension(nx)     :: rad
      real, dimension(ny)     :: tht
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
      real    :: r0,rn,theta0,theta1,dr,dth,dr1,dth1
      real    :: theta,radius,k2,xp,yp
      real    :: delr,delp,fr,fp,delx,dely,fx,fy
      real    :: p1,p2,p3,p4,interp_pot
!
      integer :: ix1,ix2,iy1,iy2,ir1,ir2,ip1,ip2
      integer :: i,j,ikx,iky,ir,im,ido,iup
      integer :: nr,nth,nnx,nny,nnghost
      integer :: nrgrid,nthgrid,nnxgrid,nnygrid
!
! Keep the notation consistent.
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
        crossbig(:,1:ny)=nphi(:,1:ny)
      endif
!
!  Convert back to cylindrical
!
      yserial=xc
      do ip=1,Nth
        do ir=1,Nr
!
          xp=rad(ir)*cos(tht(ip))
          yp=rad(ir)*sin(tht(ip))
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
            phi(ir,ip,n)=interp_pot
          enddo
!
        enddo
      enddo
!
    endsubroutine inverse_laplacian_fft_cyl
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
        "$Id: poisson.f90,v 1.44 2008-02-28 11:03:53 wlyra Exp $")
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
