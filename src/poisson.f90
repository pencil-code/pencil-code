! $Id: poisson.f90,v 1.35 2007-12-03 21:19:41 wlyra Exp $

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
        if (lcylindrical_coords) then
          call inverse_laplacian_semispec_cyl(phi)
        else
          call inverse_laplacian_semispectral(phi)
        endif
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
        "$Id: poisson.f90,v 1.35 2007-12-03 21:19:41 wlyra Exp $")
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
      use Mpicomm,only:stop_it
!
      real, dimension (nx,ny,nz) :: phi
      real, dimension (nx,ny) :: sigmaxy,phixy
!     
      real, dimension (2*nx,2*ny) :: nb1,nsigmaxy,nphi
!
      real, dimension(nx) :: rad,tht,xc,yc
      real, dimension(2*nx) :: xf,yf,rr,kkx_fft,kky_fft
      real :: r0,rn,x0,xn,y0,yn,dxc,dyc,dr,dth,radius
      real :: delr,delp,theta,fr,fp,p1,p2,p3,p4,interp_pot
      integer :: nr,nth,ir,im,ix1,ix2,iy1,iy2
      integer :: i,ir1,ir2,ip1,ip2,nnx,nny

      real :: distx,disty,fx,fy,delx,dely,xp,yp,k2
      real :: Lxn,Lyn,th0,thn
!
      integer :: ikx, iky,ng
!
! transform the grid to cartesian 
!
      !keep the notation consistent. also use nx=nx,ny=ny
      nr=nx;rad=x(l1:l2) 
      nth=ny;tht=y(m1:m2)
!
      rn=rad(nr);r0=rad(1)
      thn=tht(nth);th0=tht(1)
      xn=rad(nr);x0=-xn
      yn=rad(nr);y0=-yn
      dr =dx
      dth=dy
!
      do i=1,nx
        xc(i)=1.*(i-1)/(nx-1)*(xn-x0)+x0
        yc(i)=xc(i)
      enddo
      dxc=xc(2)-xc(1)
      dyc=dxc
!
      do m=1,ny
        im=m1-1+m
        do i=1,nx
          radius=sqrt(xc(i)**2+yc(m)**2)
          if ((radius.ge.r0).and.(radius.le.rn)) then 
!       
            theta=atan2(yc(m),xc(i))
!              
            ir1=floor((radius-r0)/dr ) +1;ir2=ir1+1
            ip1=floor((theta -th0)/dth)+1;ip2=ip1+1
!
            if (ip1==0) then
              ip1=nth
              delp=theta-tht(ip1) + 2*pi
            else
              delp=theta -tht(ip1)
            endif
            if (ip2==nth+1) ip2=1
!
            if (ir1.lt.1     ) call stop_it("cyl2cart: ir1<1")
            if (ir2.gt.nr    ) call stop_it("cyl2cart: ir2>nr")
            if (ip1.lt.0     ) call stop_it("cyl2cart: ip1<0")
            if (ip2.gt.nth+1 ) call stop_it("cyl2cart: ip2>nth+1")
!
            delr=radius-rad(ir1)

            fr=delr/dr
            fp=delp/dth
!              
            p1=phi(ir2,ip1,npoint-nghost)
            p2=phi(ir1,ip1,npoint-nghost)
            p3=phi(ir2,ip2,npoint-nghost)
            p4=phi(ir1,ip2,npoint-nghost)
!
            interp_pot=fr*fp*(p3+p2-p1-p4) + fr*(p1-p2) + fp*(p4-p2) + p2
!
            sigmaxy(i,m)=interp_pot
          else 
            sigmaxy(i,m)=0.
          endif
        enddo
      enddo
      
      do n=1,nz
        phi(:,:,n)=sigmaxy
      enddo
!
!  Now expand the grid and prepare the fourier transform
!  The right-hand-side of the Poisson equation is purely real.
!
      nb1=0.0
!
!  Expand the grid
!
      nny=nint(2.*ny)
      nnx=nint(2.*nx)
      ng=nint(nx/2.)
      do i=1,ng
        xf(i)      =xc(1) -(ng+1-i)*dxc
        xf(nnx+1-i)=xc(nx)+(ng+1-i)*dxc
      enddo
      do i=1,nx
        xf(ng+i)=xc(i)
      enddo
      yf=xf
!
      do m=1,nny
        rr=sqrt(xf**2+yf(m)**2)
        do i=1,nnx
          if (rr(i) .gt. rn) then 
!zero mass reservoir
            nsigmaxy(i,m)=0.
          else
            !copy original density
!              
            ix = i-nx/2
            iy = m-ny/2
!
            nsigmaxy(i,m)=sigmaxy(ix,iy)
!
          endif
        enddo
      enddo
!
!  Forward transform (to k-space).
!
      nb1=0.
      call fourier_transform_other(nsigmaxy,nb1)
      nphi=nsigmaxy
!
!  Solve Poisson equation
!
      Lxn=2*xf(nnx)
      Lyn=2*yf(nny)
!
      kkx_fft=cshift((/(i-(nnx+1)/2,i=0,nnx-1)/),+(nnx+1)/2)*2*pi/Lxn
      kky_fft=cshift((/(i-(nny+1)/2,i=0,nny-1)/),+(nny+1)/2)*2*pi/Lyn
!
      do iky=1,nny
        do ikx=1,nnx
          if ((kkx_fft(ikx)==0.0) .and. (kky_fft(iky)==0.0)) then
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
            
            k2 = (kkx_fft(ikx)**2+kky_fft(iky)**2)
            nphi(ikx,iky) = -.5*nphi(ikx,iky) / sqrt(k2)
            nb1(ikx,iky) = -.5*nb1(ikx,iky) / sqrt(k2)
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
      call fourier_transform_other(nphi,nb1,linv=.true.)
!
!  Map newphi back into phi
!
      do m=1,ny
        do i=1,nx
!
          ix=i+nx/2
          iy=m+ny/2
!
          phixy(i,m)=nphi(ix,iy)
!
        enddo
      enddo
!
!  Convert back to cylindrical 
!
      do n=1,Nz;do ir=1,Nr; do ip=1,Nth
    
        xp=rad(ir)*cos(tht(ip))
        yp=rad(ir)*sin(tht(ip))
          
        distx = xp - x0
        disty = yp - y0

        ix1 = floor(distx/dxc)+1 ; ix2 = ix1+1
        iy1 = floor(disty/dyc)+1 ; iy2 = iy1+1
!
        if (ix1 .lt.  1)  call stop_it("ix1 lt 1")
        if (iy1 .lt.  1)  call stop_it("iy1 lt 1")
        if (ix2 .gt. Nx)  call stop_it("ix2 gt Nx")
        if (iy2 .gt. Ny)  call stop_it("iy2 gy Ny")
!          
        delx = xp - xc(ix1)
        dely = yp - yc(iy1) 
!
        fx=delx/dxc
        fy=dely/dyc
!         
        p1=phixy(ix2,iy1)
        p2=phixy(ix1,iy1)
        p3=phixy(ix2,iy2)
        p4=phixy(ix1,iy2)
!
        interp_pot=fx*fy*(p3+p2-p1-p4) + fx*(p1-p2) + fy*(p4-p2) + p2
!          
        phi(ir,ip,1:nz)=interp_pot
!
      enddo;enddo;enddo
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
        "$Id: poisson.f90,v 1.35 2007-12-03 21:19:41 wlyra Exp $")
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
    subroutine inverse_laplacian_semispec_cyl(phi)
!
! Solve the cylindrical Poisson equation by Fourier-transforming 
! the azimuthal direction and solving the tridiagonal system in 
! the radial direction.
!
! Currently works only for 2D cases
!
! 16-nov-07/wlad: coded
!
      use General, only: tridag
      use Mpicomm, only: transp_xy
!
      real, dimension (nx,ny,nz) :: phi,b1
!      
      real, dimension (nygrid,nx/nprocy) :: tmp
      real, dimension (nygrid,nx/nprocy,nz) :: phit,b1t
      real, dimension (nx) :: a_tri,b_tri,c_tri,rad
      real, dimension (nx) :: re_tri,im_tri,u_re_tri,u_im_tri
      real    :: alpha,dr,r0,rn
      integer :: i,m,n,ikx,iky,ikz
      logical :: err
!
      intent(inout) :: phi
!
! From x to r, for the sake of clarity
!
      dr=dx;rad=x(l1:l2);r0=xyz0(1);rn=xyz1(1)
!
! The right-hand-side of the Poisson equation is purely real
!
      b1=0.0
!
      do n=1,nz
!
! Transpose prior to Fourier transform 
!
        tmp=phi(:,:,n);call transp_xy(tmp);phit(:,:,n)=tmp
!
! Fourier transform x (transposed y) to k-space
!
        call fourier_transform_x(phit,b1t)
!
! Transpose back to solve the tridiagonal matrix
!
        tmp=phit(:,:,n);call transp_xy(tmp);phi(:,:,n)=tmp
        tmp= b1t(:,:,n);call transp_xy(tmp); b1(:,:,n)=tmp
!
        do iky=1,ny
!
          do i=2,nx-1
            alpha= .5/((i-1)+r0/dr)
            a_tri(i) = (1 - alpha)/dr**2
            b_tri(i) =-2/dr**2 - (ky_fft(iky)/rad(i))**2
            c_tri(i) = (1 + alpha)/dr**2
          enddo
!
! Symmetric boundary condition using L'Hospital rule
!  lim  f'/f  =  lim  f"/f'
! r-->0         r-->0  
!
          b_tri(1) =-4./dr**2 !- (ky_fft(iky)/rad(1))**2
          c_tri(1) = 4./dr**2
!
! Vanishing second order derivative in the outer boundary
!
          !b_tri(nx)=-2./dr**2 - 2.*dr/rn
          !a_tri(nx)= 1./dr    + 1.
          c_tri(nx)=1/dr**2
          b_tri(nx)=-2/dr**2
          a_tri(nx)=1/dr**2
!
          re_tri = phi(:,iky,n)
          im_tri =  b1(:,iky,n)
!
          re_tri(nx)=0.
          im_tri(nx)=0.
!
          call tridag(a_tri,b_tri,c_tri,re_tri,u_re_tri,err)
          call tridag(a_tri,b_tri,c_tri,im_tri,u_im_tri,err)
!      
          phi(:,iky,n)=u_re_tri
          b1 (:,iky,n)=u_im_tri
!
        enddo
!
! Transpose to perform the fourier transform back to real space
!
        tmp=phi(:,:,n);call transp_xy(tmp);phit(:,:,n)=tmp
        tmp= b1(:,:,n);call transp_xy(tmp); b1t(:,:,n)=tmp
!
! Transform it back to real space
!
        call fourier_transform_x(phit,b1t,linv=.true.)
!
! Transpose the result
!
        tmp=phit(:,:,n);call transp_xy(tmp);phi(:,:,n)=tmp
!
      enddo
!
    endsubroutine inverse_laplacian_semispec_cyl
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
