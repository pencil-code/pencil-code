! $Id$
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
!
  use Cdata
  use Cparam
  use Fourier
  use Messages
!
  implicit none
!
  real :: kmax=0.0
  logical :: lrazor_thin=.false., lsemispectral=.false., lklimit_shear=.false.
  logical :: lexpand_grid=.false.
  logical :: lisoz=.false.
!
  include 'poisson.h'
!
  namelist /poisson_init_pars/ &
      lsemispectral, kmax, lrazor_thin, lklimit_shear, lexpand_grid, lisoz
!
  namelist /poisson_run_pars/ &
      lsemispectral, kmax, lrazor_thin, lklimit_shear, lexpand_grid, lisoz
!
  logical :: luse_fourier_transform = .false.
!
  contains
!***********************************************************************
    subroutine initialize_poisson
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  18-oct-07/anders: adapted
!
      if (lrazor_thin) then
        if (nzgrid/=1) then
          if (lroot) print*, 'inverse_laplacian_fft: razor-thin approximation only works with nzgrid==1'
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
      if (lklimit_shear) kmax = kx_nyq/sqrt(2.0)
!
!  Dimensionality
!
      call decide_fourier_routine
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
      real, dimension(nx,ny,nz), intent(inout) :: phi
!
      if (lcylindrical_coords) then
        if (lroot) print*,'You are using cylindrical coordinates. '//&
             'Use poisson_multigrid.f90 instead'
        call fatal_error("inverse_laplacian","")
      endif
!
      if (lspherical_coords) then
        if (lroot) then
          print*,'There is no poisson solver for spherical '
          print*,'coordinates yet. Please feel free to implement it. '
          print*,'Many people will thank you for that.'
          call fatal_error("inverse_laplacian","")
        endif
      endif
!
      if (lsemispectral) then
        call inverse_laplacian_semispectral(phi)
      else if (lisoz) then
        call inverse_laplacian_isoz(phi)
      else if (lexpand_grid) then
        if (lroot) then
          print*,'The poisson solver for global disk was moved to an '
          print*,'experimental module because it uses too much memory. '
          print*,'If your memory requirements are modest ( < 500 proc) '
          print*,'use POISSON=experimental/poisson_expand.f90 in your '
          print*,'Makefile.local file'
          call fatal_error("inverse_laplacian","")
        endif
        !call inverse_laplacian_expandgrid(phi)
      else if (nprocx==1.and.(nxgrid/=1.and.nygrid/=1.and.nzgrid/=1.and.nxgrid/=nzgrid)) then
        call inverse_laplacian_fft_z(phi)
      else
        call inverse_laplacian_fft(phi)
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
      if (lroot .and. ip<10) call svn_id( &
        "$Id$")
!
!  The right-hand-side of the Poisson equation is purely real.
!
      b1 = 0.0
!
!  Forward transform (to k-space).
!
      if (luse_fourier_transform) then
        if (lshear) then
          call fourier_transform_shear(phi,b1)
        else
          call fourier_transform(phi,b1)
        endif
      else
        if (.not.(nxgrid/=1.and.nzgrid/=1.and.nxgrid/=nzgrid)) then
          call fft_xyz_parallel(phi,b1)
        else
          call fatal_error("inverse_laplacian_fft",&
               "not yet coded for the general case nxgrid/=nzgrid")
        endif
      endif
!
!  Solve Poisson equation.
!
      do ikz=1,nz; do iky=1,ny; do ikx=1,nx
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
                    ky_fft(iky+ipy*ny)**2 + kz_fft(ikx+ipx*nx)**2
            else                ! Order (kx,ky',kz)
              k2 = (kx_fft(ikx+ipx*nx)+deltay/Lx*ky_fft(iky+ipy*ny))**2 + &
                    ky_fft(iky+ipy*ny)**2 + kz_fft(ikz+ipz*nz)**2
            endif
!  The ordering of the array is not important here, because there is no shear!
          else
            k2 = kx_fft(ikx+ipx*nx)**2 + ky_fft(iky+ipy*ny)**2 + kz_fft(ikz+ipz*nz)**2
          endif
!
!  Solution of Poisson equation.
!
          if (k2 > 0.0) then
            phi(ikx,iky,ikz) = -phi(ikx,iky,ikz) / k2
            b1(ikx,iky,ikz)  = - b1(ikx,iky,ikz) / k2
          else
            phi(ikx,iky,ikz) = 0.0
            b1(ikx,iky,ikz) = 0.0
          endif
!
        else
!
!  Razor-thin approximation. Here we solve the equation
!    del2Phi=4*pi*G*Sigma(x,y)*delta(z)
!  The solution at scale k=(kx,ky) is
!    Phi(x,y,z)=-(2*pi*G/|k|)*Sigma(x,y)*exp[i*(kx*x+ky*y)-|k|*|z|]
!
          if (lshear) then
            k2 = (kx_fft(ikx+ipx*nx)+deltay/Lx*ky_fft(iky+ipy*ny))**2+ky_fft(iky+ipy*ny)**2
          else
            k2 = kx_fft(ikx+ipx*nx)**2+ky_fft(iky+ipy*ny)**2
          endif
!
          if (k2 > 0.0) then
            phi(ikx,iky,ikz) = -.5*phi(ikx,iky,ikz) / sqrt(k2)
            b1(ikx,iky,ikz)  = -.5* b1(ikx,iky,ikz) / sqrt(k2)
          else
            phi(ikx,iky,ikz) = 0.0
            b1(ikx,iky,ikz) = 0.0
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
      if (luse_fourier_transform) then
        if (lshear) then
          call fourier_transform_shear(phi,b1,linv=.true.)
        else
          call fourier_transform(phi,b1,linv=.true.)
        endif
      else
        call fft_xyz_parallel(phi,b1,linv=.true.)
      endif
!
    endsubroutine inverse_laplacian_fft
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
      integer :: ikx, iky
      logical :: err
!
!  Identify version.
!
      if (lroot .and. ip<10) call svn_id( &
          '$Id$')
!
!  The right-hand-side of the Poisson equation is purely real.
!
      b1 = 0.0
!
!  Forward transform (to k-space).
!
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
!
!  Inverse transform (to real space).
!
      if (lshear) then
        call fourier_transform_shear_xy(phi,b1,linv=.true.)
      else
        call fourier_transform_xy(phi,b1,linv=.true.)
      endif
!
    endsubroutine inverse_laplacian_semispectral
!***********************************************************************
    subroutine inverse_laplacian_fft_z(phi)
!
!  Solve the Poisson equation with nxgrid = nygrid /= nzgrid.
!
!  10-sep-2009/ccyang: coded
!
      use Mpicomm, only: transp_xz, transp_zx
!
      real, dimension(nx,ny,nz), intent(inout) :: phi
!
      real, dimension(nx,ny,nz) :: b1
!
      integer, parameter :: nxt = nx / nprocz
      real, dimension(nzgrid,nxt) :: phirt, b1t
!
      logical :: l1st = .true.
      real, save, dimension(4*nzgrid+15) :: wsave
      real, save, dimension(nzgrid)      :: kz2
      integer, save :: ikx0, iky0
!
      complex, dimension(nzgrid) :: cz
!
      integer :: ix, iy
      real    :: kx2, ky2, a0, a1
!
!  Initialize the array wsave and other constants for future use.
!
      if (l1st) then
        call cffti(nzgrid, wsave)
        kz2 = kz_fft**2
        ikx0 = ipz * nxt
        iky0 = ipy * ny
        l1st = .false.
      endif
!
!  Forward transform in xy
!
      b1 = 0.
      if (nprocx==1) then
        if (nxgrid==nygrid) then
          if (lshear) then
            call fourier_transform_shear_xy(phi, b1)
          else
            call fourier_transform_xy(phi, b1)
            !#ccyang: Note that x and y are swapped in fourier_transform_xy.
          endif
        else
          call fft_xy_parallel(phi,b1)
        endif
      else
        call fatal_error("inverse_laplacian_fft_z",&
             "Does not work for x-parallelization")
      endif
!
!  Convolution in z
!
      if (lshear) a0 = deltay / Lx
      do iy = 1, ny
!
        call transp_xz(phi(:,iy,:),  phirt)
        call transp_xz(b1(:,iy,:), b1t)
!
        if (lshear) then
          ky2 = ky_fft(iky0+iy)**2
          a1  = a0 * ky_fft(iky0+iy)
        else
          ky2 = kx_fft(iky0+iy)**2
        endif
!
        do ix = 1, nxt
!
          if (lshear) then
            kx2 = (kx_fft(ikx0+ix) + a1)**2
          else
            kx2 = ky_fft(ikx0+ix)**2
          endif
!
          cz = cmplx(phirt(:,ix), b1t(:,ix))
          call cfftf (nzgrid, cz, wsave)
          where (kx2 /= 0. .or. ky2 /= 0. .or. kz2 /= 0)
            cz = -cz / (kx2 + ky2 + kz2) / nzgrid
          elsewhere
            cz = 0.
          endwhere
          call cfftb (nzgrid, cz, wsave)
!
          phirt(:,ix) = real(cz)
          b1t(:,ix) = aimag(cz)
!
        enddo
!
        call transp_zx(phirt, phi(:,iy,:))
        call transp_zx(b1t, b1(:,iy,:))
!
      enddo
!
!  Inverse transform in xy
!
      if (nprocx==1) then
        if (nxgrid==nygrid) then
          if (lshear) then
            call fourier_transform_shear_xy(phi,b1,linv=.true.)
          else
            call fourier_transform_xy(phi,b1,linv=.true.)
            !#ccyang: Note that x and y are swapped in fourier_transform_xy.
          endif
        else
          call fft_xy_parallel(phi,b1,linv=.true.)
        endif
      else
        call fatal_error("inverse_laplacian_fft_z",&
             "Does not work for x-parallelization")
      endif
!
      return
!
    endsubroutine inverse_laplacian_fft_z
!***********************************************************************
    subroutine inverse_laplacian_isoz(phi)
!
!  Solve the Poisson equation with isolated boundary condition in the z-
!  direction and periodic in the x- and y- directions.
!
!  11-jan-2008/ccyang: coded
!
!  Reference: Yang, Mac Low, & Menou 2011 (arXiv:1103.3268)
!
      use Mpicomm, only: transp_xz, transp_zx
!
      real, dimension(nx,ny,nz) :: phi
!
      real, dimension(nx,ny,nz) :: b1
!
      integer, parameter :: nzg2 = 2 * nzgrid
      integer, parameter :: nxt = nx / nprocz
      real, dimension(nzgrid,nxt) :: phirt, b1t
!
      logical, save :: l1st = .true.
      real, save, dimension(4*nzg2+15) :: wsave
      real, save, dimension(nzg2)      :: kz2
      real,    save :: dz2h
      integer, save :: ikx0, iky0
!
      complex, dimension(nzg2) :: cz
!
      integer :: ix, iy, iz, iz1
      real    :: kx2, ky2, a0, a1
!
!  Initialize the array wsave and other constants for future use.
!
      if (l1st) then
        call cffti(nzg2, wsave)
        kz2 = (/(iz**2, iz = 0, nzgrid), &
                ((iz - nzgrid)**2, iz = 1, nzgrid - 1)/) * (pi / Lz)**2
        dz2h = .5 * dz * dz
        ikx0 = ipz * nxt
        iky0 = ipy * ny
        l1st = .false.
      endif
!
!  Forward transform in xy
!
      b1 = 0.
      if (lshear) then
        call fourier_transform_shear_xy(phi, b1)
      else
        call fourier_transform_xy(phi, b1)
!#ccyang: Note that x and y are swapped in fourier_transform_xy.
      endif
!
!  Convolution in z
!
      if (lshear) a0 = deltay / Lx
      do iy = 1, ny
        call transp_xz(phi(:,iy,:),  phirt)
        call transp_xz(b1(:,iy,:), b1t)
        if (lshear) then
          ky2 = ky_fft(iky0+iy)**2
          a1  = a0 * ky_fft(iky0+iy)
        else
          ky2 = kx_fft(iky0+iy)**2
        endif
        do ix = 1, nxt
          if (lshear) then
            kx2 = (kx_fft(ikx0+ix) + a1)**2
          else
            kx2 = ky_fft(ikx0+ix)**2
          endif
          if (kx2 /= 0. .or. ky2 /= 0.) then
            cz(1:nzgrid) = (/cmplx(phirt(:,ix), b1t(:,ix))/)
            cz(nzgrid+1:nzg2) = 0.
            call cfftf (nzg2, cz, wsave)
            cz = -cz / (kx2 + ky2 + kz2) / nzg2
            call cfftb (nzg2, cz, wsave)
          else
            cz = 0.
            do iz = 1, nzgrid; do iz1 = 1, nzgrid
              cz(iz) = cz(iz) + cmplx(phirt(iz1,ix), b1t(iz1,ix)) * &
                                abs(iz - iz1) * dz2h
            enddo; enddo
          endif
          phirt(:,ix) = real(cz(1:nzgrid))
          b1t(:,ix) = aimag(cz(1:nzgrid))
        enddo
        call transp_zx(phirt, phi(:,iy,:))
        call transp_zx(b1t, b1(:,iy,:))
      enddo
!
!  Inverse transform in xy
!
      if (lshear) then
        call fourier_transform_shear_xy(phi,b1,linv=.true.)
      else
        call fourier_transform_xy(phi,b1,linv=.true.)
      endif
!
      return
!
    endsubroutine inverse_laplacian_isoz
!***********************************************************************
    subroutine inverse_laplacian_z_2nd_neumann(f)
!
!  19-mar-2018/MR: coded
!  Second-order version in the vertical direction that uses tridag_neumann.
!  On input: phi=div(U), on output: phi=potential of the irrotational flow.
!
      use Fourier, only: fourier_transform_xy, kx_fft2, ky_fft2
      use Mpicomm, only: transp_xz, transp_zx
      !!!use Sub, only: tridag_neumann
!
      real, dimension(:,:,:,:), intent(INOUT)   :: f

      real, dimension(nx,ny,0:nz+1) :: phi, b1
      real, dimension(nx,ny), save :: k2dz2
      real, dimension(nx,ny,1) :: uzhatlow_re, uzhatlow_im

      real, dimension (nz-3) :: a, b
      real, dimension (nz-2) :: c
!
      integer, parameter :: nxt = nx / nprocz
      logical, save :: l1st = .true.
      real,    save :: dz2
      integer, save :: ikx0, iky0
!
      complex, dimension(nx,ny,nz) :: rhs
!
      integer :: ikx, iky, i1
      real    :: ky2
!
      call fatal_error('inverse_laplacian_z_2nd_neumann',"don't call, not yet functioning")
!
!  Initialize the array k2dz2 and other constants for future use.
!
      if (l1st) then
        dz2 = dz**2
        ikx0 = ipz * nxt
        iky0 = ipy * ny
        do iky = 1, ny
          ky2 = ky_fft2(iky0+iky)
          do ikx = 1, nx
            k2dz2(ikx,iky) = (kx_fft2(ikx0+ikx)+ky2)*dz2
          enddo
        enddo
        l1st = .false.
      endif
!
!  Forward transform of div(U) and U on lower boundary in xy.
!
      phi(:,:,1:nz)=f(l1:l2,m1:m2,n1:n2,iphiuu); b1 = 0.
      call fourier_transform_xy(phi(:,:,1:nz), b1(:,:,1:nz))    ! Fourier transform of div(U) in cmplx(phi,b1)
print*, 'nach fourier(div)'
      uzhatlow_re=f(l1:l2,m1:m2,n1,iuz:iuz); uzhatlow_im=0
      call fourier_transform_xy(uzhatlow_re, uzhatlow_im)       ! Fourier transform of U_z on lower boundary
 
      if (lfirst_proc_z) then
        rhs(:,:,1:1)=2.*dz*cmplx(uzhatlow_re,uzhatlow_im)         ! auxiliary to compute phi on boundary
        rhs(:,:,2)=3.*dz2*cmplx(phi(:,:,2),b1(:,:,2))+rhs(:,:,1)
        i1=3
      else
        i1=1
      endif
      rhs(:,:,i1:nz)=dz2*cmplx(phi(:,:,i1:nz),b1(:,:,i1:nz))
!
      !do iky = 1, ny
        !call transp_xz(phi(:,iky,:),  phit)
        !call transp_xz(b1(:,iky,:), b1t)
!
!  Solution of z dependent problem.
!
      a=1.; b=-2.; c=1.; c(1)=2.
      !!!call tridag_neumann(a,b,c,k2dz2,rhs,phi,b1)
        !call transp_zx(phit, phi(:,iky,:))
        !call transp_zx(b1t, b1(:,iky,:))
!
!  Inverse transform in xy
!
      call fourier_transform_xy(phi(:,:,1:nz),b1(:,:,1:nz),linv=.true.)   ! potential in physical space in phi
      f(l1:l2,m1:m2,n1:n2,iphiuu)=phi(:,:,1:nz)
!
    endsubroutine inverse_laplacian_z_2nd_neumann
!***********************************************************************
    subroutine decide_fourier_routine
!
! Decide, based on the dimensionality and on the geometry
! of the grid, which fourier transform routine is to be
! used. "fourier_transform" and "fourier_tranform_shear"
! are functional only without x-parallelization, and for
! cubic, 2D square, and 1D domains only. Everything else
! should use fft_parallel instead.
!
! 05-dec-2011/wlad: coded
!
      logical :: lxpresent,lypresent,lzpresent
      logical :: l3D,l2Dxy,l2Dyz,l2Dxz
      logical :: l1Dx,l1Dy,l1Dz
!
! Check the dimensionality and store it in logicals
!
      lxpresent=(nxgrid/=1)
      lypresent=(nygrid/=1)
      lzpresent=(nzgrid/=1)
!
      l3D  =     lxpresent.and.     lypresent.and.     lzpresent
      l2Dxy=     lxpresent.and.     lypresent.and..not.lzpresent
      l2Dyz=.not.lxpresent.and.     lypresent.and.     lzpresent
      l2Dxz=     lxpresent.and..not.lypresent.and.     lzpresent
      l1Dx =     lxpresent.and..not.lypresent.and..not.lzpresent
      l1Dy =.not.lxpresent.and.     lypresent.and..not.lzpresent
      l1Dz =.not.lxpresent.and..not.lypresent.and.     lzpresent
!
      if (ldebug.and.lroot) then
        if (l3D)   print*,"This is a 3D simulation"
        if (l2Dxy) print*,"This is a 2D xy simulation"
        if (l2Dyz) print*,"This is a 2D yz simulation"
        if (l2Dxz) print*,"This is a 2D xz simulation"
        if (l1Dx)  print*,"This is a 1D x simulation"
        if (l1Dy)  print*,"This is a 1D y simulation"
        if (l1Dz)  print*,"This is a 1D z simulation"
      endif
!
! The subroutine "fourier_transform" should only be used
! for 1D, square 2D or cubic domains without x-parallelization.
! Everything else uses fft_parallel.
!
      luse_fourier_transform=(nprocx==1.and.&
           (l1dx                                       .or.&
            l1dy                                       .or.&
            l1dz                                       .or.&
            (l2dxy.and.nxgrid==nygrid)                 .or.&
            (l2dxz.and.nxgrid==nzgrid)                 .or.&
            (l2dyz.and.nygrid==nzgrid)                 .or.&
            (l3d.and.nxgrid==nygrid.and.nygrid==nzgrid)))
!
    endsubroutine decide_fourier_routine
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
