! $Id: poisson.f90,v 1.22 2007-03-28 18:10:27 dobler Exp $

!
!  This module solves the Poisson equation
!    (d^2/dx^2 + d^2/dy^2 + d^2/dz^2 - c) f = RHS(x,y,z)
!  for the function f(x,y,z).
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

  include 'poisson.h'

  contains

!***********************************************************************
    subroutine inverse_laplacian_fft(phi,kmax,c)
!
!  Solve the Poisson equation by Fourier transforming on a periodic grid.
!  This method works both with and without shear.
!
!  15-may-2006/anders+jeff: coded
!
      real, dimension (nx,ny,nz) :: phi
      real, optional             :: kmax,c
!
      real, dimension (nx,ny,nz) :: b1
      real :: k2
      integer :: ikx, iky, ikz
!
!  identify version
!
      if (lroot .and. ip<10) call cvs_id( &
        "$Id: poisson.f90,v 1.22 2007-03-28 18:10:27 dobler Exp $")
!
!  The right-hand-side of the Poisson equation is purely real.
!
      b1 = 0.0
!  Forward transform (to k-space).
      if (lshear) then
        call fourier_transform_shear(phi,b1)
      else
        call fourier_transform(phi,b1)
      endif
!
!  FT(PHI) = -rhs_poisson_const*FT(rho)/k^2
!
      do ikz=1,nz; do iky=1,ny; do ikx=1,nx
        if ((kx_fft(ikx)==0.0) .and. &
            (ky_fft(iky)==0.0) .and. (kz_fft(ikz)==0.0) ) then
          phi(ikx,iky,ikz) = 0.0
          b1(ikx,iky,ikz) = 0.0
        else
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
          if (present(c)) k2 = k2+c
!
!  Solution of Poisson equation.
!
          phi(ikx,iky,ikz) = -phi(ikx,iky,ikz) / k2
          b1(ikx,iky,ikz) = -b1(ikx,iky,ikz) / k2
!
!  Limit |k| < kmax
!
          if (present(kmax)) then
            if (sqrt(k2)>=kmax) then
              phi(ikx,iky,ikz) = 0.
              b1(ikx,iky,ikz) = 0.
            endif
          endif
!
        endif
      enddo; enddo; enddo
!  Inverse transform (to real space).
      if (lshear) then
        call fourier_transform_shear(phi,b1,linv=.true.)
      else
        call fourier_transform(phi,b1,linv=.true.)
      endif
!
    endsubroutine inverse_laplacian_fft
!***********************************************************************
    subroutine inverse_laplacian_semispectral(phi,c)
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
      real, optional             :: c
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
        "$Id: poisson.f90,v 1.22 2007-03-28 18:10:27 dobler Exp $")
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
          if (present(c)) k2 = k2 + c
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
          if (present(c)) k2 = k2 + c
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

endmodule Poisson

