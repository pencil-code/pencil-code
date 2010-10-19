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
!
  include 'poisson.h'
!
  namelist /poisson_init_pars/ &
      lsemispectral, kmax, lrazor_thin, lklimit_shear, lexpand_grid
!
  namelist /poisson_run_pars/ &
      lsemispectral, kmax, lrazor_thin, lklimit_shear, lexpand_grid
!
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
      if (nprocx/=1.and.nzgrid/=1) then
        if (lroot) print*, 'inverse_laplacian_fft: fourier transforms in x-parallel are still restricted to the 2D case.'
        call fatal_error('inverse_laplacian_fft','')
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
    subroutine inverse_laplacian(f,phi)
!
!  Dispatch solving the Poisson equation to inverse_laplacian_fft
!  or inverse_laplacian_semispectral, based on the boundary conditions
!
!  17-jul-2007/wolf: coded wrapper
!
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny,nz) :: phi
!
      intent(inout) :: phi
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
      else
        if (lexpand_grid) then
          if (lroot) then 
            print*,'The poisson solver for global disk was moved to an '
            print*,'experimental module because it uses too much memory. '
            print*,'If your memory requirements are modest ( < 500 proc) '
            print*,'use POISSON=experimental/poisson_expand.f90 in your '
            print*,'Makefile.local file'
            call fatal_error("inverse_laplacian","")
          endif
          !call inverse_laplacian_expandgrid(phi)
        else
          call inverse_laplacian_fft(phi)
        endif
      endif
!
      call keep_compiler_quiet(f)
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
      if (nprocx==1) then 
        if (lshear) then
          call fourier_transform_shear(phi,b1)
        else
          call fourier_transform(phi,b1)
        endif
      else 
        call fft_xy_parallel(phi,b1)
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
              k2 = (kx_fft(ikx+ipx*nx)+deltay/Lx*ky_fft(iky+ipy*ny))**2+ky_fft(iky+ipy*ny)**2
            else
              k2 = kx_fft(ikx+ipx*nx)**2+ky_fft(iky+ipy*ny)**2
            endif
!            
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
      if (nprocx==1) then 
        if (lshear) then
          call fourier_transform_shear(phi,b1,linv=.true.)
        else
          call fourier_transform(phi,b1,linv=.true.)
        endif
      else
        call fft_xy_parallel(phi,b1,linv=.true.)
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
!
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
!
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
