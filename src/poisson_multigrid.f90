! $Id$

!
!  This module solves the Poisson equation
!    (d^2/dx^2 + d^2/dy^2 + d^2/dz^2 - h) f = RHS(x,y,z)
!  [which for h/=0 could also be called inhomogenous nonuniform Helmholtz
!  equation] for the function f(x,y,z), starting from the second-order
!  accurate 7-point discretization of that equation.

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

  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages

  implicit none

  real :: dummy=0.0

  include 'poisson.h'

  namelist /poisson_init_pars/ &
      niter_poisson

  namelist /poisson_run_pars/ &
      niter_poisson
!
!  Number of iterations for multigrid solver.
!
  integer :: niter_poisson=30

  contains

!***********************************************************************
    subroutine initialize_poisson()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  18-oct-07/anders: adapted
!
    endsubroutine initialize_poisson
!***********************************************************************
    subroutine inverse_laplacian(phi)
!
!  Solve the Poisson (or inhomogeneous Helmholtz) equation
!    (\Laplace - h) f = rhs
!  using a multigrid approach. On entry, phi is rhs, on exit phi contains
!  the solution f.
!
!  18-may-2007/wolf: adapted from IDL prototype
!
      real, dimension(nx,ny,nz), intent(inout) :: phi
!
      real, dimension(nx,ny,nz) :: rhs
      real, dimension(3) :: dxyz
      integer :: iter
!
!  identify version
!
      if (lroot .and. ip<10) call svn_id( &
        "$Id$")
!
      if (lshear) then
        call fatal_error("inverse_laplacian", &
            "Never thought of shear in poisson_multigrid")
      endif
!
      rhs = phi
      dxyz = (/ dx, dy, dz /)
      do iter=1,niter_poisson
        if (lroot .and. ip<=12) then
          if (mod(iter,10) == 0) &
              write(*,'(A,I4,"/",I4)') 'Iteration ', iter, niter_poisson
        endif
        call v_cycle(phi,rhs,dxyz)
      enddo
!
    endsubroutine inverse_laplacian
! !***********************************************************************
!     subroutine inverse_curl2(phi, h, kmax)
! !
! !  Solve the elliptic equation
! !    (-\curl\curl - h) f = rhs
! !  using a  multigrid approach.
! !  On entry, phi is rhs, on exit phi contains the solution f.
! !
! !  17-jul-2007/wolf: adapted from inverse_laplacian
! !
!       real, dimension(nx,ny,nz,3)           :: phi,rhs
!       real, dimension(nx,ny,nz,3), optional :: h
!       real, optional                        :: kmax
!       real, dimension(3)                    :: dxyz
!       integer                               :: iter
! !
!       intent(in)    :: h, kmax
!       intent(inout) :: phi
! !
! !  identify version
! !
!       if (lroot .and. ip<10) call svn_id( &
!         "$Id$")
! !
!       if (present(kmax)) then
!         call warning("inverse_curl2", &
!             "kmax ignored by multigrid solver")
!       endif
! !
!       if (lshear) then
!         call fatal_error("inverse_curl2", &
!             "Never thought of shear in poisson_multigrid")
!       endif
! !
!       rhs = phi
!       dxyz = (/ dx, dy, dz /)
!       do iter=1,niter_poisson
!         if (lroot .and. ip<=12) then
!           if (mod(iter,10) == 0) &
!               write(*,'(A,I4,"/",I4)') 'Iteration ', iter, niter_poisson
!         endif
!         if (present(h)) then
!           call v_cycle_c2(phi,rhs,dxyz,h)
!         else
!           call v_cycle_c2(phi,rhs,dxyz)
!         endif
!       enddo
! !
!     endsubroutine inverse_curl2
!***********************************************************************
    recursive subroutine v_cycle(f, rhs, dxyz, h)
!
!  Do one full multigrid V cycle, i.e. go from finest to coarsest grid
!  and back:
!
!    h            h
!     \          /
!     2h        2h
!       \      /
!        4h   4h
!         \  /
!          8h
!
!  to get an approximate solution f of (\Laplace - h) f = rhs.
!
!  This subroutine calls itself recursively, so from the viewpoint of the
!  numerical grid the V cycle in question may as well be
!
!     2h        2h
!       \      /
!        4h   4h
!         \  /
!          8h
!
!  or
!
!        4h   4h
!         \  /
!          8h
!
!
!  18-may-2007/wolf: coded
!
      real, dimension(:,:,:)           :: f,rhs
      real, dimension(:,:,:), optional :: h
      real, dimension(size(f,1),size(f,2),size(f,3)) :: res,df
      real, allocatable, dimension(:,:,:) :: res_coarse, c_coarse, df_coarse
      real, dimension(3)                  :: dxyz,dxyz_coarse
      integer, dimension(3)               :: nxyz, nxyz_coarse
!
      intent(in)    :: rhs, dxyz, h
      intent(inout) :: f
!
      nxyz = (/ size(f,1), size(f,2), size(f,3) /)

      !
      ! Residual
      !
      if (present(h)) then
        call residual(f, rhs, dxyz, res, h)
      else
        call residual(f, rhs, dxyz, res)
      endif
!
! For convergence plots. Use in connection with
!   start.csh > start.log    # (or the same with run.csh)
!   ${PENCIL_HOME}/utils/wolf/extract_residuals start.log
!   idl ${PENCIL_HOME}/doc/multigrid/convergence.pro
!
      if (lroot .and. ip <= 12) then
        if (nxyz(1) == nx) write(*,*) '------------------------------'
        write(*,'(A,I4,": ",3(", ",1pG10.3)," ",L1)') &
            'nxyz(1), residual (min,max,rms), pres(h): ', &
            nxyz(1), &
            minval(res), maxval(res), sqrt(sum(res**2)/size(res)), &
            present(h)
      endif
      !
      ! Set up coarse grid
      !
      nxyz_coarse = (nxyz+1)/2  ! so 3 -> 2, not 3 -> 1
      nxyz_coarse = max(nxyz_coarse, 2)
      dxyz_coarse = dxyz * (nxyz - 1) / (nxyz_coarse - 1)
      !
      allocate(res_coarse(nxyz_coarse(1), nxyz_coarse(2), nxyz_coarse(3)))
      allocate(c_coarse(  nxyz_coarse(1), nxyz_coarse(2), nxyz_coarse(3)))
      allocate(df_coarse( nxyz_coarse(1), nxyz_coarse(2), nxyz_coarse(3)))
      !
      ! ..and restrict to coarse grid
      !
      if (.not. allocated(res_coarse)) &
          call fatal_error("v_cycle", "res_coarse is unallocated")
      call restrict(res, nxyz, nxyz_coarse, res_coarse)
      !
      if (present(h)) then
        if (.not. allocated(c_coarse)) &
            call fatal_error("v_cycle", "c_coarse is unallocated")
        call restrict(h,   nxyz, nxyz_coarse, c_coarse  )
      endif

      !
      ! Do one relaxation step (i.e. damp larger-scale errors on coarse grid)
      !
      if (.not. allocated(df_coarse)) &
          call fatal_error("v_cycle", "df_coarse is unallocated")
      df_coarse = 0.            ! initial value

      if (any(nxyz > 2)) then
        if (present(h)) then
          call v_cycle(df_coarse, res_coarse, dxyz_coarse, c_coarse)
        else
          call v_cycle(df_coarse, res_coarse, dxyz_coarse)
        endif
      endif

      !
      ! Refine (interpolate) coarse-grid correction to finer grid, and
      ! apply
      !
      ! TODO:
      !   Can be memory-optimized to immediately add df(i,j,k) to f(i,j,k)
      !   -> no df needed
      call refine(df_coarse, nxyz_coarse, nxyz, df)
      f = f - df

      !
      ! Do one smoothing iteration (i.e. damp small-scale errors)
      !
      if (present(h)) then
        call gauss_seidel_iterate(f, rhs, dxyz, h)
      else
        call gauss_seidel_iterate(f, rhs, dxyz)
      endif
!
      deallocate(res_coarse)
      deallocate(c_coarse)
      deallocate(df_coarse)
!
    endsubroutine v_cycle
!***********************************************************************
    subroutine residual(f, rhs, dxyz, res, h)
!
!  Calculate the residual r = L*f - rhs
!
!  19-may-2007/wolf: coded
!
      real, dimension(:,:,:)           :: f,rhs,res
      real, dimension(:,:,:), optional :: h
      real, dimension(3)               :: dxyz
      integer                          :: nnx,nny,nnz
!
      intent(in)  :: f, rhs, dxyz, h
      intent(out) :: res
!
      nnx = size(f,1)
      nny = size(f,2)
      nnz = size(f,3)
      !
      ! Second-order Laplacian
      !
      res = 0.                  ! avoid floating exception below
      res(2:nnx-1, 2:nny-1, 2:nnz-1) &
          = (    f(1:nnx-2, 2:nny-1, 2:nnz-1) &
             - 2*f(2:nnx-1, 2:nny-1, 2:nnz-1) &
             +   f(3:nnx-0, 2:nny-1, 2:nnz-1) &
            ) / dxyz(1)**2 &
      !
          + (    f(2:nnx-1, 1:nny-2, 2:nnz-1) &
             - 2*f(2:nnx-1, 2:nny-1, 2:nnz-1) &
             +   f(2:nnx-1, 3:nny-0, 2:nnz-1) &
            ) / dxyz(2)**2 &
      !
          + (    f(2:nnx-1, 2:nny-1, 1:nnz-2) &
             - 2*f(2:nnx-1, 2:nny-1, 2:nnz-1) &
             +   f(2:nnx-1, 2:nny-1, 3:nnz-0) &
            ) / dxyz(3)**2
      !
      res = res - rhs

      if (present(h)) then
        res = res - h*f
      endif

      !
      !  Make sure we do not change f at the boundaries:
      !
      res(1,:,:) = 0.;  res(nnx,  :,  :) = 0.
      res(:,1,:) = 0.;  res(  :,nny,  :) = 0.
      res(:,:,1) = 0.;  res(  :,  :,nnz) = 0.
!
    endsubroutine residual
!***********************************************************************
    subroutine restrict(var, nxyz, nxyz_coarse, var_coarse)
!
!  Restrict, (i.e. coarse-grain) var from the original to coarser grid.
!  We use `full weighting' on the original grid, followed by trilinear
!  interpolation to the (cell-centered) coarse grid we really want.
!
!  This could be `chunked' (i.e. applied in blocks) to save some memory.
!
!  19-may-2007/wolf: coded
!
      real, dimension(:,:,:)          :: var, var_coarse
      integer, dimension(3)           :: nxyz, nxyz_coarse
      real, dimension(size(var,1),size(var,2),size(var,3)) :: tmp
!
      intent(in)  :: var, nxyz, nxyz_coarse
      intent(out) :: var_coarse
      !
      ! Smooth, then interpolate
      !
      tmp = var
      call smooth_full_weight(tmp)
      call trilinear_interpolate(tmp, nxyz, nxyz_coarse, var_coarse)
!
    endsubroutine restrict
!***********************************************************************
    subroutine smooth_full_weight(f)
!
!  Apply `full weighting' smoothing ([1,2,1]/4 weighting for each
!  direction) to var.
!
!  Full-weighting smoothing filters out the Nyquist frequency in x, y, z,
!  diagonal and space diagonal directions. For vertex-centered multigrid,
!  this corresponds directly to the restricton operator (if downsampled to
!  those points that are part of the coarser grid); for cell-centered
!  multigrid (our case), we will need to interpolate.
!
!  19-may-2007/wolf: coded
!
      real, dimension(:,:,:) :: f
      integer                :: nnx, nny, nnz
!
      intent(inout) :: f
!
      nnx = size(f,1)
      nny = size(f,2)
      nnz = size(f,3)
!
      if (nnx > 2) then
        f(2:nnx-1, :, :)  &
            = 0.25 * (    f(1:nnx-2, :, :) &
                      + 2*f(2:nnx-1, :, :) &
                      +   f(3:nnx-0, :, :) &
                     )
      endif
      !
      if (nny > 2) then
        f(:, 2:nny-1, :) &
            = 0.25 * (    f(: ,1:nny-2, :) &
                      + 2*f(: ,2:nny-1, :) &
                      +   f(: ,3:nny-0, :) &
                     )
      endif
      !
      if (nnz > 2) then
        f(:, :, 2:nnz-1) &
            = 0.25 * (    f(:, :, 1:nnz-2) &
                      + 2*f(:, :, 2:nnz-1) &
                      +   f(:, :, 3:nnz-0) &
                     )
      endif
!
    endsubroutine smooth_full_weight
!***********************************************************************
    subroutine trilinear_interpolate(var1, nxyz1, nxyz2, var2)
!
!  Trilinear interpolation from grid with nxyz1 grid points to one with
!  nxyz2 points, assuming either periodicity (NOT IMPLEMENTED YET), or
!  identical end points for the grids (e.g. we can assume both grids to
!  cover [0,1] x [0,1] x [0,1])
!
!  20-may-2007/wolf: coded
!
      real, dimension(:,:,:) :: var1, var2
      integer, dimension(3)  :: nxyz1, nxyz2
      real                   :: idx,idy,idz, rx,ry,rz, safety_factor
      integer                :: ix,iy,iz
      integer                :: jx,jy,jz
!
      intent(in)  :: var1, nxyz1, nxyz2
      intent(out) :: var2
!
! wd: TO DO: Precalculate idxx(:), rxx1(:), etc.
!
      do jz=1,nxyz2(3)
        do jy=1,nxyz2(2)
          do jx=1,nxyz2(1)
            !
            ! Construct floating-point indices on grid1
            !
            safety_factor = 1 - 2*epsilon(1.) ! ensure ix<nxyz(1), etc
            idx = 1 + (jx-1) * (nxyz1(1)-1.) / (nxyz2(1)-1.) * safety_factor
            idy = 1 + (jy-1) * (nxyz1(2)-1.) / (nxyz2(2)-1.) * safety_factor
            idz = 1 + (jz-1) * (nxyz1(3)-1.) / (nxyz2(3)-1.) * safety_factor
            !
            ! Split into integer (position of lower left neighbour) and
            ! fractional part (relative distance from that neighbour)
            !
            ix = floor(idx);  rx = idx - ix
            iy = floor(idy);  ry = idy - iy
            iz = floor(idz);  rz = idz - iz
            !
            if (ix >= nxyz1(1)) call fatal_error("trilin", &
                "I told you this would give trouble 1")
            if (iy >= nxyz1(2)) call fatal_error("trilin", &
                "I told you this would give trouble 2")
            if (iz >= nxyz1(3)) call fatal_error("trilin", &
                "I told you this would give trouble 3")
            !
            ! Interpolation formula
            !
            var2(jx,jy,jz) &
                =   var1(ix  , iy  , iz  ) * (1-rx) * (1-ry) * (1-rz) &
                !
                  + var1(ix+1, iy  , iz  ) * rx     * (1-ry) * (1-rz) &
                  + var1(ix  , iy+1, iz  ) * (1-rx) * ry     * (1-rz) &
                  + var1(ix  , iy  , iz+1) * (1-rx) * (1-ry) * rz     &
                  !
                  + var1(ix+1, iy+1, iz  ) * rx     * ry     * (1-rz) &
                  + var1(ix  , iy+1, iz+1) * (1-rx) * ry     * rz     &
                  + var1(ix+1, iy  , iz+1) * rx     * (1-ry) * rz     &
                  !
                  + var1(ix+1, iy+1, iz+1) * rx     * ry     * rz
          enddo
        enddo
      enddo
!
    endsubroutine trilinear_interpolate
!***********************************************************************
    subroutine refine(var, nxyz, nxyz_fine, var_fine)
!
!  Refine (i.e. fine-grain by interpolating) f from original to finer grid
!
!  20-may-2007/wolf:coded
!
      real, dimension(:,:,:) :: var, var_fine
      integer, dimension(3)  :: nxyz, nxyz_fine
!
      intent(in)  :: var, nxyz, nxyz_fine
      intent(out) :: var_fine
!

      !
      ! We simply use trilinear interpolation, which works fine
      !
      call trilinear_interpolate(var, nxyz, nxyz_fine, var_fine)

!
    endsubroutine refine
!***********************************************************************
    subroutine gauss_seidel_iterate(f, rhs, dxyz, h)
!
!  One Gauss-Seidel iteration sweep for the Laplace equation (with h)
!
!  20-may-2007/wolf: coded
!
      real, dimension(:,:,:)           :: f, rhs
      real, dimension(:,:,:), optional :: h
      real, dimension(3)     :: dxyz
      real                   :: dx_2,dy_2,dz_2, sum
      real                   :: denom1, denom
      integer                :: nnx,nny,nnz, ix,iy,iz
!
      intent(in)    :: rhs,dxyz,h
      intent(inout) :: f
!
      nnx = size(f,1)
      nny = size(f,2)
      nnz = size(f,3)
!
      dx_2 = 1./dxyz(1)**2
      dy_2 = 1./dxyz(2)**2
      dz_2 = 1./dxyz(3)**2

      ! Be clever about one (or more) of nnx, nny, nnz = 2: just discard the
      ! corresponding second derivative from the Laplacian.
      if (nnx <= 2)  dx_2 = 0.
      if (nny <= 2)  dy_2 = 0.
      if (nnz <= 2)  dz_2 = 0.

      call apply_boundcond(f)

      denom1 = 2*(dx_2+dy_2+dz_2)
      if (any( (/ nnx, nny, nnz /) > 2)) then
        do ix=2,nnx-1
          do iy=2,nny-1
            do iz=2,nnz-1
              sum = -rhs(ix,iy,iz)
              if (nnx > 2) then
                sum = sum + (f(ix-1,iy  ,iz  ) + f(ix+1,iy  ,iz  )) * dx_2
              endif
              if (nny > 2) then
                sum = sum + (f(ix  ,iy-1,iz  ) + f(ix  ,iy+1,iz  )) * dy_2
              endif
              if (nnz > 2) then
                sum = sum + (f(ix  ,iy  ,iz-1) + f(ix  ,iy  ,iz+1)) * dz_2
              endif
              !
              if (present(h)) then
                denom = denom1 + h(ix,iy,iz)
              else
                denom = denom1
              endif
              !
              f(ix,iy,iz) = sum / denom
            enddo
          enddo
        enddo
      endif
!
    endsubroutine gauss_seidel_iterate
!***********************************************************************
    subroutine apply_boundcond(f)
!
!  Apply boundary conditions -- currently Dirichlet cond. f=0
!
!  20-may-2007/wolf: coded
!
      real, dimension(:,:,:) :: f
      integer                :: nnx,nny,nnz
!
      intent(out) :: f
!
      nnx = size(f,1)
      nny = size(f,2)
      nnz = size(f,3)
!
      f(1  ,:  ,:) = 0
      f(nnx,:  ,:) = 0
      f(:  ,1  ,:) = 0
      f(:  ,nny,:) = 0
      f(:  ,:,1  ) = 0
      f(:  ,:,nnz) = 0
!
    endsubroutine apply_boundcond
!***********************************************************************
    subroutine inverse_laplacian_semispectral(phi,h)
!
!  Solve the Poisson equation by Fourier transforming in the xy-plane and
!  solving the discrete matrix equation in the z-direction.
!
      real, dimension (nx,ny,nz) :: phi
      real, optional             :: h
!
      call fatal_error('inverse_laplacian_semispectral', &
          'Cheating! This is not multigrid we are using')
!
      call keep_compiler_quiet(phi)
      call keep_compiler_quiet(h)
!
    endsubroutine inverse_laplacian_semispectral
!***********************************************************************
    subroutine inverse_laplacian_fft_z(phi)
!
!  15-may-2006/anders+jeff: dummy
!
      real, dimension(nx,ny,nz), intent(in) :: phi
!
      call keep_compiler_quiet(phi)
!
    endsubroutine inverse_laplacian_fft_z
!***********************************************************************
    subroutine inverse_laplacian_z_2nd_neumann(f)
!
!  15-may-2006/anders+jeff: dummy
!
      real, dimension(:,:,:,:), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine inverse_laplacian_z_2nd_neumann
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
      real, dimension(nx,ny,nz,3), intent(out) :: acceleration           !should I (CAN I?) make this allocatable?
!
      call keep_compiler_quiet(acceleration)
!
    endsubroutine get_acceleration
!***********************************************************************
endmodule Poisson
