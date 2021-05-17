! $Id$
!
!  This module provides general facilities to implicitly solve a
!  diffusion equation.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: limplicit_diffusion = .true.
!
!***************************************************************
module ImplicitDiffusion
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: fatal_error, warning
!
  implicit none
!
  include 'implicit_diffusion.h'
!
!  Run-time parameters
!
  character(len=6) :: implicit_method = 'full'
!
  namelist /implicit_diffusion_run_pars/ implicit_method
!
!  Module-specific variables
!
  real, dimension(nxgrid) :: ax_imp = 0.0, bx_imp = 0.0, cx_imp = 0.0
  real, dimension(nygrid) :: ay_imp = 0.0, by_imp = 0.0, cy_imp = 0.0
  real, dimension(nzgrid) :: az_imp = 0.0, bz_imp = 0.0, cz_imp = 0.0
!
  integer, parameter :: BOT=1, TOP=2
!
  contains
!***********************************************************************
!***********************************************************************
!  PUBLIC ROUTINES GO BELOW HERE.
!***********************************************************************
!***********************************************************************
    subroutine read_implicit_diff_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=implicit_diffusion_run_pars, iostat=iostat)
!
    endsubroutine read_implicit_diff_run_pars
!***********************************************************************
    subroutine write_implicit_diff_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, nml=implicit_diffusion_run_pars)
!
    endsubroutine write_implicit_diff_run_pars
!***********************************************************************
    subroutine integrate_diffusion(get_diffus_coeff, f, ivar1, ivar2)
!
! Integrate the diffusion term for components ivar1 to ivar2 according
! to implicit_method.  ivar2 can be omitted if only operating on one
! component.
!
! 05-sep-14/ccyang: coded.
!
      external :: get_diffus_coeff
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: ivar1
      integer, intent(in), optional :: ivar2
!
      integer :: iv1, iv2
!
!  Check the indices.
!
      iv1 = ivar1
      comp: if (present(ivar2)) then
        iv2 = ivar2
      else comp
        iv2 = ivar1
      endif comp
!
!  Dispatch.
!
      method: select case (implicit_method)
!
      case ('full') method
        call integrate_diffusion_full(get_diffus_coeff, f, iv1, iv2)
!
      case ('fft') method
        call integrate_diffusion_fft(get_diffus_coeff, f, iv1, iv2)
!
      case ('zonly') method
        call integrate_diffusion_zonly(get_diffus_coeff, f, iv1, iv2)
!
      case default method
        call fatal_error('integrate_diffusion', 'unknown implicit_method = ' // implicit_method)
!
      endselect method
!
    endsubroutine integrate_diffusion
!***********************************************************************
    subroutine integrate_diffusion_full(get_diffus_coeff, f, ivar1, ivar2)
!
! Implicitly integrate the diffusion term for components ivar1 to ivar2
! by full dimensional splitting.
!
! 04-sep-13/ccyang: coded.
!
      use Shear, only: sheared_advection_fft
!
      external :: get_diffus_coeff
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: ivar1, ivar2
!
      real :: dth
!
!  Sanity check
!
      if (nprocx /= 1) call fatal_error('integrate_diffusion_full', 'nprocx /= 1')
      if (nygrid > 1 .and. nxgrid /= nygrid) &
        call fatal_error('integrate_diffusion_full', 'nxgrid /= nygrid')
      if (nzgrid > 1 .and. mod(nxgrid, nprocz) /= 0) &
        call fatal_error('integrate_diffusion_full', 'mod(nxgrid,nprocz) /= 0')
!
!  Get half the time step.
!
      dth = 0.5 * dt
!
!  Integrate.
!
      shearing: if (lshear) then
!       Shear is present: Work in unsheared frame.
        call zsweep(f, bcz12, ivar1, ivar2, dth, get_diffus_coeff)
        call ysweep(f, spread(spread('p', 1, mcom),1,2), ivar1, ivar2, dth, get_diffus_coeff)
!
        call sheared_advection_fft(f, ivar1, ivar2, real(-t))
        !!call xsweep(f, spread('p', 1, mcom), spread('p', 1, mcom), iv1, iv2, dth, get_diffus_coeff)
        !!call xsweep(f, spread('p', 1, mcom), spread('p', 1, mcom), iv1, iv2, dth, get_diffus_coeff)
        call xsweep(f, spread(spread('p', 1, mcom),1,2), ivar1, ivar2, dth, get_diffus_coeff)
        call xsweep(f, spread(spread('p', 1, mcom),1,2), ivar1, ivar2, dth, get_diffus_coeff)
        call sheared_advection_fft(f, ivar1, ivar2, real(t))
!
        !!call ysweep(f, spread('p', 1, mcom), spread('p', 1, mcom), iv1, iv2, dth, get_diffus_coeff)
        call ysweep(f, spread(spread('p', 1, mcom),1,2), ivar1, ivar2, dth, get_diffus_coeff)
        call zsweep(f, bcz12, ivar1, ivar2, dth, get_diffus_coeff)
      else shearing
!       No shear is present.
        call xsweep(f, bcx12, ivar1, ivar2, dth, get_diffus_coeff)
        call ysweep(f, bcy12, ivar1, ivar2, dth, get_diffus_coeff)
        call zsweep(f, bcz12, ivar1, ivar2, dth, get_diffus_coeff)
!
        call zsweep(f, bcz12, ivar1, ivar2, dth, get_diffus_coeff)
        call ysweep(f, bcy12, ivar1, ivar2, dth, get_diffus_coeff)
        call xsweep(f, bcx12, ivar1, ivar2, dth, get_diffus_coeff)
      endif shearing
!
    endsubroutine integrate_diffusion_full
!***********************************************************************
    subroutine integrate_diffusion_fft(get_diffus_coeff, f, ivar1, ivar2)
!
! Integrate the diffusion term exactly for components ivar1 to ivar2 by
! Fourier decomposition.  A constant diffusion coefficient is assumed.
!
! 25-sep-14/ccyang: coded.
!
      use Fourier, only: fft_xyz_parallel, kx_fft, ky_fft, kz_fft
!
      interface
        subroutine get_diffus_coeff(ndc, dc, iz)
          integer, intent(in) :: ndc
          real, dimension(ndc), intent(out) :: dc
          integer, intent(in), optional :: iz
        endsubroutine
      endinterface
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: ivar1, ivar2
!
      real, dimension(nx,ny,nz) :: a_re, a_im
      real, dimension(nx) :: decay, k2dt
      real, dimension(1) :: dc
      integer :: ll1, ll2, m0, n0
      integer :: iv, j, k
      real :: ky2, kz2, c
!
! Get the diffusion coefficient.
!
      call get_diffus_coeff(1, dc)
!
! Integrate.
!
      ll1 = ipx * nx + 1
      ll2 = (ipx + 1) * nx
      m0 = ipy * ny
      n0 = ipz * nz
      if (lshear) c = deltay / Lx
      comp: do iv = ivar1, ivar2
        a_re = f(l1:l2,m1:m2,n1:n2,iv)
        a_im = 0.0
        call fft_xyz_parallel(a_re, a_im)
        zscan: do k = 1, nz
          kz2 = kz_fft(n0+k)**2
          yscan: do j = 1, ny
            ky2 = ky_fft(m0+j)**2
            if (lshear) then
              k2dt = dt * (ky2 + kz2 + (kx_fft(ll1:ll2) + c * ky_fft(m0+j))**2)
            else
              k2dt = dt * (ky2 + kz2 + kx_fft(ll1:ll2)**2)
            endif
            decay = exp(-dc(1) * k2dt)
            a_re(:,j,k) = decay * a_re(:,j,k)
            a_im(:,j,k) = decay * a_im(:,j,k)
          enddo yscan
        enddo zscan
        call fft_xyz_parallel(a_re, a_im, linv=.true.)
        f(l1:l2,m1:m2,n1:n2,iv) = a_re
      enddo comp
!
    endsubroutine integrate_diffusion_fft
!***********************************************************************
    subroutine integrate_diffusion_zonly(get_diffus_coeff, f, ivar1, ivar2)
!
! Integrate the diffusion term for components ivar1 to ivar2 by Fourier
! decomposition horizontally and implicit solution vertically.  It is
! assumed that the diffusion coefficient does not depend on x or y.
!
! 05-sep-14/ccyang: coded.
!
      external :: get_diffus_coeff
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: ivar1, ivar2
!
      call zsweep(f, bcz12, ivar1, ivar2, 0.5 * dt, get_diffus_coeff)
      call integrate_diffusion_fft_xy(get_diffus_coeff, f, ivar1, ivar2)
      call zsweep(f, bcz12, ivar1, ivar2, 0.5 * dt, get_diffus_coeff)
!
    endsubroutine integrate_diffusion_zonly
!***********************************************************************
!
!  LOCAL ROUTINES GO BELOW HERE.
!***********************************************************************
!***********************************************************************
    subroutine integrate_diffusion_fft_xy(get_diffus_coeff, f, ivar1, ivar2)
!
! Integrate the diffusion term exactly for components ivar1 to ivar2 by
! Fourier decomposition in the x and y directions.  A horizontally
! constant diffusion coefficient is assumed.
!
! 05-sep-14/ccyang: coded.
!
      use Fourier, only: fft_xy_parallel, kx_fft, ky_fft
!
      interface
        subroutine get_diffus_coeff(ndc, dc, iz)
          integer, intent(in) :: ndc
          real, dimension(ndc), intent(out) :: dc
          integer, intent(in), optional :: iz
        endsubroutine
      endinterface
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: ivar1, ivar2
!
      real, dimension(nx,ny,nz) :: a_re, a_im
      real, dimension(nx) :: k2dt, decay
      real, dimension(nzgrid) :: dc
      integer :: ll1, ll2, m0, n0
      integer :: iv, j, k
      real :: c
!
! Get the diffusion coefficient.
!
      call get_diffus_coeff(nzgrid, dc)
!
! Integrate.
!
      ll1 = ipx * nx + 1
      ll2 = (ipx + 1) * nx
      m0 = ipy * ny
      n0 = ipz * nz
      if (lshear) c = deltay / Lx
      comp: do iv = ivar1, ivar2
        a_re = f(l1:l2,m1:m2,n1:n2,iv)
        a_im = 0.0
        call fft_xy_parallel(a_re, a_im)
        yscan: do j = 1, ny
          if (lshear) then
            k2dt = dt * ((kx_fft(ll1:ll2) + c * ky_fft(m0+j))**2 + ky_fft(m0+j)**2)
          else
            k2dt = dt * (kx_fft(ll1:ll2)**2 + ky_fft(m0+j)**2)
          endif
          zscan: do k = 1, nz
            decay = exp(-dc(n0+k) * k2dt)
            a_re(:,j,k) = decay * a_re(:,j,k)
            a_im(:,j,k) = decay * a_im(:,j,k)
          enddo zscan
        enddo yscan
        call fft_xy_parallel(a_re, a_im, linv=.true.)
        f(l1:l2,m1:m2,n1:n2,iv) = a_re
      enddo comp
!
    endsubroutine integrate_diffusion_fft_xy
!***********************************************************************
    subroutine set_diffusion_equations(get_diffus_coeff, direction, a, b, c, iz)
!
! Determines the coefficients a_i, b_i, and c_i of the discretized
! diffusion equations:
!
!   d q_i / dt = psi(q_{i-1},q_i,q_{i+1})
!
!              = a_i * q_{i-1} + b_i * q_i + c_i * q_{i+1},
!
! where psi is a spatially discretized diffusion operator and i runs
! from 1 to nxgrid, nygrid, or nzgrid according to specified direction.
! Input argument direction is an integer indicating the direction of the
! diffusive operation.  Input procedure argument get_diffus_coeff is
! a subroutine returning the diffusion coefficient with the interface
!
!   subroutine get_diffus_coeff(ndc, diffus_coeff, iz)
!
!     integer, intent(in) :: ndc
!     real, dimension(ndc), intent(out) :: diffus_coeff
!     integer, intent(in), optional :: iz
!
! where diffus_coeff(i) is the diffusion coefficient at i.  If present,
! the optional argument iz is the global z-index for diffus_coeff.  If
! not present, diffus_coeff is assumed to be vertical.
!
! 16-aug-13/ccyang: coded
! 23-jul-14/ccyang: introduced positional dependence
!
      interface
        subroutine get_diffus_coeff(ndc, dc, iz)
          integer, intent(in) :: ndc
          real, dimension(ndc), intent(out) :: dc
          integer, intent(in), optional :: iz
        endsubroutine
      endinterface
!
      integer, intent(in) :: direction
      real, dimension(:), intent(out) :: a, b, c
      integer, intent(in), optional :: iz
!
      real, dimension(:), allocatable :: dc
!
! Find the coefficients.
!
      dir: select case (direction)
!
      case (1) dir
        allocate(dc(nxgrid))
        call get_diffus_coeff(nxgrid, dc, iz)
        a = 0.5 * dc * dx1grid * (dx1grid - 0.5 * dxtgrid)
        b = -dc * dx1grid**2
        c = 0.5 * dc * dx1grid * (dx1grid + 0.5 * dxtgrid)
!
      case (2) dir
        allocate(dc(nygrid))
        call get_diffus_coeff(nygrid, dc, iz)
        a = 0.5 * dc * dy1grid * (dy1grid - 0.5 * dytgrid)
        b = -dc * dy1grid**2
        c = 0.5 * dc * dy1grid * (dy1grid + 0.5 * dytgrid)
!
      case (3) dir
        allocate(dc(nzgrid))
        call get_diffus_coeff(nzgrid, dc)
        a = 0.5 * dc * dz1grid * (dz1grid - 0.5 * dztgrid)
        b = -dc * dz1grid**2
        c = 0.5 * dc * dz1grid * (dz1grid + 0.5 * dztgrid)
!
      case default dir
        call fatal_error('set_up_equations', 'unknown direction')
!
      endselect dir
!
      deallocate(dc)
!
    endsubroutine set_diffusion_equations
!***********************************************************************
    subroutine xsweep(f, bcx, ivar1, ivar2, dt, get_diffus_coeff)
!
! Implicitly integrate the diffusion term in the x-direction.
!
! 16-aug-13/ccyang: coded
! 23-jul-14/ccyang: introduced positional dependence
!
! Input Arguments
!   bcx1: array of the lower boundary conditions for each component
!   bcx2: array of the upper boundary conditions for each component
!   ivar1: beginning f component
!   ivar2: ending f component
!   dt: time step
! Input Procedure
!   subroutine get_diffus_coeff(ndc, diffus_coeff): see set_diffusion_equations
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      character(len=*), dimension(mcom,2), intent(in) :: bcx
      integer, intent(in) :: ivar1, ivar2
      real, intent(in) :: dt
      external :: get_diffus_coeff
!
      real, dimension(nxgrid) :: a, b, c
      real, dimension(nxgrid) :: adt, opbdt, ombdt, cdt
      integer :: j, k, l
!
      int_x: if (nxgrid > 1) then
        zscan: do k = n1, n2
          call set_diffusion_equations(get_diffus_coeff, 1, a, b, c, iz=ipz*nz+k-nghost)
          call get_tridiag(a, b, c, dt, adt, opbdt, ombdt, cdt)
          yscan: do j = m1, m2
            do l = ivar1, ivar2
              call implicit_pencil( f(l1-1:l2+1,j,k,l), nxgrid, adt, opbdt, ombdt, cdt, &
                                    bcx(l,:), dx2_bound(-1:1), xgrid((/1,nxgrid/)) )
            enddo
          enddo yscan
        enddo zscan
      endif int_x
!
    endsubroutine xsweep
!***********************************************************************
    subroutine ysweep(f, bcy, ivar1, ivar2, dt, get_diffus_coeff)
!
! Implicitly integrate the diffusion term in the y-direction.
!
! 16-aug-13/ccyang: coded
! 23-jul-14/ccyang: introduced positional dependence
!
! Input Arguments
!   bcy1: array of the lower boundary conditions for each component
!   bcy2: array of the upper boundary conditions for each component
!   ivar1: beginning f component
!   ivar2: ending f component
!   dt: time step
! Input Procedure
!   subroutine get_diffus_coeff(ndc, diffus_coeff): see set_diffusion_equations
!
      use Mpicomm, only: transp_xy
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      character(len=*), dimension(mcom,2), intent(in) :: bcy
      integer, intent(in) :: ivar1, ivar2
      real, intent(in) :: dt
      external :: get_diffus_coeff
!
      real, dimension(nx,ny) :: axy      ! assuming nxgrid = nygrid and nprocx = 1
      real, dimension(0:nx+1) :: penc    ! assuming nxgrid = nygrid and nprocx = 1
      real, dimension(nygrid) :: a, b, c
      real, dimension(nygrid) :: adt, opbdt, ombdt, cdt
      real, dimension(2) :: bound, boundl
      real, dimension(-1:1) :: d2_bound

      integer :: j, k, l
      real :: fac
!
      int_y: if (nygrid > 1) then
!
        bound = ygrid((/1,nygrid/))
        d2_bound(-1)= 2.*(ygrid(2)-ygrid(1))
        d2_bound( 1)= 2.*(ygrid(nygrid)-ygrid(nygrid-1))
!
        comp: do l = ivar1, ivar2

          boundl=bound
          if ( bcy(l,1)(2:)=='fr' ) boundl(1) = tan(bound(1))
          if ( bcy(l,2)(2:)=='fr' ) boundl(2) = tan(bound(2))

          zscan: do k = n1, n2
            call set_diffusion_equations(get_diffus_coeff, 2, a, b, c, iz=ipz*nz+k-nghost)
            if (.not.lspherical_coords) &
              call get_tridiag(a, b, c, dt, adt, opbdt, ombdt, cdt)
            axy = f(l1:l2,m1:m2,k,l)
            call transp_xy(axy)  ! assuming nxgrid = nygrid
            xscan: do j = 1, ny
              if (lspherical_coords) then
                fac=1./xgrid(j)**2
                call get_tridiag(fac*a, fac*b, fac*c, dt, adt, opbdt, ombdt, cdt)
              endif
              penc(1:nx) = axy(:,j)
              call implicit_pencil( penc, nygrid, adt, opbdt, ombdt, cdt, bcy(l,:), &
                                    d2_bound, boundl )
              axy(:,j) = penc(1:nx)
            enddo xscan
            call transp_xy(axy)
            f(l1:l2,m1:m2,k,l) = axy
          enddo zscan
        enddo comp
      endif int_y
!
    endsubroutine ysweep
!***********************************************************************
    subroutine zsweep(f, bcz, ivar1, ivar2, dt, get_diffus_coeff)
!
! Implicitly integrate the diffusion term in the z-direction.
!
! 30-sep-14/ccyang: coded
!
! Input Arguments
!   bcz1: array of the lower boundary conditions for each component
!   bcz2: array of the upper boundary conditions for each component
!   ivar1: beginning f component
!   ivar2: ending f component
!   dt: time step
! Input Procedure
!   subroutine get_diffus_coeff(ndc, diffus_coeff): see set_diffusion_equations
!
      use Mpicomm, only: remap_to_pencil_yz, unmap_from_pencil_yz
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      character(len=*), dimension(mcom,2), intent(in) :: bcz
      integer, intent(in) :: ivar1, ivar2
      real, intent(in) :: dt
      external :: get_diffus_coeff
!
      integer, parameter :: nyt = ny / nprocz    !!! dangerous
      real, dimension(nx,ny,nz) :: ff
      real, dimension(nx,nyt,nzgrid) :: ft
      real, dimension(0:nzgrid+1) :: penc
      real, dimension(nzgrid) :: a, b, c
      real, dimension(nzgrid) :: adt, opbdt, ombdt, cdt
      integer :: i, j, iv
      real :: fac, facy
      real, dimension(2) :: bound
      real, dimension(-1:1) :: d2_bound
!
      int_z: if (nzgrid > 1) then
!
        bound = zgrid((/1,nzgrid/))
        d2_bound(-1)= 2.*(zgrid(2)-zgrid(1))
        d2_bound( 1)= 2.*(zgrid(nzgrid)-zgrid(nzgrid-1))
!
        call set_diffusion_equations(get_diffus_coeff, 3, a, b, c)
        if (.not.lspherical_coords) &
          call get_tridiag(a, b, c, dt, adt, opbdt, ombdt, cdt)
        comp: do iv = ivar1, ivar2
          call remap_to_pencil_yz(f(l1:l2,m1:m2,n1:n2,iv), ft)
          yscan: do j = 1, nyt
            if (lspherical_coords) facy=1./sin(ygrid(j))**2
            xscan: do i = 1, nx
              if (lspherical_coords) then
                fac=facy/xgrid(i)**2
                call get_tridiag(fac*a, fac*b, fac*c, dt, adt, opbdt, ombdt, cdt)
              endif
              penc(1:nzgrid) = ft(i,j,:)
              call implicit_pencil( penc, nzgrid, adt, opbdt, ombdt, cdt, bcz(iv,:),  &
                                    dz2_bound, bound )
              ft(i,j,:) = penc(1:nzgrid)
            enddo xscan
          enddo yscan
          call unmap_from_pencil_yz(ft, ff)
          f(l1:l2,m1:m2,n1:n2,iv) = ff
        enddo comp
      endif int_z
!
    endsubroutine zsweep
!***********************************************************************
    subroutine get_tridiag(a, b, c, dt, adt, opbdt, ombdt, cdt)
!
! Prepare the tridiagonal elements of the linear system for the Crank-
! Nicolson algorithm.
!
! 22-may-12/ccyang: coded
!
! Input Arguments
!   a, b, c: a_i, b_i, and c_i, respectively, defined in the
!     subroutine implicit_pencil except a factor of time step dt
!   dt: time step
!
! Output Arguments
!   adt  : array of dt * a_i
!   opbdt: array of 1 + dt * b_i
!   ombdt: array of 1 - dt * b_i
!   cdt  : array of dt * c_i
!
      real, dimension(:), intent(in) :: a, b, c
      real, dimension(:), intent(out) :: adt, opbdt, ombdt, cdt
      real, intent(in) :: dt
!
      adt = dt * a
      cdt = dt * b
      opbdt = 1.0 + cdt
      ombdt = 1.0 - cdt
      cdt = dt * c
!
    endsubroutine get_tridiag
!***********************************************************************
    subroutine implicit_pencil(q, n, a, opb, omb, c, bc12, d2_bound, bound)
!
! Implicitly integrate a state variable along one pencil using the
! Crank-Nicolson method.  The (linear) system of equations read
!
!   q_i^{n+1} - q_i^n = psi(q_{i-1}^n, q_i^n, q_{i+1}^n)
!                     + psi(q_{i-1}^{n+1}, q_i^{n+1}, q_{i+1}^{n+1})
!
! where q_i^n is the value of q at x = x_i and t = t_n and
!
!   psi(q_{i-1},q_i,q_{i+1}) = a_i * q_{i-1} + b_i * q_i + c_i * q_{i+1}.
!
! 30-may-12/ccyang: coded
!  1-apr-15/MR: parameters d2_bound and bound added
! 
! Comments:
!   Although this subroutine is currently only used by Magnetic module, it is
! general enough such that it can be considered moving to Sub or General
! module.
!
! Input/Output Arguments
!   q - state variable along one pencil including one ghost cell on each side
!
! Input Arguments
!   n - number of active cells
!   a - array of a_i
!   opb - array of 1 + b_i
!   omb - array of 1 - b_i
!   c - array of c_i
!   bc1 - lower boundary conditions
!   bc2 - upper boundary conditions
!   d2_bound - doubled cumulative cell sizes at boundary
!              (d2_bound(-nghost:-1) - at lower, d2_bound(1:nghost) - at upper)
!   bound - boundary coordinates
!
      use Boundcond, only: bc_pencil
      use General, only: cyclic, tridag
!
      integer, intent(in) :: n
      real, dimension(0:n+1), intent(inout) :: q
      real, dimension(n), intent(in) :: a, opb, omb, c
      character(len=*), dimension(2), intent(in) :: bc12
      real, dimension(-1:1), optional :: d2_bound
      real, dimension(2), optional :: bound
!
      real, dimension(n) :: r, omb1, a1, c1
      character(len=255) :: msg
      logical :: lcyclic, err
      real :: alpha, beta
!
! Prepare the RHS of the linear system.
!
      call bc_pencil(q, n, 1, bc12, d2_bound, bound)
      r = a * q(0:n-1) + opb * q(1:n) + c * q(2:n+1)
!
! Assume non-periodic boundary conditions first.
!
      lcyclic = .false.
      alpha = 0.0
      beta = 0.0
      omb1 = omb; a1=a; c1=c
!
! Check the lower boundary conditions.
!
      lower: select case (bc12(1))
!     Periodic
      case ('p') lower
        beta = -a(1)
        lcyclic = .true.
!     Zero: do nothing
      case ('0', '') lower
!     Zeroth-order extrapolation
      case ('cop') lower
        omb1(1) = omb(1) - a(1)
      case ('s') lower
        c1(1) = c(1) + a(1)
      case ('a') lower
        c1(1) = c(1) - a(1)
      case ('sfr') lower
        c1(1) = c(1) + a(1)*(bound(BOT)-d2_bound(-1)/2.)/(bound(BOT)+d2_bound(-1)/2.)
!
! Alternative formulation
!
        !c1(1) = c(1) + a(1)
        !omb1(1) = omb(1) + a(1)*(d2_bound(-1)/bound(BOT))
      case ('nfr') lower
        c1(1) = c(1) + a(1)
        omb1(1) = omb(1) - a(1)*(d2_bound(-1)/bound(BOT))
!     Unknown
      case default lower
        call fatal_error('implicit_pencil', 'generalize this subroutine to deal with your boundary conditions')
      endselect lower
!
! Check the upper boundary condition.
!
      upper: select case (bc12(2))
!     Periodic
      case ('p') upper
        alpha = -c(n)
        lcyclic = .true.
!     Zero: do nothing
      case ('0', '') upper
!     Zeroth-order extrapolation
      case ('cop') upper
        omb1(n) = omb(n) - c(n)
      case ('s') upper
        a1(n) = a(n) + c(n)
      case ('a') upper
        a1(n) = a(n) - c(n)
      case ('sfr') upper
        a1(1) = a(n) + c(n)*(bound(TOP)+d2_bound(1)/2.)/(bound(TOP)-d2_bound(1)/2.)
!
! Alternative formulation
!
        !a1(n) = a(n) + c(n)
        !omb1(n) = omb(n) - c(n)*(d2_bound(1)/bound(TOP))
      case ('nfr') upper
        a1(n) = a(n) + c(n)
        omb1(n) = omb(n) + c(n)*(d2_bound(1)/bound(TOP))
!     Unknown
      case default upper
        call fatal_error('implicit_pencil', 'generalize this subroutine to deal with your boundary conditions')
      endselect upper
!
! Solve the linear system.
!
      if (lcyclic) then
        call cyclic(-a, omb1, -c, alpha, beta, r, q(1:n), n)
      else
        call tridag(-a1, omb1, -c1, r, q(1:n), err, msg)
        if (err) call warning('implicit_pencil', trim(msg))
      endif
!
    endsubroutine implicit_pencil
!***********************************************************************
endmodule ImplicitDiffusion
