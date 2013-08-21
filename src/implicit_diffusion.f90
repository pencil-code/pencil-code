! $Id$
!
!  This module provides general facilities to implicitly solve a
!  diffusion equation.
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
!  Module-specific variables
!
  real, dimension(nxgrid) :: ax_imp = 0.0, bx_imp = 0.0, cx_imp = 0.0
  real, dimension(nygrid) :: ay_imp = 0.0, by_imp = 0.0, cy_imp = 0.0
  real, dimension(nzgrid) :: az_imp = 0.0, bz_imp = 0.0, cz_imp = 0.0
!
  contains
!***********************************************************************
!***********************************************************************
!  PUBLIC ROUTINES GO BELOW HERE.
!***********************************************************************
!***********************************************************************
    subroutine integrate_diffusion(get_diffus_coeff, f, ivar1, ivar2)
!
! Implicitly integrate the diffusion term for components ivar1 to ivar2.
! ivar2 can be omitted if only operating on one component.
!
! 16-aug-13/ccyang: coded.
!
      use Shear, only: sheared_advection_fft
!
      external :: get_diffus_coeff
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: ivar1
      integer, intent(in), optional :: ivar2
!
      integer :: iv1, iv2
      real :: dth
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
!  Get half the time step.
!
      dth = 0.5 * dt
!
!  Integrate.
!
      shearing: if (lshear) then
!       Shear is present: Work in unsheared frame.
        call zsweep(f, bcz1, bcz2, iv1, iv2, dth, get_diffus_coeff)
        call ysweep(f, spread('p', 1, mcom), spread('p', 1, mcom), iv1, iv2, dth, get_diffus_coeff)
!
        call sheared_advection_fft(f, iv1, iv2, real(-t))
        call xsweep(f, spread('p', 1, mcom), spread('p', 1, mcom), iv1, iv2, dth, get_diffus_coeff)
        call xsweep(f, spread('p', 1, mcom), spread('p', 1, mcom), iv1, iv2, dth, get_diffus_coeff)
        call sheared_advection_fft(f, iv1, iv2, real(t))
!
        call ysweep(f, spread('p', 1, mcom), spread('p', 1, mcom), iv1, iv2, dth, get_diffus_coeff)
        call zsweep(f, bcz1, bcz2, iv1, iv2, dth, get_diffus_coeff)
      else shearing
!       No shear is present.
        call xsweep(f, bcx1, bcx2, iv1, iv2, dth, get_diffus_coeff)
        call ysweep(f, bcy1, bcy2, iv1, iv2, dth, get_diffus_coeff)
        call zsweep(f, bcz1, bcz2, iv1, iv2, dth, get_diffus_coeff)
!
        call zsweep(f, bcz1, bcz2, iv1, iv2, dth, get_diffus_coeff)
        call ysweep(f, bcy1, bcy2, iv1, iv2, dth, get_diffus_coeff)
        call xsweep(f, bcx1, bcx2, iv1, iv2, dth, get_diffus_coeff)
      endif shearing
!
    endsubroutine integrate_diffusion
!***********************************************************************
!***********************************************************************
!  LOCAL ROUTINES GO BELOW HERE.
!***********************************************************************
!***********************************************************************
    subroutine set_diffusion_equations(get_diffus_coeff, direction, a, b, c)
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
!   subroutine get_diffus_coeff(diffus_coeff)
!
!     real, dimension(:), intent(out) :: diffus_coeff
!
! where diffus_coeff(i) is the diffusion coefficient at position i.
!
! 16-aug-13/ccyang: coded
!
      external :: get_diffus_coeff
      integer, intent(in) :: direction
      real, dimension(:), intent(out) :: a, b, c
!
      real, dimension(max(nxgrid,nygrid,nzgrid)) :: dc
!
! Find the coefficients.
!
      dir: select case (direction)
!
      case (1) dir
        call get_diffus_coeff(dc(1:nxgrid))
        a(1:nxgrid) = 0.5 * dc(1:nxgrid) * dx1grid * (dx1grid - 0.5 * dxtgrid)
        b(1:nxgrid) = -dc(1:nxgrid) * dx1grid**2
        c(1:nxgrid) = 0.5 * dc(1:nxgrid) * dx1grid * (dx1grid + 0.5 * dxtgrid)
!
      case (2) dir
        call get_diffus_coeff(dc(1:nygrid))
        a(1:nygrid) = 0.5 * dc(1:nygrid) * dy1grid * (dy1grid - 0.5 * dytgrid)
        b(1:nygrid) = -dc(1:nygrid) * dy1grid**2
        c(1:nygrid) = 0.5 * dc(1:nygrid) * dy1grid * (dy1grid + 0.5 * dytgrid)
!
      case (3) dir
        call get_diffus_coeff(dc(1:nzgrid))
        a(1:nzgrid) = 0.5 * dc(1:nzgrid) * dz1grid * (dz1grid - 0.5 * dztgrid)
        b(1:nzgrid) = -dc(1:nzgrid) * dz1grid**2
        c(1:nzgrid) = 0.5 * dc(1:nzgrid) * dz1grid * (dz1grid + 0.5 * dztgrid)
!
      case default dir
        call fatal_error('set_up_equations', 'unknown direction')
!
      endselect dir
!
    endsubroutine set_diffusion_equations
!***********************************************************************
    subroutine xsweep(f, bcx1, bcx2, ivar1, ivar2, dt, get_diffus_coeff)
!
! Implicitly integrate the diffusion term in the x-direction.
!
! 16-aug-13/ccyang: coded
!
! Input Arguments
!   bcx1: array of the lower boundary conditions for each component
!   bcx2: array of the upper boundary conditions for each component
!   ivar1: beginning f component
!   ivar2: ending f component
!   dt: time step
! Input Procedure
!   subroutine get_diffus_coeff(diffus_coeff): see set_diffusion_equations
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      character(len=*), dimension(mcom), intent(in) :: bcx1, bcx2
      integer, intent(in) :: ivar1, ivar2
      real, intent(in) :: dt
      external :: get_diffus_coeff
!
      real, dimension(nxgrid) :: a, b, c
      real, dimension(nxgrid) :: adt, opbdt, ombdt, cdt
      integer :: j, k, l
!
      int_x: if (nxgrid > 1) then
        call set_diffusion_equations(get_diffus_coeff, 1, a, b, c)
        call get_tridiag(a, b, c, dt, adt, opbdt, ombdt, cdt)
        zscan: do k = n1, n2
          yscan: do j = m1, m2
            do l = ivar1, ivar2
              call implicit_pencil(f(l1-1:l2+1,j,k,l), nxgrid, adt, opbdt, ombdt, cdt, bcx1(l), bcx2(l))
            enddo
          enddo yscan
        enddo zscan
      endif int_x
!
    endsubroutine xsweep
!***********************************************************************
    subroutine ysweep(f, bcy1, bcy2, ivar1, ivar2, dt, get_diffus_coeff)
!
! Implicitly integrate the diffusion term in the y-direction.
!
! 16-aug-13/ccyang: coded
!
! Input Arguments
!   bcy1: array of the lower boundary conditions for each component
!   bcy2: array of the upper boundary conditions for each component
!   ivar1: beginning f component
!   ivar2: ending f component
!   dt: time step
! Input Procedure
!   subroutine get_diffus_coeff(diffus_coeff): see set_diffusion_equations
!
      use Mpicomm, only: transp_xy
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      character(len=*), dimension(mcom), intent(in) :: bcy1, bcy2
      integer, intent(in) :: ivar1, ivar2
      real, intent(in) :: dt
      external :: get_diffus_coeff
!
      real, dimension(nx,ny) :: axy      ! assuming nxgrid = nygrid and nprocx = 1
      real, dimension(0:nx+1) :: penc    ! assuming nxgrid = nygrid and nprocx = 1
      real, dimension(nygrid) :: a, b, c
      real, dimension(nygrid) :: adt, opbdt, ombdt, cdt
      integer :: j, k, l
!
      int_y: if (nygrid > 1) then
        call set_diffusion_equations(get_diffus_coeff, 2, a, b, c)
        call get_tridiag(a, b, c, dt, adt, opbdt, ombdt, cdt)
        comp: do l = ivar1, ivar2
          zscan: do k = n1, n2
            axy = f(l1:l2,m1:m2,k,l)
            call transp_xy(axy)  ! assuming nxgrid = nygrid
            xscan: do j = 1, ny
              penc(1:nx) = axy(:,j)
              call implicit_pencil(penc, nygrid, adt, opbdt, ombdt, cdt, bcy1(l), bcy2(l))
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
    subroutine zsweep(f, bcz1, bcz2, ivar1, ivar2, dt, get_diffus_coeff)
!
! Implicitly integrate the diffusion term in the z-direction.
!
! 16-aug-13/ccyang: coded
!
! Input Arguments
!   bcz1: array of the lower boundary conditions for each component
!   bcz2: array of the upper boundary conditions for each component
!   ivar1: beginning f component
!   ivar2: ending f component
!   dt: time step
! Input Procedure
!   subroutine get_diffus_coeff(diffus_coeff): see set_diffusion_equations
!
      use Mpicomm, only: transp_xz, transp_zx
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      character(len=*), dimension(mcom), intent(in) :: bcz1, bcz2
      integer, intent(in) :: ivar1, ivar2
      real, intent(in) :: dt
      external :: get_diffus_coeff
!
      integer, parameter :: nxt = nx / nprocz
      real, dimension(nx,nz) :: axz
      real, dimension(nzgrid,nxt) :: azx
      real, dimension(0:nzgrid+1) :: penc
      real, dimension(nzgrid) :: a, b, c
      real, dimension(nzgrid) :: adt, opbdt, ombdt, cdt
      integer :: j, k, l
!
      int_z: if (nzgrid > 1) then
        call set_diffusion_equations(get_diffus_coeff, 3, a, b, c)
        call get_tridiag(a, b, c, dt, adt, opbdt, ombdt, cdt)
        comp: do l = ivar1, ivar2
          yscan: do j = m1, m2
            axz = f(l1:l2,j,n1:n2,l)
            call transp_xz(axz, azx)
            xscan: do k = 1, nxt
              penc(1:nzgrid) = azx(:,k)
              call implicit_pencil(penc, nzgrid, adt, opbdt, ombdt, cdt, bcz1(l), bcz2(l))
              azx(:,k) = penc(1:nzgrid)
            enddo xscan
            call transp_zx(azx, axz)
            f(l1:l2,j,n1:n2,l) = axz
          enddo yscan
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
!   adt: array of dt * a_i
!   opbdt: array of 1 + dt * b_i
!   ombdt: array of 1 - dt * b_i
!   cdt: array of dt * c_i
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
    subroutine implicit_pencil(q, n, a, opb, omb, c, bc1, bc2)
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
!
      use Boundcond, only: bc_pencil
      use General, only: cyclic, tridag
!
      integer, intent(in) :: n
      real, dimension(0:n+1), intent(inout) :: q
      real, dimension(n), intent(in) :: a, opb, omb, c
      character(len=*), intent(in) :: bc1, bc2
!
      real, dimension(n) :: r, omb1
      character(len=255) :: msg
      logical :: lcyclic, err
      real :: alpha, beta
!
! Prepare the RHS of the linear system.
!
      call bc_pencil(q, n, 1, bc1, bc2)
      r = a * q(0:n-1) + opb * q(1:n) + c * q(2:n+1)
!
! Assume non-periodic boundary conditions first.
!
      lcyclic = .false.
      alpha = 0.0
      beta = 0.0
      omb1 = omb
!
! Check the lower boundary conditions.
!
      lower: select case (bc1)
!     Periodic
      case ('p') lower
        beta = -a(1)
        lcyclic = .true.
!     Zero: do nothing
      case ('0', '') lower
!     Zeroth-order extrapolation
      case ('cop') lower
        omb1(1) = omb(1) - a(1)
!     Unknown
      case default lower
        call fatal_error('implicit_pencil', 'generalize this subroutine to deal with your boundary conditions')
      endselect lower
!
! Check the upper boundary condition.
!
      upper: select case (bc2)
!     Periodic
      case ('p') upper
        alpha = -c(n)
        lcyclic = .true.
!     Zero: do nothing
      case ('0', '') upper
!     Zeroth-order extrapolation
      case ('cop') upper
        omb1(n) = omb(n) - c(n)
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
        call tridag(-a, omb1, -c, r, q(1:n), err, msg)
        if (err) call warning('implicit_pencil', trim(msg))
      endif
!
    endsubroutine implicit_pencil
!***********************************************************************
endmodule ImplicitDiffusion
