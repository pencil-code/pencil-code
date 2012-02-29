! $Id$
!
module Timestep
!
  use Cparam
  use Cdata
!
  implicit none
!
  private
!
  public :: rk_2n
!
  ! Parameters for adaptive time stepping
  real, parameter :: safety           = 0.9
  real, parameter :: dt_decrease      = -0.25
  real, parameter :: dt_increase      = -0.20
  real            :: errcon
!
  contains
!
!***********************************************************************
    subroutine rk_2n(f,df,p)
!
!  Runge-Kutta-Fehlberg accurate to 5th order
!
!  22-jun-06/tony: coded
!
      use Mpicomm
      use Cdata
      use Messages
!!      use Particles_main
!!      use Interstellar, only: calc_snr_damp_int
!!      use Shear, only: advance_shear
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real :: errmax, tnew
      real :: dt_temp, dt_next, dt_did
      integer :: j,i
!
      ldt=.false.
!
      ! General error condition
      errcon = (5.0/safety)**(1.0/dt_increase)
!
      if (itorder/=5) &
        call fatal_error('rk_2n','itorder must be 5 for Runge-Kutta-Fehlberg')
!
!
!  dt_beta_ts may be needed in other modules (like Dustdensity) for fixed dt
!
!      if (.not. ldt) dt_beta_ts=dt*beta_ts
      !
!
      if (linterstellar.or.lshear.or.lparticles) &
            call fatal_error("rk_2n", &
                   "Shear, interstallar and particles are not" // &
                   " yet supported by the adaptive rkf scheme")
!
      lfirst=.true.
      do i=1,10
        ! Do a Runge-Kutta step
        call rkck(f, df, p, errmax)
        ! Step succeeded so exit
        if (errmax <= safety) exit
        ! Step didn't succeed so decrease the time step
!        print*,"Decreasing"
        dt_temp = safety*dt*(errmax**dt_decrease)
        ! Don't decrease the time step by more than a factor of ten
        dt = sign(max(abs(dt_temp), 0.1*abs(dt)), dt)
        ! New time
        tnew = t+dt
        if (tnew == t) then
          ! Guard against infinitesimal time steps
          print*, 'WARNING: Timestep underflow in rkqs()'
        endif
      enddo
!
!      print*,"errmax, errcon", errmax,errcon
      if (errmax > errcon) then
        ! Increase the time step
        dt_next = safety*dt*(errmax**dt_increase)
      else
        ! But not by more than a factor of 5
        dt_next = 5.0*dt
      endif
!
      ! Time step that was actually performed
      dt_did = dt
!
      if (ip<=6) print*,'TIMESTEP: iproc,dt=',iproc,dt  !(all have same dt?)
      ! Increase time
      t = t+dt
      ! Time step to try next time
      dt = dt_next
!
!  Time evolution of grid variables
!  (do this loop in pencils, for cache efficiency)
!
      do j=1,mvar; do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,j)=f(l1:l2,m,n,j)+df(l1:l2,m,n,j)
      enddo; enddo; enddo
!
    endsubroutine rk_2n
!***********************************************************************
    subroutine rkck(f, df, p, errmax)
    ! Explicit fifth order Runge--Kutta--Fehlberg time stepping
      use Cdata
      use Mpicomm, only: mpiallreduce_max
      use Messages
      use Equ
    ! RK parameters by Cash and Karp
      real, parameter :: b21      = 0.2
      real, parameter :: b31      = 0.075
      real, parameter :: b32      = 0.225
      real, parameter :: b41      = 0.3
      real, parameter :: b42      = -0.9
      real, parameter :: b43      = 1.2
      real, parameter :: b51      = -11.0 / 54.0
      real, parameter :: b52      = 2.5
      real, parameter :: b53      = -70.0 / 27.0
      real, parameter :: b54      = 35.0 / 27.0
      real, parameter :: b61      = 1631.0 / 55296.0
      real, parameter :: b62      = 175.0 / 512.0
      real, parameter :: b63      = 575.0 / 13824.0
      real, parameter :: b64      = 44275.0 / 110592.0
      real, parameter :: b65      = 253.0 / 4096.0
      real, parameter :: c1       = 37.0 / 378.0
      real, parameter :: c2       = 0.0
      real, parameter :: c3       = 250.0 / 621.0
      real, parameter :: c4       = 125.0 / 594.0
      real, parameter :: c5       = 0.0
      real, parameter :: c6       = 512.0 / 1771.0
      real, parameter :: dc1      = c1 - 2825.0 / 27648.0
      real, parameter :: dc2      = c2 - 0.0
      real, parameter :: dc3      = c3 - 18575.0 / 48384.0
      real, parameter :: dc4      = c4 - 13525.0 / 55296.0
      real, parameter :: dc5      = c5 - 277.0 / 14336.0
      real, parameter :: dc6      = c6 - 0.25
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,mvar), intent(out) :: df
      type (pencil_case), intent(inout) :: p
      real, dimension(mx,my,mz,mfarray,5) :: k
      real, save, dimension(mx,my,mz,mfarray) :: tmp1
      real, dimension(mx,my,mz,mfarray) :: tmp2
      real, dimension(nx) :: scal, err
      real, intent(inout) :: errmax
      real :: errmaxs
      integer :: j
      logical, save :: first_call=.true.
!
      if (ny /= 1 .or. nz /= 1) then
        call fatal_error("rkck", "timestep_rkf_1d only works for the 1D case")
      endif
!
      errmax=0.
!
      if (first_call) then
        ! Initialize tmp1 to arbitrary value /= 0, so chemistry.f90 (which
        ! operates in the ghost zones where it most probably shouldn't)
        ! doesn't divide by 0.
        tmp1 = real(0.577215664901532860606512090082402431042159335)
        first_call = .false.
      endif
!
      tmp2(:,m1:m2,n1:n2,:) = 0.
      call pde(f,tmp2,p)
      k(:,m1:m2,n1:n2,:,1) = tmp2(:,m1:m2,n1:n2,:)
      do j=1,mvar; do n=n1,n2; do m=m1,m2
          k(l1:l2,m,n,j,1) = dt*k(l1:l2,m,n,j,1)
!
      enddo; enddo; enddo
!
      lfirst=.false.
!
      tmp1(:,m1:m2,n1:n2,:) = f(:,m1:m2,n1:n2,:) &
          + b21*k(:,m1:m2,n1:n2,:,1)
      tmp2(:,m1:m2,n1:n2,:) = 0.
      call pde(tmp1,tmp2,p)
      k(:,m1:m2,n1:n2,:,2) = tmp2(:,m1:m2,n1:n2,:)
      do j=1,mvar; do n=n1,n2; do m=m1,m2
          k(l1:l2,m,n,j,2) = dt*k(l1:l2,m,n,j,2)
      enddo; enddo; enddo
!
      tmp1(:,m1:m2,n1:n2,:) = f(:,m1:m2,n1:n2,:) &
          + b31*k(:,m1:m2,n1:n2,:,1) &
          + b32*k(:,m1:m2,n1:n2,:,2)
      tmp2(:,m1:m2,n1:n2,:) = 0.
      call pde(tmp1,tmp2,p)
      k(:,m1:m2,n1:n2,:,3) = tmp2(:,m1:m2,n1:n2,:)
      do j=1,mvar; do n=n1,n2; do m=m1,m2
          k(l1:l2,m,n,j,3) = dt*k(l1:l2,m,n,j,3)
      enddo; enddo; enddo
!
      tmp1(:,m1:m2,n1:n2,:) = f(:,m1:m2,n1:n2,:) &
          + b41*k(:,m1:m2,n1:n2,:,1) &
          + b42*k(:,m1:m2,n1:n2,:,2) &
          + b43*k(:,m1:m2,n1:n2,:,3)
      tmp2(:,m1:m2,n1:n2,:) = 0.
      call pde(tmp1,tmp2,p)
      k(:,m1:m2,n1:n2,:,4) = tmp2(:,m1:m2,n1:n2,:)
      do j=1,mvar; do n=n1,n2; do m=m1,m2
          k(l1:l2,m,n,j,4) = dt*k(l1:l2,m,n,j,4)
      enddo; enddo; enddo
!
      tmp1(:,m1:m2,n1:n2,:) = f(:,m1:m2,n1:n2,:) &
          + b51*k(:,m1:m2,n1:n2,:,1) &
          + b52*k(:,m1:m2,n1:n2,:,2) &
          + b53*k(:,m1:m2,n1:n2,:,3) &
          + b54*k(:,m1:m2,n1:n2,:,4)
      tmp2(:,m1:m2,n1:n2,:) = 0.
      call pde(tmp1,tmp2,p)
      k(:,m1:m2,n1:n2,:,5) = tmp2(:,m1:m2,n1:n2,:)
      do j=1,mvar; do n=n1,n2; do m=m1,m2
          k(l1:l2,m,n,j,5) = dt*k(l1:l2,m,n,j,5)
      enddo; enddo; enddo
!
      errmaxs=0.
      tmp1(:,m1:m2,n1:n2,:) = f(:,m1:m2,n1:n2,:) &
          + b61*k(:,m1:m2,n1:n2,:,1) &
          + b62*k(:,m1:m2,n1:n2,:,2) &
          + b63*k(:,m1:m2,n1:n2,:,3) &
          + b64*k(:,m1:m2,n1:n2,:,4) &
          + b65*k(:,m1:m2,n1:n2,:,5)
      df(:,m1:m2,n1:n2,:)=0.
      call pde(tmp1,df,p)
      do j=1,mvar; do n=n1,n2; do m=m1,m2
          df(l1:l2,m,n,j) = dt*df(l1:l2,m,n,j)
!
          err = dc1*k(l1:l2,m,n,j,1) + dc2*k(l1:l2,m,n,j,2) + &
                dc3*k(l1:l2,m,n,j,3) + dc4*k(l1:l2,m,n,j,4) + &
                dc5*k(l1:l2,m,n,j,5) + dc6*df(l1:l2,m,n,j)
!
          df(l1:l2,m,n,j) = c1*k(l1:l2,m,n,j,1) + c2*k(l1:l2,m,n,j,2) + &
                            c3*k(l1:l2,m,n,j,3) + c4*k(l1:l2,m,n,j,4) + &
                            c5*k(l1:l2,m,n,j,5) + c6*df(l1:l2,m,n,j)
!
          ! Get the maximum error over the whole field
          !
          select case (timestep_scaling(j))
          case ('per_var_err')
            !
            ! Per variable error
            !
            scal=  ( &
                 sqrt(f(l1:l2,m,n,1)**2+f(l1:l2,m,n,2)**2)  + &
                 sqrt(k(l1:l2,m,n,1,1)**2 + k(l1:l2,m,n,2,1)**2) + &
                 1e-30)
            errmaxs = max(maxval(abs(err/scal)),errmaxs)
            !scal=  ( &
            !     abs(f(l1:l2,m,n,j))  + abs(k(l1:l2,m,n,j,1)) + 1e-30)
            !errmaxs = max(maxval(abs(err/scal)),errmaxs)
          case ('cons_frac_err')
            !
            ! Constant fractional error
            !
            errmaxs = max(maxval(abs(err/f(l1:l2,m,n,j))),errmaxs)
          case ('cons_err')
            !
            ! Constant error
            !
            scal = max(abs(f(l1:l2,n,m,j)), 1e-8)
            errmaxs = max(maxval(abs(err/scal)),errmaxs)
            !
          case ('none')
            !
            ! No error check
            !
            errmaxs = 0
            !
          endselect
          !
        enddo; enddo; enddo
        !
        ! Divide your maximum error by the required accuracy
        !
        errmaxs=errmaxs/eps_rkf
        !
      call mpiallreduce_max(errmaxs,errmax)
!
    endsubroutine rkck
!***********************************************************************
endmodule Timestep
