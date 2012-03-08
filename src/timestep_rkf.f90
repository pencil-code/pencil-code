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
  public :: time_step
!
! Parameters for adaptive time stepping
  real, parameter :: safety      =  0.9
  real, parameter :: dt_decrease = -0.25
  real, parameter :: dt_increase = -0.20
  real            :: errcon
!
  contains
!
!***********************************************************************
    subroutine time_step(f,df,p)
!
!  Cash-Karp variant of Runge-Kutta-Fehlberg accurate to 5th order
!  To use this, set itorder to 5.
!
!  22-jun-06/tony: coded
!
      use Messages, only: fatal_error
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real :: errmax, tnew
      real :: dt_temp, dt_next
      integer :: j,i
!
      ldt=.false.
!
      ! General error condition
      errcon = (5.0/safety)**(1.0/dt_increase)
!
      if (itorder/=5) &
        call fatal_error('time_step','itorder must be 5 for Runge-Kutta-Fehlberg')
!
!
!  dt_beta_ts may be needed in other modules (like Dustdensity) for fixed dt
!
!      if (.not. ldt) dt_beta_ts=dt*beta_ts
!
!
      if (linterstellar .or. lshear .or. lparticles) &
            call fatal_error("time_step", &
                   "Shear, interstellar and particles are not" // &
                   " yet supported by the adaptive rkf scheme")
!
      lfirst=.true.
      do i=1,10
        ! Do a Runge-Kutta step
        call rkck(f, df, p, errmax)
        ! Step succeeded so exit
        if (errmax <= safety) exit
        ! Step didn't succeed so decrease the time step
!        print*,"Decreasing",errmax,dt
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
      if (errmax > errcon) then
        ! Increase the time step
        dt_next = safety*dt*(errmax**dt_increase)
      else
        ! But not by more than a factor of 5
        dt_next = 5.0*dt
      endif
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
    endsubroutine time_step
!***********************************************************************
    subroutine rkck(f, df, p, errmax)
!
! Explicit fifth order Runge--Kutta--Fehlberg time stepping
!
      use Mpicomm, only: mpiallreduce_max
      use Equ, only: pde
!
! RK parameters by Cash and Karp
!
      real, parameter :: b21 = 0.2
      real, parameter :: b31 = 0.075
      real, parameter :: b32 = 0.225
      real, parameter :: b41 = 0.3
      real, parameter :: b42 = -0.9
      real, parameter :: b43 = 1.2
      real, parameter :: b51 = -11.0 / 54.0
      real, parameter :: b52 = 2.5
      real, parameter :: b53 = -70.0 / 27.0
      real, parameter :: b54 = 35.0 / 27.0
      real, parameter :: b61 = 1631.0 / 55296.0
      real, parameter :: b62 = 175.0 / 512.0
      real, parameter :: b63 = 575.0 / 13824.0
      real, parameter :: b64 = 44275.0 / 110592.0
      real, parameter :: b65 = 253.0 / 4096.0
      real, parameter :: c1  = 37.0 / 378.0
      real, parameter :: c2  = 0.0
      real, parameter :: c3  = 250.0 / 621.0
      real, parameter :: c4  = 125.0 / 594.0
      real, parameter :: c5  = 0.0
      real, parameter :: c6  = 512.0 / 1771.0
      real, parameter :: dc1 = c1 - 2825.0 / 27648.0
      real, parameter :: dc2 = c2 - 0.0
      real, parameter :: dc3 = c3 - 18575.0 / 48384.0
      real, parameter :: dc4 = c4 - 13525.0 / 55296.0
      real, parameter :: dc5 = c5 - 277.0 / 14336.0
      real, parameter :: dc6 = c6 - 0.25
!
! Note: f is intent(inout), as pde may change it due to
!   lshift_datacube_x, density floor, or velocity ceiling.
!   None of those will do exctly what is intended, because they are
!   only really modifying f during the first substep.
!
      intent(inout) :: f
      intent(out)   :: df, p, errmax
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension(mx,my,mz,mfarray,5) :: k
! Note: The tmp array will not use more memory than the temporary
!   array that would be implicitly created with calls like
!   call pde(f + b21*k(:,:,:,:,1), k(:,:,:,:,2), p)
      real, dimension (mx,my,mz,mfarray) :: tmp
      real, dimension(nx) :: scal, err
      real :: errmax, errmaxs
      integer :: j
!
      df=0.
      errmax=0.
      errmaxs=0.
      k=0.
!
      call pde(f, k(:,:,:,:,1), p)
      do j=1,mvar; do n=n1,n2; do m=m1,m2
        k(l1:l2,m,n,j,1) = dt*k(l1:l2,m,n,j,1)
      enddo; enddo; enddo
!
      lfirst=.false.
      tmp = f + b21*k(:,:,:,:,1)
!
      call pde(tmp, k(:,:,:,:,2), p)
      do j=1,mvar; do n=n1,n2; do m=m1,m2
        k(l1:l2,m,n,j,2) = dt*k(l1:l2,m,n,j,2)
      enddo; enddo; enddo
!
      tmp = f + b31*k(:,:,:,:,1) &
              + b32*k(:,:,:,:,2)
!
      call pde(tmp, k(:,:,:,:,3), p)
      do j=1,mvar; do n=n1,n2; do m=m1,m2
        k(l1:l2,m,n,j,3) = dt*k(l1:l2,m,n,j,3)
      enddo; enddo; enddo
!
      tmp = f + b41*k(:,:,:,:,1) &
              + b42*k(:,:,:,:,2) &
              + b43*k(:,:,:,:,3)
!
      call pde(tmp, k(:,:,:,:,4), p)
      do j=1,mvar; do n=n1,n2; do m=m1,m2
        k(l1:l2,m,n,j,4) = dt*k(l1:l2,m,n,j,4)
      enddo; enddo; enddo
!
      tmp = f + b51*k(:,:,:,:,1) &
              + b52*k(:,:,:,:,2) &
              + b53*k(:,:,:,:,3) &
              + b54*k(:,:,:,:,4)
!
      call pde(tmp, k(:,:,:,:,5), p)
      do j=1,mvar; do n=n1,n2; do m=m1,m2
        k(l1:l2,m,n,j,5) = dt*k(l1:l2,m,n,j,5)
      enddo; enddo; enddo
!
      tmp = f + b61*k(:,:,:,:,1) &
              + b62*k(:,:,:,:,2) &
              + b63*k(:,:,:,:,3) &
              + b64*k(:,:,:,:,4) &
              + b65*k(:,:,:,:,5)
!
      call pde(tmp, df, p)
!
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
          ! Per variable error
          !          scal=  ( &
          !              sqrt(f(l1:l2,m,n,1)**2+f(l1:l2,m,n,2)**2)  + &
          !              sqrt(k(l1:l2,m,n,1,1)**2 + k(l1:l2,m,n,2,1)**2) + 1e-30)
          !          errmaxs = max(maxval(abs(err/scal)),errmaxs)
          scal=  sqrt(f(l1:l2,m,n,j)**2  + k(l1:l2,m,n,j,1)**2 + 1e-30)
          errmaxs = max(maxval(abs(err/scal)),errmaxs)
        case ('cons_frac_err')
          ! Constant fractional error
          errmaxs = max(maxval(abs(err/f(l1:l2,m,n,j))),errmaxs)
        case ('cons_err')
          ! Constant error
          scal = max(abs(f(l1:l2,m,n,j)), 1e-8)
          errmaxs = max(maxval(abs(err/scal)),errmaxs)
          !
        case ('none')
          ! No error check
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
