! $Id$
!
module Timestep
!
  use Cdata
!
  implicit none
!
  private
!
  include 'timestep.h'
!
! Parameters for adaptive time stepping
  real, parameter :: safety      =  0.95
  real            :: errcon, dt_next, dt_increase, dt_decrease
  real, dimension(mvar) :: farraymin
!
  contains
!
!***********************************************************************
    subroutine initialize_timestep
!
      use Messages, only: fatal_error, warning
      use General, only: rtoa
!
      if (lparticles) call fatal_error("initialize_timestep", "Particles are"// &
                                       " not yet supported by the adaptive rkf scheme")
      if (itorder/=5.and.itorder/=3) then
        call warning('initialize_timestep','itorder set to 5 or 3 for Runge-Kutta-Fehlberg')
        itorder=5
      else if (itorder==3) then
        call warning('initialize_timestep',&
                     'Runge-Kutta-Fehlberg itorder is 3: set to 5 for higher accuracy')
      endif
!
      ldt = (dt==0.)
      if (ldt) then
        if (dt0==0.) then
          dt = dt_epsi
        else
          dt = dt0
        endif
      endif
      lcourant_dt=.false.
      dt0 = 0.
!
      if (eps_rkf0/=0.) eps_rkf=eps_rkf0
!
      dt_next=dt
      dt_increase=-1./(itorder+dtinc)
      dt_decrease=-1./(itorder-dtdec)
      num_substeps = itorder
!
    endsubroutine initialize_timestep
!***********************************************************************
    subroutine time_step(f,df,p)
!
!  Cash-Karp variant of Runge-Kutta-Fehlberg accurate to 5th order
!  To use this, set itorder to 5.
!
!  22-jun-06/tony: coded
!  08-feb-24/fred: revisions based on high resolution simulations of SN-driven turbulence
!     notes: Courant time is insufficient to safeguard high Mach number turbulence with
!     significant sources and sinks, and strong viscous stresses that are beyond the Courant
!     analysis, hence tried RKF. Fifth order produces large time step than previous ISM
!     methods used but the increased iterations over the pde results in longer total
!     integration time. 3rd order has fewer overheads than default CFT method and some savings in
!     redundant timestep calculations, providing the number of iterations of rkck can be
!     minimised. Cash-Karp 5th order method uses dt_increase=-1/5 and dt_devrease=-1/4.
!     dtinc and dtdec optimal around 0.5 for 3rd order scheme, but resolution sensitve, so worth
!     testing on a new physical setup to optimise algorithm. Cash-Karp saftey=0.9 revised here to
!     0.95 as 0.9 overshoots, saftey<1 for dt_temp can actually reduce the timestep, so now omitted.
!     "cons_frac_err" for relative error normalisation only method verified and sensitive to choices
!     dt_epsi and dt_ratio. Sensitivity to eps_rkf very nonlinear and resolution dependent.
!  14-oct-24/fred: "rel_err" now reliable option for accuracy constraint.
!                  "lreiterate=F": implement solution with current dt and instead only adjust dt
!                  for the next time step. Timestep very similar and stable without reiterating.
!
      use Messages, only: warning

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real :: errmax, dt_temp, dt_last
      real(KIND=rkind8) :: tnew, told
      integer :: j,i
!
!  dt_ratio is lower bound on a processor for the denominator of each variable
!  dt_epsi is the lower bound for the denominator on any variable or processor
!
      do j=1,mvar
        farraymin(j) = max(dt_ratio*maxval(abs(f(l1:l2,m1:m2,n1:n2,j))),dt_epsi)
      enddo
      if (lroot.and.it==1) print*,"farraymin",farraymin
!
      lfirst=.true.
      told=t
      dt_last=dt
      dt=dt_next
      do i=1,20
        ! Do a Runge-Kutta step
        if (itorder==5) then
          call rkck(f, df, p, errmax)
        else
          call rkck3(f, df, p, errmax)
        endif
        ! Step succeeded so exit
        if (errmax <= 1.or..not.ldt) exit
        ! If f not to be stored for reiteration errmax constraint below must be removed TBA
        if (.not.lreiterate.and.errmax<=1.75) exit
        ! Step didn't succeed so decrease the time step
        dt_temp = safety*dt*errmax**dt_decrease
        ! Don't decrease the time step by more than a factor of ten
        dt = sign(max(abs(dt_temp), 0.1*abs(dt)), dt)
        if (lroot.and.ip==6787) print*,"time_step: dt",dt,"to dt_temp",dt_temp,"at errmax",errmax
        exit
        tnew=told+dt
        if (tnew == told) then
          ! Guard against infinitesimal time steps
          call warning('time_step','Timestep underflow in Runge-Kutta-Fehlberg')
          ! Increase dt???
          exit
        endif
        t=told
      enddo
!
!  The algorithm is optimised if the number of iterations is mainly 1 and occasionally 2
!  errmax should not be much less than 1, otherwise dt_next is more likely to overshoot
!  and there will be more iterations. The ratio of dt/dt_next should not be small if dt_next
!  is a reasonable try.
!
      if (lroot.and.ip==7.and.i>1) then
        print*,"time_step: rkck",i,"iterations"
        print*,"time_step: rkck",errmax,"converged errmax"
        print*,"time_step: rkck",dt_next,"tried dt_next"
        print*,"time_step: rkck",dt,"used dt"
        print*,"time_step: rkck",dt/dt_next,"ratio reduction"
      endif
      call update_after_substep(f,df,dt,.true.)
!
! Time step to try next time
!
      if (ldt) then
        if (lreiterate) then
          dt_next = dt*errmax**dt_increase
        else
          if (errmax <= 1) then
            dt_next = dt*errmax**dt_increase
          else
            dt_temp = safety*dt*errmax**dt_decrease
            ! Don't decrease the time step by more than a factor of ten
            dt_next = sign(max(abs(dt_temp), 0.1*abs(dt)), dt)
          endif
        endif
      endif
!
      if (ip<=6) print*,'TIMESTEP: iproc,dt=',iproc_world,dt
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
      use Mpicomm, only: mpiallreduce_max,MPI_COMM_PENCIL
      use Equ, only: pde, impose_floors_ceilings
      use Shear, only: advance_shear
      use Boundcond, only: update_ghosts
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
      real, parameter :: substep1 = 0.2
      real, parameter :: substep2 = 0.3
      real, parameter :: substep3 = 0.6
      real, parameter :: substep4 = 1.0
      real, parameter :: substep5 = 0.875
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
      real, dimension(mx,my,mz,mvar,5) :: k
! Note: The tmp array will not use more memory than the temporary
!   array that would be implicitly created with calls like
!   call pde(f + b21*k(:,:,:,:,1), k(:,:,:,:,2), p)
      real, dimension (mx,my,mz,mfarray) :: tmp
      real, dimension(nx) :: scal, err
      real :: errmax, errmaxs, dtsub, told
      integer :: j
!
      df=0.
      errmax=0.
      errmaxs=0.
      k=0.
      told=t
!
!  FIRST SUBSTEP
!
      call pde(f, k(:,:,:,:,1), p)
      dtsub = dt*substep1
      if (lshear) then
        tmp = f
        call impose_floors_ceilings(tmp)
        call update_ghosts(tmp)  ! Necessary for non-FFT advection but unnecessarily overloading FFT advection
        call advance_shear(tmp, k(:,:,:,:,1), dtsub)
      endif
      do j=1,mvar; do n=n1,n2; do m=m1,m2
        k(l1:l2,m,n,j,1) = dt*k(l1:l2,m,n,j,1)
      enddo; enddo; enddo
!
      lfirst=.false.
!
!  Transferring auxiliaries from f to tmp.
!
      if (mfarray>mvar) tmp(:,:,:,mvar+1:mfarray) = f(:,:,:,mvar+1:mfarray)
      tmp(:,:,:,1:mvar) = f(:,:,:,1:mvar) + b21*k(:,:,:,:,1)
!
      call update_after_substep(tmp,k(:,:,:,:,1),dtsub,.false.)
      t = told + dtsub
!
!  SECOND SUBSTEP
!
      call pde(tmp, k(:,:,:,:,2), p)
      dtsub = dt*substep2
      do j=1,mvar; do n=n1,n2; do m=m1,m2
        k(l1:l2,m,n,j,2) = dt*k(l1:l2,m,n,j,2)
      enddo; enddo; enddo
!
      tmp(:,:,:,1:mvar) = f(:,:,:,1:mvar) + b31*k(:,:,:,:,1) &
                                          + b32*k(:,:,:,:,2)
      call update_after_substep(tmp,k(:,:,:,:,2),0.,.false.)
      t = told + dtsub
!
!
!  THIRD SUBSTEP
!
      call pde(tmp, k(:,:,:,:,3), p)
      dtsub = dt*substep3
      if (lshear) then
        call impose_floors_ceilings(tmp)
        call update_ghosts(tmp)
        call advance_shear(tmp, k(:,:,:,:,3), dtsub)
      endif
      do j=1,mvar; do n=n1,n2; do m=m1,m2
        k(l1:l2,m,n,j,3) = dt*k(l1:l2,m,n,j,3)
      enddo; enddo; enddo
!
      tmp(:,:,:,1:mvar) = f(:,:,:,1:mvar) + b41*k(:,:,:,:,1) &
                                          + b42*k(:,:,:,:,2) &
                                          + b43*k(:,:,:,:,3)
!
      call update_after_substep(tmp,k(:,:,:,:,3),dtsub,.false.)
      t = told + dtsub
!
!  FOURTH SUBSTEP
!
      call pde(tmp, k(:,:,:,:,4), p)
      dtsub = dt*substep4
      if (lshear) then
        call impose_floors_ceilings(tmp)
        call update_ghosts(tmp)
        call advance_shear(tmp, k(:,:,:,:,4), dtsub)
      endif
      do j=1,mvar; do n=n1,n2; do m=m1,m2
        k(l1:l2,m,n,j,4) = dt*k(l1:l2,m,n,j,4)
      enddo; enddo; enddo
!
      tmp(:,:,:,1:mvar) = f(:,:,:,1:mvar) + b51*k(:,:,:,:,1) &
                                          + b52*k(:,:,:,:,2) &
                                          + b53*k(:,:,:,:,3) &
                                          + b54*k(:,:,:,:,4)
!
      call update_after_substep(tmp,k(:,:,:,:,4),dtsub,.false.)
      t = told + dtsub
!
!  FIFTH SUBSTEP
!
      call pde(tmp, k(:,:,:,:,5), p)
      dtsub = dt*substep5
      if (lshear) then
        call impose_floors_ceilings(tmp)
        call update_ghosts(tmp)
        call advance_shear(tmp, k(:,:,:,:,5), dtsub)
      endif
      do j=1,mvar; do n=n1,n2; do m=m1,m2
        k(l1:l2,m,n,j,5) = dt*k(l1:l2,m,n,j,5)
      enddo; enddo; enddo
!
      tmp(:,:,:,1:mvar) = f(:,:,:,1:mvar) + b61*k(:,:,:,:,1) &
                                          + b62*k(:,:,:,:,2) &
                                          + b63*k(:,:,:,:,3) &
                                          + b64*k(:,:,:,:,4) &
                                          + b65*k(:,:,:,:,5)
!
      call update_after_substep(tmp,k(:,:,:,:,5),dtsub,.false.)
      t = told + dtsub
!
!  FINALIZE STEP
!
      call pde(tmp, df, p)
      t = told + dt
!
      do j=1,mvar; do n=n1,n2; do m=m1,m2
        df(l1:l2,m,n,j) = dt*df(l1:l2,m,n,j)
!
        err = dc1*k(l1:l2,m,n,j,1) + &!dc2*k(l1:l2,m,n,j,2) + &
              dc3*k(l1:l2,m,n,j,3) + dc4*k(l1:l2,m,n,j,4) + &
              dc5*k(l1:l2,m,n,j,5) + dc6*df(l1:l2,m,n,j)
!
        df(l1:l2,m,n,j) = c1*k(l1:l2,m,n,j,1) + &!c2*k(l1:l2,m,n,j,2) + &
                          c3*k(l1:l2,m,n,j,3) + c4*k(l1:l2,m,n,j,4) + &
                          !c5*k(l1:l2,m,n,j,5) +
                          c6*df(l1:l2,m,n,j)
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
            scal=  sqrt(f(l1:l2,m,n,j)**2  + k(l1:l2,m,n,j,1)**2 + farraymin(j)**2)
            errmaxs = max(maxval(abs(err/scal)),errmaxs)
          case ('cons_frac_err')
            ! Constrained fractional error: variable lower bound to avoid division by zero
            scal = max(abs(f(l1:l2,m,n,j)) + abs(df(l1:l2,m,n,j)),farraymin(j))
            errmaxs = max(maxval(abs(err/scal)),errmaxs)
          case ('rel_err')
            ! Relative error to f constrained with farraymin floor
            scal=max(abs(f(l1:l2,m,n,j)+df(l1:l2,m,n,j)),farraymin(j))
            errmaxs = max(maxval(abs(err)/scal),errmaxs)
          case ('abs_err')
            errmaxs = max(maxval(abs(err)),errmaxs)
          case ('cons_err')
            ! Relative error constraint for abs(f) > 1e-8
            scal=max(abs(f(l1:l2,m,n,j)),1e-8)
            errmaxs = max(maxval(abs(err/scal)),errmaxs)
          case ('none')
            ! No error check
        endselect
!
      enddo; enddo;

      enddo
!
!  Transferring auxiliaries from tmp to f.
!
      if (mfarray>mvar) f(:,:,:,mvar+1:mfarray) = tmp(:,:,:,mvar+1:mfarray)
!
! Divide your maximum error by the required accuracy
!
      errmaxs=errmaxs/eps_rkf
!
      call mpiallreduce_max(errmaxs,errmax,MPI_COMM_PENCIL)
!
    endsubroutine rkck
!***********************************************************************
    subroutine rkck3(f, df, p, errmax)
!
! Explicit third order Runge--Kutta--Fehlberg time stepping
!
      use Mpicomm, only: mpiallreduce_max,MPI_COMM_PENCIL
      use Equ, only: pde, impose_floors_ceilings
      use Shear, only: advance_shear
      use Boundcond, only: update_ghosts
!
! RK parameters by Cash and Karp
!
      real, parameter :: b21 = 0.2
      real, parameter :: b31 = 0.075
      real, parameter :: b32 = 0.225
      real, parameter :: b41 = 0.3
      real, parameter :: b42 = -0.9
      real, parameter :: b43 = 1.2
      real, parameter :: c1  = 0.35185185185185186 !19.0/54.0
      real, parameter :: c2  = 0.0
      real, parameter :: c3  = -0.37037037037037035 !-10.0/27.0
      real, parameter :: c4  = 1.0185185185185186 !55.0/54.0
      real, parameter :: dc1 = c1 + 1.5
      real, parameter :: dc2 = c2 - 2.5
      real, parameter :: dc3 = c3
      real, parameter :: dc4 = c4
      real, parameter :: substep1 = 0.2
      real, parameter :: substep2 = 0.3
      real, parameter :: substep3 = 0.6
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
      real, dimension(mx,my,mz,mvar,3) :: k
! Note: The tmp array will not use more memory than the temporary
!   array that would be implicitly created with calls like
!   call pde(f + b21*k(:,:,:,:,1), k(:,:,:,:,2), p)
      real, dimension (mx,my,mz,mfarray) :: tmp
      real, dimension(nx) :: scal, err
      real :: errmax, errmaxs, dtsub, told
      integer :: j
!
      df=0.
      errmax=0.
      errmaxs=0.
      k=0.
      told=t
!
!  FIRST SUBSTEP
!
      call pde(f, k(:,:,:,:,1), p)
      dtsub = dt*substep1
      if (lshear) then
        tmp = f
        call impose_floors_ceilings(tmp)
        call update_ghosts(tmp)  ! Necessary for non-FFT advection but unnecessarily overloading FFT advection
        call advance_shear(tmp, k(:,:,:,:,1), dtsub)
      endif
      do j=1,mvar; do n=n1,n2; do m=m1,m2
        k(l1:l2,m,n,j,1) = dt*k(l1:l2,m,n,j,1)
      enddo; enddo; enddo
!
      lfirst=.false.
!
!  Transferring auxiliaries from f to tmp.
!
      if (mfarray>mvar) tmp(:,:,:,mvar+1:mfarray) = f(:,:,:,mvar+1:mfarray)
      tmp(:,:,:,1:mvar) = f(:,:,:,1:mvar) + b21*k(:,:,:,:,1)
!
      call update_after_substep(tmp,k(:,:,:,:,1),dtsub,.false.)
      t = told + dtsub
!
!  SECOND SUBSTEP
!
      call pde(tmp, k(:,:,:,:,2), p)
      dtsub = dt*substep2
      do j=1,mvar; do n=n1,n2; do m=m1,m2
        k(l1:l2,m,n,j,2) = dt*k(l1:l2,m,n,j,2)
      enddo; enddo; enddo
!
      tmp(:,:,:,1:mvar) = f(:,:,:,1:mvar) + b31*k(:,:,:,:,1) &
                                          + b32*k(:,:,:,:,2)
      call update_after_substep(tmp,k(:,:,:,:,2),0.,.false.)
      t = told + dtsub
!
!
!  THIRD SUBSTEP
!
      call pde(tmp, k(:,:,:,:,3), p)
      dtsub = dt*substep3
      if (lshear) then
        call impose_floors_ceilings(tmp)
        call update_ghosts(tmp)
        call advance_shear(tmp, k(:,:,:,:,3), dtsub)
      endif
      do j=1,mvar; do n=n1,n2; do m=m1,m2
        k(l1:l2,m,n,j,3) = dt*k(l1:l2,m,n,j,3)
      enddo; enddo; enddo
!
      tmp(:,:,:,1:mvar) = f(:,:,:,1:mvar) + b41*k(:,:,:,:,1) &
                                          + b42*k(:,:,:,:,2) &
                                          + b43*k(:,:,:,:,3)
!
      call update_after_substep(tmp,k(:,:,:,:,3),dtsub,.false.)
      t = told + dtsub
!
!  FINALIZE STEP
!
      call pde(tmp, df, p)
      t = told + dt
!
      do j=1,mvar; do n=n1,n2; do m=m1,m2
        df(l1:l2,m,n,j) = dt*df(l1:l2,m,n,j)
!
        err = dc1*k(l1:l2,m,n,j,1) + dc2*k(l1:l2,m,n,j,2) + &
              dc3*k(l1:l2,m,n,j,3) + dc4*df(l1:l2,m,n,j)
!
        df(l1:l2,m,n,j) = c1*k(l1:l2,m,n,j,1) + c2*k(l1:l2,m,n,j,2) + &
                          c3*k(l1:l2,m,n,j,3) + c4*df(l1:l2,m,n,j)
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
            scal=  sqrt(f(l1:l2,m,n,j)**2  + k(l1:l2,m,n,j,1)**2 + farraymin(j)**2)
            errmaxs = max(maxval(abs(err/scal)),errmaxs)
          case ('cons_frac_err')
            ! Constrained fractional error: variable lower bound to avoid division by zero
            scal = max(abs(f(l1:l2,m,n,j)) + abs(df(l1:l2,m,n,j)),farraymin(j))
            errmaxs = max(maxval(abs(err/scal)),errmaxs)
          case ('rel_err')
            ! Relative error to f constrained with farraymin floor
            scal=max(abs(f(l1:l2,m,n,j)+df(l1:l2,m,n,j)),farraymin(j))
            errmaxs = max(maxval(abs(err)/scal),errmaxs)
          case ('abs_err')
            errmaxs = max(maxval(abs(err)),errmaxs)
          case ('cons_err')
            ! Relative error constraint for abs(f) > 1e-8
            scal=max(abs(f(l1:l2,m,n,j)),1e-8)
            errmaxs = max(maxval(abs(err/scal)),errmaxs)
          case ('none')
            ! No error check
        endselect
!
      enddo; enddo;

      enddo
!
!  Transferring auxiliaries from tmp to f.
!
      if (mfarray>mvar) f(:,:,:,mvar+1:mfarray) = tmp(:,:,:,mvar+1:mfarray)
!
! Divide your maximum error by the required accuracy
!
      errmaxs=errmaxs/eps_rkf
!
      call mpiallreduce_max(errmaxs,errmax,MPI_COMM_PENCIL)
!
    endsubroutine rkck3
!***********************************************************************
    subroutine update_after_substep(f,df,dtsub,llast)
!
!   Hooks for modifying f and df after the timestep is performed.
!
!  12-03-17/wlyra: coded
!  28-03-17/MR: removed update_ghosts; checks for already communicated variables enabled.
!
      use Hydro,    only: hydro_after_timestep
      use Energy,   only: energy_after_timestep
      use Magnetic, only: magnetic_after_timestep
      use Special,  only: special_after_timestep
      use Particles_main, only: particles_special_after_dtsub
!
      logical, intent(in) :: llast
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: dtsub
!
!  Enables checks to avoid unnecessary communication
!
      ighosts_updated=0
!
!  Dispatch to respective modules. The module which communicates
!  the biggest number of variables should come first here.
!
      if (lhydro)    call hydro_after_timestep   (f,df,dtsub)
      if (lmagnetic) call magnetic_after_timestep(f,df,dtsub)
      if (lenergy)   call energy_after_timestep  (f,df,dtsub)
      if (lspecial) then
        call special_after_timestep(f, df, dtsub, llast)
        if (lparticles) call particles_special_after_dtsub(f, dtsub)
      endif
!
!  Disables checks.
!
      ighosts_updated=-1
!
    endsubroutine update_after_substep
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Messages, only: fatal_error

    integer, parameter :: n_pars=0
    integer(KIND=ikind8), dimension(:) :: p_par

    call fatal_error('timestep_rkf','alpha_ts, beta_ts not defined')

    endsubroutine pushpars2c
!***********************************************************************
endmodule Timestep
