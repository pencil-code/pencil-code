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
  real, parameter :: safety      =  0.9
!  maxerr typically ~1, so Cash-Karp index of -1/5 too slow to converge
  real, parameter :: dt_decrease = -2.5!-0.2
  real, parameter :: dt_increase = -0.1!-0.25
  real, parameter :: dt_increase_factor = 1.91!5.0
  real            :: errcon, farraymin, dt_next
!
  contains
!
!***********************************************************************
    subroutine initialize_timestep
!
      use Messages, only: fatal_error, warning

      if (lparticles) call fatal_error("initialize_timestep", "Particles are"// &
                                                 " not yet supported by the adaptive rkf scheme")
      if (itorder/=5) then
        call warning('initialize_timestep','itorder set to 5 for Runge-Kutta-Fehlberg')
        itorder=5
      endif
!
      if (dt==0.) then
        call warning('initialize_timestep','dt=0 not appropriate for Runge-Kutta-Fehlberg'// &
                     'set to 1e-6')
        dt=1e-6
      endif
!
      !General error condition: errcon ~1e-4 redundant with maxerr ~1. 0.1 more effective accelerator
      errcon = 0.1!(5.0/safety)**(1.0/dt_increase)
      !prefactor for absolute floor on value of farray in calculating err
      farraymin = 1e-8
!
      !ldt is set after read_persistent in rsnap, so dt0==0 used to read, but ldt=F for write
      ldt=.false.
      !overwrite the persistent time_step from dt0 in run.in if dt too high to initialize run
      if (dt0/=0.) dt=dt0
      dt_next=dt

    endsubroutine initialize_timestep
!***********************************************************************
    subroutine time_step(f,df,p)
!
!  Cash-Karp variant of Runge-Kutta-Fehlberg accurate to 5th order
!  To use this, set itorder to 5.
!
!  22-jun-06/tony: coded
!
      use Messages, only: warning

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real :: errmax, dt_temp
      real(KIND=rkind8) :: tnew, told!, time1, time2
      integer :: j,i
!
!  dt_beta_ts may be needed in other modules (like Dustdensity) for fixed dt
!
!      if (.not. ldt) dt_beta_ts=dt*beta_ts
!
      lfirst=.true.
      told=t
      dt=dt_next
      do i=1,20
        ! Do a Runge-Kutta step
        call rkck(f, df, p, errmax)
        ! Step succeeded so exit
        if (errmax <= safety) exit
        ! Step didn't succeed so decrease the time step
        dt_temp = safety*dt*((errmax/safety)**dt_decrease)
        !if (lroot) print*,"Decreasing dt, dt_temp, errmax",dt,dt_temp,errmax
        ! Don't decrease the time step by more than a factor of ten
        dt = sign(max(abs(dt_temp), 0.3/eps_rkf**0.05*abs(dt)), dt)
        tnew=told+dt
        if (tnew == told) then
          ! Guard against infinitesimal time steps
          call warning('time_step','Timestep underflow in Runge-Kutta-Fehlberg')
          ! Increase dt???
          exit
        endif
        t=told
      enddo
      if (lroot.and.ip==6) then
        print*,"time_step: rkck took",i,"iterations to converge to errmax",errmax
        print*,"time_step: tried time",dt_next,"reducing to",dt,"reduction rate",dt/dt_next
      endif
      call update_after_substep(f,df,dt,.true.)
!
! Time step to try next time
!
!
      if (errmax < errcon) then
        ! Increase the time step > dt_increase_factor
        dt_next = safety*dt*((errmax/2.)**dt_increase)
      else
        ! But not by more than a factor of 1.9 (dt_increase_factor)
        ! Previous try 5.0 inefficient: excessive decrease interations
        dt_next = dt_increase_factor*dt
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
      use Mpicomm, only: mpiallreduce_max,MPI_COMM_WORLD
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
      logical :: llogarithmic
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
        llogarithmic=((j==ilnrho.and..not.ldensity_nolog).or.(j==iss).or.(j==ilnTT))
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
          scal=  sqrt(f(l1:l2,m,n,j)**2  + k(l1:l2,m,n,j,1)**2 + 1e-30)
          errmaxs = max(maxval(abs(err/scal)),errmaxs)
        case ('cons_frac_err')
          ! Constant fractional error, prone to division by zero
          scal = abs(f(l1:l2,m,n,j)) + abs(df(l1:l2,m,n,j))
          errmaxs = max(maxval(abs(err/scal)),errmaxs)
        case ('rel_err')
          ! Relative error constraint with farraymin=1e-8 floor
          scal=max(abs(f(l1:l2,m,n,j)),farraymin)
          errmaxs = max(maxval(abs(err)/scal),errmaxs)
        case ('cons_err')
          ! Relative error constraint for f>1 or abs error
          scal = abs(f(l1:l2,m,n,j))
          scal = max(scal,1.)
          errmaxs = max(maxval(abs(err)/scal),errmaxs)
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
      call mpiallreduce_max(errmaxs,errmax,MPI_COMM_WORLD)
!
    endsubroutine rkck
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
