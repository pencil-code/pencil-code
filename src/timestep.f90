! $Id$
!
! MODULE_DOC: Runge-Kutta time advance, accurate to order itorder.
! MODULE_DOC: At the moment, itorder can be 1, 2, or 3.
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
  contains
!***********************************************************************
    subroutine rk_2n(f,df,p)
!
!   2-apr-01/axel: coded
!  14-sep-01/axel: moved itorder to cdata
!
      use BorderProfiles, only: border_quenching
      use Equ, only: pde
      use Interstellar, only: calc_snr_damp_int
      use Mpicomm, only: mpifinalize, mpiallreduce_max
      use Particles_main, only: particles_timestep_first, &
          particles_timestep_second
      use Shear, only: advance_shear
      use Special, only: special_after_timestep
      use Snapshot, only: shift_dt
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real :: ds
      real :: dt1,dt1_local
      integer :: j
!
!  Coefficients for up to order 3.
!
      if (itorder==1) then
        alpha_ts=(/ 0.0, 0.0, 0.0 /)
        beta_ts =(/ 1.0, 0.0, 0.0 /)
      elseif (itorder==2) then
        alpha_ts=(/   0.0, -1/2.0, 0.0 /)
        beta_ts =(/ 1/2.0,    1.0, 0.0 /)
      elseif (itorder==3) then
        !alpha_ts=(/   0.0, -2/3.0, -1.0   /)
        !beta_ts =(/ 1/3.0,    1.0,  1/2.0 /)
        !  use coefficients of Williamson (1980)
        alpha_ts=(/   0.0, -5/9.0 , -153/128.0 /)
        beta_ts =(/ 1/3.0, 15/16.0,    8/15.0  /)
      else
        if (lroot) print*,'Not implemented: itorder=',itorder
        call mpifinalize
      endif
!
!  dt_beta_ts may be needed in other modules (like Dustdensity) for fixed dt.
!
      if (.not. ldt) dt_beta_ts=dt*beta_ts
!
!  Set up df and ds for each time sub.
!
      do itsub=1,itorder
        lfirst=(itsub==1)
        llast=(itsub==itorder)
        if (lfirst) then
          df=0.0
          ds=0.0
        else
          df=alpha_ts(itsub)*df !(could be subsumed into pde, but is dangerous!)
          ds=alpha_ts(itsub)*ds
        endif
!
!  Set up particle derivative array.
!
        if (lparticles) call particles_timestep_first()
!
!  Change df according to the chosen physics modules.
!
        call pde(f,df,p)
!
        ds=ds+1.0
!
!  If we are in the first time substep we need to calculate timestep dt.
!  Only do it on the root processor, then broadcast dt to all others.
!
        if (lfirst.and.ldt) then
          dt1_local=maxval(dt1_max(1:nx))
          ! Timestep growth limiter
          if (real(ddt) > 0.) dt1_local=max(dt1_local,dt1_last)
          call mpiallreduce_max(dt1_local,dt1)
          dt=1.0/dt1
          if (loutput_varn_at_exact_tsnap) call shift_dt(dt)
          ! Timestep growth limiter
          if (ddt/=0.) dt1_last=dt1_local/ddt
        endif
!
!  Calculate dt_beta_ts.
!
        if (ldt) dt_beta_ts=dt*beta_ts
        if (ip<=6) print*, 'rk_2n: iproc, dt=', iproc, dt  !(all have same dt?)
!
!  Add artificial damping at the location of SN explosions for a short time
!  after insertion.
!
        if (linterstellar) call calc_snr_damp_int(dt_beta_ts(itsub))
!
!  Time evolution of grid variables.
!  (do this loop in pencils, for cache efficiency)
!
        do j=1,mvar; do n=n1,n2; do m=m1,m2
!ajwm Note to self... Just how much overhead is there in calling
!ajwm a sub this often...
          if (lborder_profiles) call border_quenching(f,df,j,dt_beta_ts(itsub))
          f(l1:l2,m,n,j)=f(l1:l2,m,n,j)+dt_beta_ts(itsub)*df(l1:l2,m,n,j)
        enddo; enddo; enddo
!
!  Time evolution of particle variables.
!
        if (lparticles) call particles_timestep_second()
!
!  Advance deltay of the shear (and, optionally, perform shear advection
!  by shifting all variables and their derivatives).
!
        if (lshear) call advance_shear(f,df,dt_beta_ts(itsub)*ds)
!
        if (lspecial) &
            call special_after_timestep(f,df,dt_beta_ts(itsub)*ds)
!
!  Increase time.
!
        t = t + dt_beta_ts(itsub)*ds
!
      enddo
!
    endsubroutine rk_2n
!***********************************************************************
endmodule Timestep
