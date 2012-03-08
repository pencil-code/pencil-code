! $Id$
!
!  Timestepping routine corresponding to the use of LSODE to solve chemistry.
!  The transport equations are solved as usual using RK methods but the
!  chemistry ODEs are separated and solved implicitly using LSODE either
!  following a sequential (1 chemistry step) or symmetric (2 chemistry steps)
!  splitting scheme.
!
module Timestep
!
  use Cparam
  use Cdata
  use LsodeForChemistry
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
!  Runge Kutta advance, accurate to order itorder.
!  At the moment, itorder can be 1, 2, or 3.
!
!   2-apr-01/axel: coded
!  14-sep-01/axel: moved itorder to cdata
!
      use BorderProfiles, only: border_quenching
      use Equ
      use Interstellar, only: calc_snr_damp_int
      use Mpicomm
      use Particles_main
      use Shear, only: advance_shear
      use Special, only: special_after_timestep
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real :: ds
      real :: t0, t1
      real :: dt1, dt1_local, dt1_last=0.0
      integer :: j
!
      if (lroot .and. headt .and. llsode) print*, 'timestep_LSODE: '// &
          'Chemistry is solved using LSODE'
      if (lroot .and. headt .and. lsplit_second) print*, 'timestep_LSODE: '// &
          'Second order symmetric splitting procedure (Strang splitting)'
!
!  First step of the splitting procedure: chemistry
!  (only if second order symmetric splitting activated)
!
      if (llsode .and. lsplit_second .and. .not.headt) then
        lstep1=.false.
        t0=t
        t1=t+dt/2.
!
        call pde_chemistry(f,df,p)
!
        call lsode_fc(t0,t1,f,df)
!
      endif
!
!  Coefficients for up to order 3.
!
      if (itorder==1) then
        alpha_ts=(/ 0., 0., 0. /)
        beta_ts =(/ 1., 0., 0. /)
      elseif (itorder==2) then
        alpha_ts=(/ 0., -1./2., 0. /)
        beta_ts=(/ 1./2.,  1.,  0. /)
      elseif (itorder==3) then
        !alpha_ts=(/0., -2./3., -1./)
        !beta_ts=(/1./3., 1., 1./2./)
        !  use coefficients of Williamson (1980)
        alpha_ts=(/  0. ,  -5./9., -153./128. /)
        beta_ts=(/ 1./3., 15./16.,    8./15.  /)
      else
        if (lroot) print*,'Not implemented: itorder=',itorder
        call mpifinalize
      endif
!
!  dt_beta_ts may be needed in other modules (like Dustdensity) for fixed dt
!
      if (.not. ldt) dt_beta_ts=dt*beta_ts
      lstep1=.true.
!
!  Set up df and ds for each time sub.
!
      do itsub=1,itorder
        llast=(itsub==itorder)
        if (itsub==1) then
          lfirst=.true.
          df=0.
          ds=0.
        else
          lfirst=.false.
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
        ds=ds+1.
!
!  If we are in the first time substep we need to calculate timestep dt.
!  Only do it on the root processor, then broadcast dt to all others.
!
        if (lfirst.and.ldt) then
          dt1_local=maxval(dt1_max(1:nx))
          !Timestep growth limiter
          if (real(ddt) > 0.) dt1_local=max(dt1_local(1),dt1_last)
          call mpireduce_max(dt1_local,dt1,1)
          if (lroot) dt=1.0/dt1(1)
          !Timestep growth limiter
          if (ddt/=0.) dt1_last=dt1_local(1)/ddt
          call mpibcast_real(dt,1)
        endif
!
!  Calculate dt_beta_ts (e.g. for t=t+dt_beta_ts(itsub)*ds or for Dustdensity)
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
!  Second (or third) step of the splitting procedure: chemistry
!  (Done for each timestep as long as lsode is used)
!
      if (llsode) then
        lstep1=.false.
        t0=t
        if (lsplit_second) then
          t1=t+dt/2.
        else
          t1=t+dt
        endif
!
        call pde_chemistry(f,df,p)
!
        call lsode_fc(t0,t1,f,df)
!
      endif
!
    endsubroutine rk_2n
!***********************************************************************
endmodule Timestep
