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
  include 'timestep.h'
!
  real :: dt_major = 0.0
!
  contains
!***********************************************************************
    subroutine initialize_timestep
! 
!  Coefficients for up to order 3.
!    
      use Messages, only: fatal_error
      use General, only: itoa
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
        call fatal_error('initialize_timestep','Not implemented: itorder= '// &
                         trim(itoa(itorder)))
      endif

    endsubroutine initialize_timestep
!***********************************************************************
    subroutine time_step(f, df, p)
!
!  Use Strang splitting to advance one time step.
!
!  22-jun-15/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type(pencil_case), intent(inout) :: p
!
      logical :: lfirstcall = .true.
      logical :: lout1
!
!  Save the time step if specified.
!
      init: if (lfirstcall) then
        dt_major = dt
        lfirstcall = .false.
      endif init
!
!  First half time-step with RK3.
!
      if (.not. ldt) dt = 0.5 * dt_major
      call rk3(f, df, p, .true.)
!
!  Full time-step with operator-split integration.
!
      dt = dt_major
      call split_update(f)
!
!  Turn off the calculation of the diagnostic output since it has been
!  done, if any.
!
      lout1 = lout
      lout = .false.
!
!  Second half time-step with RK3.
!
      dt = 0.5 * dt_major
      call rk3(f, df, p, .false.)
!
!  Retore the state of the diagnostic output.
!
      lout = lout1
!
    endsubroutine time_step
!***********************************************************************
    subroutine rk3(f, df, p, lfirst_half)
!
!   2-apr-01/axel: coded
!  14-sep-01/axel: moved itorder to cdata
!
      use Boundcond, only: update_ghosts
      use BorderProfiles, only: border_quenching
      use Equ, only: pde, impose_floors_ceilings
      use Mpicomm, only: mpiallreduce_max
      use Particles_main, only: particles_timestep_first, &
          particles_timestep_second
      use Shear, only: advance_shear
      use Special, only: special_after_timestep
      use Sub, only: shift_dt
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type(pencil_case), intent(inout) :: p
      logical, intent(in) :: lfirst_half
!
      real :: ds, dtsub
      real :: dt1, dt1_local, dt1_last=0.0
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
        if (lparticles) call particles_timestep_first(f)
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
        if (lfirst_half .and. lfirst .and. ldt) then
          dt1_local=maxval(dt1_max(1:nx))
          ! Timestep growth limiter
          if (real(ddt) > 0.) dt1_local=max(dt1_local,dt1_last)
          call mpiallreduce_max(dt1_local,dt1)
          dt_major = 1.0 / dt1
          if (loutput_varn_at_exact_tsnap) call shift_dt(dt_major)
          dt = 0.5 * dt_major
          ! Timestep growth limiter
          if (ddt/=0.) dt1_last=dt1_local/ddt
        endif
!
!  Calculate dt_beta_ts.
!
        if (ldt) dt_beta_ts=dt*beta_ts
        if (ip<=6) print*, 'time_step: iproc, dt=', iproc, dt  !(all have same dt?)
        dtsub = ds * dt_beta_ts(itsub)
!
!  Apply border quenching.
!
        if (lborder_profiles) call border_quenching(f,df,dt_beta_ts(itsub))
!
!  Time evolution of grid variables.
!
        f(l1:l2,m1:m2,n1:n2,1:mvar) = f(l1:l2,m1:m2,n1:n2,1:mvar) + dt_beta_ts(itsub)*df(l1:l2,m1:m2,n1:n2,1:mvar)
!
!  Time evolution of particle variables.
!
        if (lparticles) call particles_timestep_second(f)
!
!  Advance deltay of the shear (and, optionally, perform shear advection
!  by shifting all variables and their derivatives).
!
        advec: if (lshear) then
          call impose_floors_ceilings(f)
          call update_ghosts(f)  ! Necessary for non-FFT advection but unnecessarily overloading FFT advection
          call advance_shear(f, df, dtsub)
        endif advec
!
        if (lspecial) call special_after_timestep(f, df, dtsub, llast)
!
!  Increase time.
!
        t = t + dtsub
!
      enddo
!
    endsubroutine rk3
!***********************************************************************
    subroutine split_update(f)
!
!  Integrate operator split terms.
!
!  14-dec-14/ccyang: coded
!
      use Density, only: split_update_density
      use Energy, only: split_update_energy
      use Magnetic, only: split_update_magnetic
      use Viscosity, only: split_update_viscosity
      use Particles_main, only: split_update_particles
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
!  Dispatch to respective modules.
!
      if (ldensity) call split_update_density(f)
      if (lenergy) call split_update_energy(f)
      if (lmagnetic) call split_update_magnetic(f)
      if (lviscosity) call split_update_viscosity(f)
!
      if (lparticles) call split_update_particles(f, dt)
!
    endsubroutine split_update
!***********************************************************************
    subroutine pushpars2c(p_par)

    integer, parameter :: n_pars=2
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr_c(alpha_ts,p_par(1))  ! (3)
    call copy_addr_c(beta_ts ,p_par(2))  ! (3)

    endsubroutine pushpars2c
!***********************************************************************
endmodule Timestep
