! $Id$
!
!  Runge-Kutta-Fehlberg Low-Storage method
!  Christopher A. Kennedy, Mark H. Carpenter, R.Michael Lewis,
!  Low-storage, explicit Runge–Kutta schemes for the compressible Navier–Stokes equations,
!  Applied Numerical Mathematics,
!  Volume 35, Issue 3,
!  2000,
!  Pages 177-219,
!  ISSN 0168-9274,
!  https://doi.org/10.1016/S0168-9274(99)00141-5.
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
  real, dimension(mvar) :: farraymin
  real, dimension (5) :: beta_hat, dt_beta_hat, dt_alpha_ts
  real            :: dt_increase, dt_decrease, errmax, errmaxs
  real            :: safety=0.95
!
  contains
!***********************************************************************
    subroutine initialize_timestep
!
!  Coefficients for up to order 3.
!  26-oct-21/wlyra: added 2N-RK4 coefficients. The coefficients are for
!                   for 4th order, but use 5 stages. Since itorder is used
!                   in the code for the number of stages, it should be 5
!                   to use it.
!
      use Messages, only: not_implemented, warning
      use General, only: itoa, rtoa
!
      if (itorder==1) then
        if (lroot) call warning('initialize_timestep','itorder 1 timestep is not adapted for itorder <3')
        alpha_ts=(/ 0.0, 0.0, 0.0, 0.0, 0.0 /)
        beta_ts =(/ 1.0, 0.0, 0.0, 0.0, 0.0 /)
      elseif (itorder==2) then
        if (lroot) call warning('initialize_timestep','itorder 2 timestep is not adapted for itorder <3')
        alpha_ts=(/   0.0, -1/2.0, 0.0, 0.0, 0.0 /)
        beta_ts =(/ 1/2.0,    1.0, 0.0, 0.0, 0.0 /)
      elseif (itorder==3) then
        ! Explicit Runge-Kutta 3rd vs 2nd order 2 Register 4-step scheme
        ! "C" indicate scheme compromises stability and accuracy criteria
        beta_ts =  (/ 1017324711453./ 9774461848756.,&
                      8237718856693./13685301971492.,&
                     57731312506979./19404895981398.,&
                   -101169746363290./37734290219643.,&
                                0.0 /)
        beta_hat =  (/15763415370699./46270243929542.,&
                        514528521746./ 5659431552419.,&
                      27030193851939./ 9429696342944.,&
                     -69544964788955./30262026368149.,&
                                  0.0 /)
        alpha_ts  =  (/0.0,11847461282814./36547543011857.,&
                                    3943225443063./ 7078155732230.,&
                                    -346793006927./ 4029903576067.,&
                                                0.0 /)
     elseif (itorder==4) then
        ! Explicit Runge-Kutta 4th vs 3rd order 2 Register 4-step scheme
        ! "C" indicate scheme compromises stability and accuracy criteria
        beta_ts=(/     1153189308089./22510343858157.,&
                       1772645290293./4653164025191., &
                      -1672844663538./4480602732383., &
                       2114624349019./3568978502595., &
                       5198255086312./14908931495163. /)
        beta_hat =  (/ 1101688804080./7410784769900., &
                      11231460423587./58533540763752.,&
                      -1563879915014./6823010717585., &
                        606302364029./971179775848.,  &
                       1097981568119./3980877426909.  /)
        alpha_ts=(/0.0, 970286171893./4311952581923., &
                       6584761158862./12103376702013.,&
                       2251764453980./15575788980749.,&
                      26877169314380./34165994151039. /)
      else
        call not_implemented('initialize_timestep','itorder= '//trim(itoa(itorder)))
      endif
      dt_increase=-1./(itorder+0.5+dtinc)
      dt_decrease=-1./(itorder-0.5-dtdec)
!
      if (dt==0.) then
        if (lroot) call warning('initialize_timestep','dt=0 not appropriate for Runge-Kutta-Fehlberg'// &
                     'set to dt_epsi = '//trim(rtoa(dt_epsi)))
        dt=dt_epsi
      endif
!
      if (eps_rkf0/=0.) eps_rkf=eps_rkf0
!
      !overwrite the persistent time_step from dt0 in run.in if dt
      !too high to initialize run
      if (dt0/=0.) dt=dt0
!
      ldt = .false.!(dt==0.)
!
    endsubroutine initialize_timestep
!***********************************************************************
    subroutine time_step(f,df,p)
!
!   2-apr-01/axel: coded
!  14-sep-01/axel: moved itorder to cdata
!
      use Boundcond, only: update_ghosts
      use BorderProfiles, only: border_quenching
      use Equ, only: pde, impose_floors_ceilings
      use Mpicomm, only: mpiallreduce_max,MPI_COMM_WORLD
      use Particles_main, only: particles_timestep_first, &
          particles_timestep_second
      use PointMasses, only: pointmasses_timestep_first, &
          pointmasses_timestep_second
      use Solid_Cells, only: solid_cells_timestep_first, &
          solid_cells_timestep_second
      use Shear, only: advance_shear
!
      real, dimension (mx,my,mz,mfarray) :: f!, tmp
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,ny,nz,mvar) :: errdf
!      real, dimension (nx,ny,nz) :: scal
      type (pencil_case) :: p
      real :: ds, dtsub, dt_temp
      integer :: j
!
!  dt_beta_ts may be needed in other modules (like Dustdensity) for fixed dt.
!
! <<<<<<<<<<<<<<  the following should be omitted at some point <<<<<<<<<<<<<<
!  There is also an option to advance the time in progressively smaller
!  fractions of the current time. This is used to model the end of inflation,
!  when lfractional_tstep_negative should be used.
!  If dt=.5 and tstart=20, then one has -20, -10, -5, -2.5, etc.
!  From radiation era onward, lfractional_tstep_positive should be used
!  to make sure the dt used in run.in is positive.
!
      do j=1,mvar
        farraymin(j) = max(dt_ratio*maxval(abs(f(l1:l2,m1:m2,n1:n2,j))),dt_epsi)
      enddo
      if (lroot.and.it==1) print*,"farraymin",farraymin
      errmax=0.
      errmaxs=0.
      errdf=0.0
!
!  Set up df and ds for each time sub.
!
      !tmp = f
      substep_loop:do itsub=1,itorder+1

        lfirst=(itsub==1)
        llast=(itsub==itorder+1)

        headtt = headt .and. lfirst .and. lroot

        if (lfirst) then
          if (.not.lgpu) df=0.0
          !ds=0.0
        else
          !if (.not.lgpu) df=alpha_ts(itsub)*df !(could be subsumed into pde, but is dangerous!)
          !ds=alpha_ts(itsub)*ds
          if (it_rmv>0) lrmv=.false.
        endif
!
!  Set up particle derivative array.
!
        if (lparticles) call particles_timestep_first(f)
!
!  Set up point masses derivative array
!
        if (lpointmasses) call pointmasses_timestep_first(f)
!
!  Set up ODE derivatives array
!
        if (lode) call ode_timestep_first
!
!  Set up solid_cells time advance
!
        if (lsolid_cells) call solid_cells_timestep_first(f)
!
!  Change df according to the chosen physics modules.
!
        call pde(f,df,p)
        if (lode) call ode

        !ds=ds+1.0
!
!  If we are in the first time substep we need to calculate timestep dt.
!  Takes minimum over and distributes to all processors.
!  With GPUs this is done on the CUDA side.
!
!        if (lfirst.and.ldt.and..not.lgpu) call set_dt(maxval(dt1_max))
!
!  Calculate dt_beta_ts.
!
        dt_alpha_ts=dt*alpha_ts
        dt_beta_ts=dt*beta_ts
        dt_beta_hat=dt*(beta_ts-beta_hat)
        if (ip<=6) print*, 'time_step: iproc, dt=', iproc_world, dt  !(all have same dt?)
        dtsub = dt_beta_ts(itsub)
!
!  Apply border quenching.
!
        if (lborder_profiles) call border_quenching(f,df,dt_beta_ts(itsub))
!
!  Time evolution of grid variables.
!
        if (.not. lgpu) then
          !tmp(l1:l2,m1:m2,n1:n2,1:mvar) = f(l1:l2,m1:m2,n1:n2,1:mvar) + dt_alpha_ts(itsub)*df(l1:l2,m1:m2,n1:n2,1:mvar)
          f(l1:l2,m1:m2,n1:n2,1:mvar) =  f(l1:l2,m1:m2,n1:n2,1:mvar) &
                                                     + dt_beta_ts(itsub)*df(l1:l2,m1:m2,n1:n2,1:mvar)
          if (itorder>2) errdf = errdf + dt_beta_hat(itsub)*df(l1:l2,m1:m2,n1:n2,1:mvar)
        endif
!
!  Time evolution of point masses.
!
        if (lpointmasses) call pointmasses_timestep_second(f)
!
!  Time evolution of particle variables.
!
        if (lparticles) call particles_timestep_second(f)
!
! Time evolution of ODE variables.
!
        if (lode) call ode_timestep_second
!
!  Time evolution of solid_cells.
!  Currently this call has only keep_compiler_quiet so set ds=1
!
        if (lsolid_cells) call solid_cells_timestep_second(f,dt_beta_ts(itsub),1.)
!
!  Advance deltay of the shear (and, optionally, perform shear advection
!  by shifting all variables and their derivatives).
!
        if (lshear) then
          call impose_floors_ceilings(f)
          call update_ghosts(f)  ! Necessary for non-FFT advection but unnecessarily overloading FFT advection
          call advance_shear(f, df, dtsub)
        endif
!
        call update_after_substep(f,df,dtsub,llast)
!
        ! [PAB] according to MR this breaks the autotest.
        ! @Piyali: there must be a reason to add an additional global communication,
        ! could this be solved differently, and, if not, why? Can you enclose this in an if-clause, like above?
        !call update_ghosts(f)  ! Necessary for boundary driving in special module
!
!  Increase time.
!
        t = t + dtsub
!
      enddo substep_loop
!
!  Integrate operator split terms.
!
      call split_update(f)
!
      if (itorder>2) then
        do j=1,mvar
!
! Get the maximum error over the whole field
!
          select case (timestep_scaling(j))
            case ('cons_frac_err')
              ! Constrained fractional error: variable lower bound to avoid division by zero
              !scal = max(abs(f(l1:l2,m1:m2,n1:n2,j)),farraymin(j))
              errmaxs = max(maxval(abs(errdf(:,:,:,j))/&
                           max(abs(f(l1:l2,m1:m2,n1:n2,j))+abs(df(l1:l2,m1:m2,n1:n2,j)),farraymin(j))),errmaxs)
            case ('none')
              ! No error check
          endselect
!
        enddo
!
! Divide your maximum error by the required accuracy
!
        errmaxs=errmaxs/eps_rkf
!
        call mpiallreduce_max(errmaxs,errmax,MPI_COMM_WORLD)
        if (errmax > 1) then
          ! Step above error threshold so decrease the next time step
          if (lroot.and.ip==6787) print*,"time_step: it",it,"dt",dt,&
               "to ",safety*dt*errmax**dt_decrease,"at errmax",errmax
          dt_temp = safety*dt*errmax**dt_decrease
          ! Don't decrease the time step by more than a factor of ten
          dt = sign(max(abs(dt_temp), 0.1*abs(dt)), dt)
        else
          dt = dt*errmax**dt_increase
        endif
      endif
!
      if (ip<=6) print*,'TIMESTEP: iproc, dt=',iproc_world,dt
!
    endsubroutine time_step
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
      use Particles_main, only: particles_special_after_dtsub, particles_write_rmv
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
!  Flush the list of removed particles to the log file.
!
      if (lparticles) call particles_write_rmv
!
!  Disables checks.
!
      ighosts_updated=-1
!
    endsubroutine update_after_substep
!***********************************************************************
    subroutine ode_timestep_first

      if (lroot) then
        if (itsub==1) then
          df_ode = 0.0
        else
          df_ode=alpha_ts(itsub)*df_ode
        endif
      endif

    endsubroutine ode_timestep_first
!***********************************************************************
    subroutine ode

      use Special, only: dspecial_dt_ode

      if (lroot) then
        call dspecial_dt_ode
      endif

    endsubroutine ode
!***********************************************************************
    subroutine ode_timestep_second

      use Mpicomm, only: mpibcast

      if (lroot) f_ode = f_ode + dt_beta_ts(itsub)*df_ode
      call mpibcast(f_ode,n_odevars)

    endsubroutine ode_timestep_second
!***********************************************************************
    subroutine pushpars2c(p_par)

    !!use Syscalls, only: copy_addr

    integer, parameter :: n_pars=0
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    !!call copy_addr(alpha_ts,p_par(1))  ! (3)
    !!call copy_addr(beta_ts ,p_par(2))  ! (3)

    endsubroutine pushpars2c
!***********************************************************************
endmodule Timestep
