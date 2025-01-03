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
  integer         :: itter
  logical         :: fixed_dt=.false.
!
  contains
!***********************************************************************
    subroutine initialize_timestep
!
!  Coefficients for up to order 4(3).
!  14-oct-24/fred: adapted from timestep.f90 using low storage adaptive
!                  Runge-Kutta-Fehlberg derived scheme.
!                  Unlike timestep_rkf.f90 intermediate substep df are not
!                  retained but overwritten inside two farray size registers
!                  which alternate for f and df at each substep. 5th order and
!                  above require more than 5 substeps and would require
!                  replacing alpha and beta din cdata.f90 so are not yet
!                  implemented.
!                  Current tests using interstellar sample yield much lower
!                  timesteps and therefore longer integration times than RKFded,
!                  but this implementation can be more easily integrated into
!                  the GPU coupled code.
!
      use Messages, only: not_implemented, warning
      use General, only: itoa, rtoa
!
      if (itorder==1) then
        if (lroot) call warning('initialize_timestep','itorder 1 timestep is not adapted for itorder <3')
        alpha_ts=(/ 0.0, 0.0, 0.0, 0.0, 0.0 /)
        beta_ts =(/ 1.0, 0.0, 0.0, 0.0, 0.0 /)
        itter=1
      elseif (itorder==2) then
        if (lroot) call warning('initialize_timestep','itorder 2 timestep is not adapted for itorder <3')
        alpha_ts=(/   0.0, -1/2.0, 0.0, 0.0, 0.0 /)
        beta_ts =(/ 1/2.0,    1.0, 0.0, 0.0, 0.0 /)
        itter=2
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
        alpha_ts  =  (/11847461282814./36547543011857.,&
                        3943225443063./ 7078155732230.,&
                        -346793006927./ 4029903576067.,&
                                                0.0,0.0 /)
        itter=4
     elseif (itorder==4) then
        ! Explicit Runge-Kutta 4th vs 3rd order 2 Register 4-step scheme
        ! "C" indicate scheme compromises stability and accuracy criteria
        beta_ts=(/     1153189308089./22510343858157.,&
                       1772645290293./4653164025191., &
                      -1672844663538./4480602732383., &
                       2114624349019./3568978502595., &
                       5198255086312./14908931495163. /)
        beta_hat =  (/ 1016888040809./7410784769900., &
                      11231460423587./58533540763752.,&
                      -1563879915014./6823010717585., &
                        606302364029./971179775848.,  &
                       1097981568119./3980877426909.  /)
        alpha_ts=(/970286171893./4311952581923., &
                       6584761158862./12103376702013.,&
                       2251764453980./15575788980749.,&
                      26877169314380./34165994151039., 0.0 /)
        itter=5
      else
        call not_implemented('initialize_timestep','itorder= '//trim(itoa(itorder)))
      endif
!
!  indices applying to error to constrain timestep adjustments near to unity
!
      dt_increase=-1./(itorder+dtinc)
      dt_decrease=-1./(itorder-dtdec)
!
      !overwrite the persistent time_step from dt0 in run.in if dt
      !too high to initialize run
      if (dt0>0.) then
        dt=dt0
      elseif (dt0<0.) then
        fixed_dt=.true.
        dt=-dt0
      else
        if (dt==0) then
          call warning('initialize_timestep','dt=0 not appropriate for Runge-Kutta-Fehlberg'// &
                     'set to dt_epsi='//trim(rtoa(dt_epsi)))
          dt=dt_epsi
        endif
      endif
!
      if (eps_rkf0/=0.) eps_rkf=eps_rkf0
!
      ldt = .false.!(dt==0.)
!
    endsubroutine initialize_timestep
!***********************************************************************
    subroutine time_step(f,df,p)
!
!  14-oct-24/fred: coded
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
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mfarray,2) :: ftmp
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,ny,nz,mvar) :: errdf
      real, dimension (nx) :: scal
      type (pencil_case) :: p
      real :: dtsub, dt_temp
      integer :: j, iR1, iR2
!
!  Determine a lower bound for each variable j by which to normalise the error
!
      do j=1,mvar
        farraymin(j) = max(dt_ratio*maxval(abs(f(l1:l2,m1:m2,n1:n2,j))),dt_epsi)
      enddo
      if (lroot.and.it==1) print*,"farraymin",farraymin
      ftmp(:,:,:,:,1) = f
      ftmp(:,:,:,:,2) = f
      dt_beta_ts=dt*beta_ts
      dt_beta_hat=dt*(beta_ts-beta_hat)
      dt_alpha_ts=dt*alpha_ts
      errmax=0.
      errmaxs=0.
      errdf=0.
!
      substep_loop:do itsub=1,itter
        lfirst=(itsub==1)
        llast=(itsub==itter)
!
!  swap which indices of ftmp store current stage f or df in pde
!
        iR1 = 1+mod(itsub+1,2)
        iR2 = 1+mod(itsub,2)
!
        headtt = headt .and. lfirst .and. lroot
!
        if (lfirst) then
          ftmp(:,:,:,1:mvar,iR2)=0.
          !f=ftmp(:,:,:,:,iR1)
        else
          ftmp(:,:,:,1:mvar,iR2)=0.
          if (it_rmv>0) lrmv=.false.
        endif
!
!  Set up particle derivative array.
!
        if (lparticles) call particles_timestep_first(ftmp(:,:,:,:,iR1))
!
!  Set up point masses derivative array
!
        if (lpointmasses) call pointmasses_timestep_first(ftmp(:,:,:,:,iR1))
!
!  Set up ODE derivatives array
!
        if (lode) call ode_timestep_first
!
!  Set up solid_cells time advance
!
        if (lsolid_cells) call solid_cells_timestep_first(ftmp(:,:,:,:,iR1))
!
!  Change df according to the chosen physics modules.
!
        call pde(ftmp(:,:,:,:,iR1),ftmp(:,:,:,1:mvar,iR2),p)
!
        if (lode) call ode
!
        dtsub = dt_beta_ts(itsub)
!
!  Apply border quenching.
!
        if (lborder_profiles) call border_quenching(ftmp(:,:,:,:,iR1),df,dtsub)
!
!  Time evolution of point masses.
!
        if (lpointmasses) call pointmasses_timestep_second(ftmp(:,:,:,:,iR1))
!
!  Time evolution of particle variables.
!
        if (lparticles) call particles_timestep_second(ftmp(:,:,:,:,iR1))
!
! Time evolution of ODE variables.
!
        if (lode) call ode_timestep_second
!
!  Time evolution of solid_cells.
!  Currently this call has only keep_compiler_quiet so set ds=1
!
        if (lsolid_cells) call solid_cells_timestep_second(ftmp(:,:,:,:,iR1),dtsub,2.)
!
!  Advance deltay of the shear (and, optionally, perform shear advection
!  by shifting all variables and their derivatives).
!
        if (lshear) then
          call impose_floors_ceilings(ftmp(:,:,:,:,iR1))
          call update_ghosts(ftmp(:,:,:,:,iR1))  ! Necessary for non-FFT advection but unnecessarily overloading FFT advection
          call advance_shear(ftmp(:,:,:,:,iR1), ftmp(:,:,:,1:mvar,iR2), dtsub)
        endif
!
        call update_after_substep(ftmp(:,:,:,:,iR1),ftmp(:,:,:,1:mvar,iR2),dtsub,llast)
!
!  Time evolution of grid variables.
!
        if (.not. lgpu) then
          if (itorder>2) then
            errdf = errdf + dt_beta_hat(itsub)*ftmp(l1:l2,m1:m2,n1:n2,1:mvar,iR2)
          endif
          ftmp(l1:l2,m1:m2,n1:n2,1:mvar,iR1) =  ftmp(l1:l2,m1:m2,n1:n2,1:mvar,iR1) &
                                             +  dt_alpha_ts(itsub) * ftmp(l1:l2,m1:m2,n1:n2,1:mvar,iR2)
          ftmp(l1:l2,m1:m2,n1:n2,1:mvar,iR2) =  ftmp(l1:l2,m1:m2,n1:n2,1:mvar,iR1) &
                                             + (dtsub - dt_alpha_ts(itsub)) * ftmp(l1:l2,m1:m2,n1:n2,1:mvar,iR2)
          if (mfarray>mvar) ftmp(:,:,:,mvar+1:mfarray,iR2) = ftmp(:,:,:,mvar+1:mfarray,iR1)
        endif
!
!  Increase time.
!
        t = t + dtsub
!
      enddo substep_loop
!
!  Integrate operator split terms.
!
      call split_update(ftmp(:,:,:,:,iR2))
!
!  Check how much the maximum error satisfies the defined accuracy threshold
!
      if (itorder>2) then
        do j=1,mvar
          do m=m1,m2; do n=n1,n2
            select case (timestep_scaling(j))
              case ('cons_frac_err')
                !requires initial f state to be stored
                scal = max(abs(f(l1:l2,m,n,j))+abs(ftmp(l1:l2,m,n,j,iR2)-f(l1:l2,m,n,j)),farraymin(j))
                errmaxs = max(maxval(abs(errdf(:,m-nghost,n-nghost,j))/scal),errmaxs)
              case ('rel_err')
                !initial f state can be overwritten
                scal = max(abs(ftmp(l1:l2,m,n,j,iR2)),farraymin(j))
                !scal = max(abs(ftmp(l1:l2,m,n,j,iR2)),dt_epsi)
                errmaxs = max(maxval(abs(errdf(:,m-nghost,n-nghost,j))/scal),errmaxs)
              case ('abs_err')
                !initial f state can be overwritten
                errmaxs = max(maxval(abs(errdf(:,m-nghost,n-nghost,j))),errmaxs)
              case ('none')
                ! No error check
            endselect
          enddo; enddo
        enddo
!
! Divide your maximum error by the required accuracy
!
        errmaxs=errmaxs/eps_rkf
!
        call mpiallreduce_max(errmaxs,errmax,MPI_COMM_WORLD)
        if (.not. fixed_dt) then
          if (errmax > 1) then
            ! Step above error threshold so decrease the next time step
            dt_temp = safety*dt*errmax**dt_decrease
            if (lroot.and.ip==6787) print*,"time_step: it",it,"dt",dt,&
                 "to ",dt_temp,"at errmax",errmax
            ! Don't decrease the time step by more than a factor of ten
            dt = sign(max(abs(dt_temp), 0.1*abs(dt)), dt)
          else
            if (lroot.and.ip==6787) print*,"time_step increased: it",it,"dt",dt,&
                 "to ",dt*errmax**dt_increase,"at errmax",errmax
            dt = dt*errmax**dt_increase
          endif
        endif
      endif
      f = ftmp(:,:,:,:,iR2)
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
