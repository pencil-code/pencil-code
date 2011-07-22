! $Id$
! Experimental module to implement super_time_stepping 
! for diffusive terms (Alexiades, V., Amiez, G., & 
! Gremaud, P. 1996, Commun. Num. Meth. Eng.,  12, 31)
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
  public :: sts, substeps, rk_2n
!
  contains
!***********************************************************************
    subroutine sts(f,df,p)
!
!  Temporal advance for a diffusion problem ussing STS,
!  itorder plays the role of N in Alexiades paper and 
!  can be 5, 10, ... , nu_sts is defined in run_pars and
!  must have vales between 0 and 1 (default 0.1).
!
!  17-march-11/gustavo: coded
!
      use Equ, only: pde
      use Mpicomm, only: mpifinalize, mpireduce_max,mpibcast_real
      use Particles_main, only: particles_timestep_first, &
          particles_timestep_second
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real :: dt_sts,super_step
      real, dimension(1) :: dt1, dt1_local
      integer :: j
      real, dimension (itorder) :: tau_sts
      logical :: lfirstcall_sts=.true.
!
!  Before going into the loop, compute tau_sts
!
        if (ldt) then
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
! Getting an array with the substeps
        super_step=0
        call substeps(dt,tau_sts)
!
! Temporal loop over N substeps      
      do itsub=1,itorder         
        lfirst=(itsub==1)
        llast=(itsub==itorder)
!
!  Change df according to the chosen physics modules.
!
     call pde(f,df,p)  
!
!  If we are in the first time substep we need to calculate timestep dt.
!  Only do it on the root processor, then broadcast dt to all others.
!
        dt_sts = tau_sts(itorder-itsub+1)
!
!  Time evolution of grid variables.
!  (do this loop in pencils, for cache efficiency)
!
        do j=1,mvar; do n=n1,n2; do m=m1,m2
!ajwm Note to self... Just how much overhead is there in calling
!ajwm a sub this often...
           f(l1:l2,m,n,j) = f(l1:l2,m,n,j) + dt_sts*df(l1:l2,m,n,j)
        enddo; enddo; enddo
!
!  Increase time.
!
        super_step = super_step + dt_sts
        t = t + dt_sts
!
      enddo
      print*,' GG original step, super_step = ',dt,super_step
!
    endsubroutine sts
!***********************************************************************
    subroutine substeps(dtdiff,tau)
!
!  Computes STS substeps
!
!  17-march-11/gustavo: coded
!
      real, dimension (itorder) :: tau
      real :: dtdiff
      integer :: it
      do it=1,itorder
         tau(it) = dtdiff / ((-1.+nu_sts)*cos(((2.*it-1.)*pi)/(2.*itorder)) &
              + 1. + nu_sts)
      enddo
!      
    endsubroutine substeps
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
      use Equ, only: pde
      use Interstellar, only: calc_snr_damp_int
      use Mpicomm, only: mpifinalize, mpireduce_max,mpibcast_real,stop_it
      use Particles_main, only: particles_timestep_first, &
          particles_timestep_second
      use Shear, only: advance_shear
      use Special, only: special_after_timestep
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real :: ds
      real, dimension(1) :: dt1, dt1_local
      integer :: j
!
!  Coefficients for up to order 3.
!
      if (lsuper_time_stepping) then
!         call stop_it ('You are ussing the UNTESTED module STS')
         call sts(f,df,p)
         return
      endif
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
!
!  Set up df and ds for each time sub.
!
      do itsub=1,itorder
        lfirst=(itsub==1)
        llast=(itsub==itorder)
        if (lfirst) then
          df=0.
          ds=0.
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
          if (lborder_profiles) call border_quenching(df,j)
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
