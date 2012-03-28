! $Id$
!
! Module for the super-time-stepping scheme (STS)
! for diffusive terms (Alexiades, V., Amiez, G., &
! Gremaud, P. 1996, Commun. Num. Meth. Eng.,  12, 31)
! which goes back to W. Gentzsch (ZAMM 58, T415-T416, 1978).
!
!
!  TODO: add i.e. particles, interstellar, shear
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
  contains
!***********************************************************************
    subroutine time_step(f,df,p)
!
!  Temporal advance for a diffusion problem using STS.
!  'itorder' plays the role of the parameter N in Alexiades paper, which
!  usually lies between 3 and 20, where lower values have higher accuracy
!  and higher values give larger speed-ups. 'nu_sts' is defined in run_pars
!  and must have vales between 0.0 and 1.0 (default 0.1).
!
!  17-march-11/gustavo: coded
!  08-march-12/bing+bourdin.kis: first stable
!
      use BorderProfiles, only: border_quenching
      use Equ, only: pde
      use Mpicomm, only: mpiallreduce_max
      use Special, only: special_after_timestep
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,mvar), intent(out) :: df
      type (pencil_case), intent(out) :: p
!
      real :: dt1, dt1_local, dt1_last=0.0
      real, dimension (itorder) :: tau_sts
      integer :: j
!
! Temporal loop over N substeps given by itorder
!
      do itsub=1,itorder
        lfirst=(itsub==1)
        llast=(itsub==itorder)
!
!  In contradiction to Runge-Kutta this scheme
!  integrates independent substeps
        df=0.
!
!  Change df according to the chosen physics modules.
!
        call pde(f,df,p)
!
        if (lfirst.and.ldt) then
          dt1_local=maxval(dt1_max(1:nx))
!  Timestep growth limiter
          if (ddt > 0.) dt1_local=max(dt1_local,dt1_last)
          call mpiallreduce_max(dt1_local,dt1)
          dt=1.0/dt1
          if (ddt/=0.) dt1_last=dt1_local/ddt
!
!  Compute the STS substeps
          call substeps(dt,tau_sts)
!
!  Set time step to the super-timestep
          dt = sum(tau_sts)
        endif
!
!  Time evolution of grid variables.
!
        do j=1,mvar; do n=n1,n2; do m=m1,m2
          if (lborder_profiles) call border_quenching(f,df,j,tau_sts(itsub))
          f(l1:l2,m,n,j)=f(l1:l2,m,n,j)+tau_sts(itsub)*df(l1:l2,m,n,j)
        enddo; enddo; enddo
!
        if (lspecial) call special_after_timestep(f,df,tau_sts(itsub))
        
!  Increase time.
!
        t = t + tau_sts(itsub)
!
      enddo
!
    endsubroutine time_step
!***********************************************************************
    subroutine substeps(dt_expl,tau)
!
!  Computes STS substeps based on the explicit time step.
!
!  17-march-11/gustavo: coded
!  28-Mar-2012/Bourdin.KIS: added reverse and permuted order
!
      use Messages, only: fatal_error
!
      real, dimension(itorder), intent(out) :: tau
      real, intent(in) :: dt_expl
      integer :: it, j
!
      do it = 1, itorder
        if (permute_sts == -1) then
          ! Reverse order:
          j = itorder + 1 - it
        elseif (permute_sts > 0) then
          if (modulo (itorder, permute_sts) == 0) call fatal_error( &
              'timestep_sts', "'permute_sts' must not be a devider or 'itorder'")
          ! Permuted order: (last substep stays last)
          j = (modulo (it * permute_sts - 1, itorder) + 1)
        elseif (permute_sts < 0) then
          if (modulo (itorder, -permute_sts) == 0) call fatal_error( &
              'timestep_sts', "'permute_sts' must not be a devider or 'itorder'")
          ! Permuted reverse order: (last substep becomes first)
          j = (modulo ((itorder + 1 - it) * (-permute_sts) - 1, itorder) + 1)
        else
          ! Regular order: (Default)
          j = it
        endif
        ! Alexiades (1996): (first: largest substep, last: smallest substep)
        tau(it) = dt_expl / ((-1+nu_sts)*cos(((2*j-1)*pi)/(2.0*itorder)) + 1+nu_sts)
        ! W. Gentzsch (1978): (first: smallest substep, last: largest substep)
        ! tau(it) = dt_expl / ((cos(((2*j-1)*pi)/(2.0*itorder)) + 1.0) + R_sts/(2.0*itorder**2))
      enddo
!
    endsubroutine substeps
!***********************************************************************
endmodule Timestep
  
