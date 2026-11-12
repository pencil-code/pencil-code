! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! CPARAM logical, parameter :: lcourant_dt = .true.
!***************************************************************
!
! Module for the super-time-stepping scheme (STS)
! for diffusive terms (Alexiades, V., Amiez, G., &
! Gremaud, P. 1996, Commun. Num. Meth. Eng.,  12, 31)
! which goes back to W. Gentzsch (ZAMM 58, T415-T416, 1978).
!
! This differs from timestep_sts.f90 module in that we make sure the superstep
! never exceeds the stability threshold of the non-diffusive terms (the idea
! being that STS works well only for parabolic equations).
!
! 05-Nov-2026/Kishore: made a copy of timestep_sts.f90
!
module Timestep
!
  use Cdata
!
  implicit none
!
  private
!
  include '../timestep.h'
!
  contains
!***********************************************************************
    subroutine initialize_timestep
!
      use Messages
!
      if (dt0 < 0.) dt = 0
      ldt = (dt==0.)
      num_substeps = itorder
!
!     These have special handling in timestep.f90
!
      if (lparticles) then
        call fatal_error('initialize_timestep', 'particles not supported')
      endif
!
      if (lpointmasses) then
        call fatal_error('initialize_timestep', 'pointmasses not supported')
      endif
!
      if (lode) then
        call fatal_error('initialize_timestep', 'ODE not supported')
      endif
!
      if (lsolid_cells) then
        call fatal_error('initialize_timestep', 'solid cells not supported')
      endif
!
      if (lborder_profiles) then
        call fatal_error('initialize_timestep', 'border profiles not supported')
      endif
!
      if (lshear) then
        call fatal_error('initialize_timestep', 'shear not supported')
      endif
!
!     TODO: split_update also not supported, but not sure how to check for that
!
    endsubroutine initialize_timestep
!***********************************************************************
    subroutine time_step(f,df,p)
!
!     Temporal advance for a diffusion problem using STS.
!     'itorder' plays the role of the parameter N in Alexiades paper, which
!     usually lies between 3 and 20, where lower values have higher accuracy
!     and higher values give larger speed-ups. 'nu_sts' is defined in run_pars
!     and must have vales between 0.0 and 1.0 (default 0.1).
!
!     17-march-11/gustavo: coded
!     08-march-12/bing+bourdin.kis: first stable
!     05-Nov-2026/Kishore: modified to limit outer timestep
!
      use BorderProfiles, only: border_quenching
      use Equ, only: pde, set_dt1_max
      use Mpicomm, only: mpiallreduce_max,MPI_COMM_PENCIL
      use Special, only: special_after_timestep
!
      real, contiguous,dimension(:,:,:,:), intent(inout) :: f
      real, contiguous,dimension(:,:,:,:), intent(out) :: df
      type (pencil_case), intent(out) :: p
!
      real :: dt1, dt1_local
      real, save :: dt1_last=0.0
      real, dimension (itorder) :: tau_sts
      real, dimension(nx) :: dt1_max_nodiffus
      real :: dt1_local_nodiffus, dt1_nodiffus
!
!     Temporal loop over N substeps given by itorder
!
      do itsub=1,itorder
        lfirst=(itsub==1)
        llast=(itsub==itorder)
!
!       Unlike Runge-Kutta, this scheme integrates independent substeps
        df=0.
!
!       Change df according to the chosen physics modules.
!
        call pde(f,df,p)
!
        if (lfirst.and.ldt) then
!
          call set_dt1_max(p, dt1_max_nodiffus, linclude_diffusion=.false.)
!
          dt1_local = maxval(dt1_max(1:nx))
          call mpiallreduce_max(dt1_local,dt1,MPI_COMM_PENCIL)
!
          dt1_local_nodiffus = maxval(dt1_max_nodiffus)
          call mpiallreduce_max(dt1_local_nodiffus,dt1_nodiffus,MPI_COMM_PENCIL)
!
          if (dt1_superstep(dt1) >= dt1_nodiffus) then
!           Diffusion is the dominant constraint, and the superstep does not exceed
!           the stability threshold for the non-diffusive terms.
!           do nothing
          else
!           Make sure the superstep never exceeds the stability threshold for
!           the non-diffusive terms
            dt1 = dt1_substep(dt1_nodiffus)
          endif
!
          !TODO: handle ddt as well?
!
!         Compute the STS substeps
          call substeps(1/dt1,tau_sts)
!
!         Set time step to the super-timestep
!         06-Nov-2026/Kishore: TODO: isn't this too late to set dt, since pde calls rhs_cpu which calls timestep_diagnostics?
          dt = sum(tau_sts)
        endif
!
!       Apply border quenching.
!
        !if (lborder_profiles) call border_quenching(f,df,dt_beta_ts(itsub)) !was in timestep_sts.f90
        if (lborder_profiles) call border_quenching(f,df,tau_sts(itsub))
!
!       Time evolution of grid variables.
!
        f(l1:l2,m1:m2,n1:n2,1:mvar) = f(l1:l2,m1:m2,n1:n2,1:mvar) + tau_sts(itsub)*df(l1:l2,m1:m2,n1:n2,1:mvar)
!
!       TODO: I suppose hydro_after_timestep and so on have to be called
!
        if (lspecial) call special_after_timestep(f,df,tau_sts(itsub),llast)

!       Increase time.
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
              'timestep_sts', "'permute_sts' must not be a divisor of 'itorder'")
          ! Permuted order: (last substep stays last)
          j = (modulo (it * permute_sts - 1, itorder) + 1)
        elseif (permute_sts < 0) then
          if (modulo (itorder, -permute_sts) == 0) call fatal_error( &
              'timestep_sts', "'permute_sts' must not be a divisor of 'itorder'")
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
    real function alexiades_fac()
!
!     The factor in equation 10 of Alexiades et al.
!
!     05-Nov-2026/Kishore: coded
!
      real :: tmp1, tmp2
      integer, save :: itorder_old=-1
      real, save :: fac, nu_sts_old=0
!
      if ((itorder/=itorder_old).or.(nu_sts/=nu_sts_old)) then
        tmp1 = (1 + sqrt(nu_sts))**(2*itorder)
        tmp2 = (1 - sqrt(nu_sts))**(2*itorder)
        fac = (2*sqrt(nu_sts)/itorder)*(tmp1 + tmp2)/(tmp1 - tmp2)
!
        itorder_old = itorder
        nu_sts_old = nu_sts
      endif
!
      alexiades_fac = fac
!
    endfunction alexiades_fac
!***********************************************************************
    real function dt1_superstep(dt1_expl)
!
!     Equation 10 of Alexiades et al.
!
!     05-Nov-2026/Kishore: coded
!
      real, intent(in) :: dt1_expl
!
      dt1_superstep = dt1_expl*alexiades_fac()
!
    endfunction dt1_superstep
!***********************************************************************
    real function dt1_substep(dt1_superstep)
!
!     Equation 10 of Alexiades et al.
!
!     05-Nov-2026/Kishore: coded
!
      real, intent(in) :: dt1_superstep
!
      dt1_substep = dt1_superstep/alexiades_fac()
!
    endfunction dt1_substep
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Messages, only: fatal_error

    integer, parameter :: n_pars=0
    integer(KIND=ikind8), dimension(:) :: p_par

    call fatal_error('timestep_sts','alpha_ts, beta_ts not defined')

    endsubroutine pushpars2c
!***********************************************************************
endmodule Timestep

