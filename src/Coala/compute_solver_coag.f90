!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


!-------------------------------------------------
! MODULE: Functions to compute coagulation solver 
!         for 1 hydro time-step
!-------------------------------------------------

module compute_solver_coag
   use precision
   use solver_DG
   use progress_bar

   implicit none
contains


!-----------------------------------------------------
! Original version to compute functions for DG scheme
!-----------------------------------------------------

!> @brief Compute coagulation solver (simple collision kernels) for 1 hydro time-step 
!!
!! DG scheme k=0, piecewise constant approximation
!!
!! @param[in]    eps           minimum value for mass distribution approximation gij
!! @param[in]    coeff_CFL     timestep coefficient for stability of the SSPRK order 3 scheme
!! @param[in]    nbins         number of dust bins
!! @param[in]    massgrid      grid of masses, borders value of mass bins
!! @param[in]    tabflux_coag  3D array to evaluate coagulation flux
!! @param[in]    dthydro       hydro timestep
!! @param[in]    gij           components of g on the polynomial basis
!! @param[out]   gijnew        evolved components of g on the polynomial basis after dthydro
!! @param[out]   nsub          number of subcycling coagulation timestep to reach dthydro
!! @param[out]   ndt           number of hydro timestep, when coagulation CFL > dthydro
subroutine compute_coag_k0(eps,coeff_CFL,nbins,massgrid,tabflux_coag,dthydro,gij,gijnew,nsub,ndt)
   implicit none
   integer,  intent(in)   :: nbins
   real(wp), intent(in)   :: massgrid(nbins+1)
   real(wp), intent(in)   :: tabflux_coag(nbins,nbins,nbins)
   real(wp), intent(in)   :: dthydro,coeff_CFL,eps
   real(wp), intent(in)   :: gij(nbins)
   real(wp), intent(out)  :: gijnew(nbins)
   integer,  intent(out)  :: nsub,ndt

   real(wp) :: dtCFLsub,dtsub,dtlast,dt,gijsub_in(nbins),gijsub_out(nbins)

   nsub = 0
   ndt = 0

   !evaluate coagulation CFL
   dtCFLsub = compute_CFL_k0(eps,nbins,massgrid,gij,tabflux_coag)
   dtCFLsub = coeff_CFL*dtCFLsub

   !compare hydro timestep and coagulation CFL
   dt = min(dtCFLsub,dthydro)

   !coagulation subcycling timesteps
   if (dt<dthydro) then
      gijsub_in = gij
      dtsub = 0._wp
      do while ( dtsub<dthydro .and. dthydro-dtsub>dtCFLsub)
         dtsub=dtsub+dtCFLsub

         nsub = nsub + 1

         call time_solver_k0(eps,nbins,massgrid,gijsub_in,tabflux_coag,dtCFLsub,gijsub_out)
         gijsub_in = gijsub_out


         dtCFLsub = compute_CFL_k0(eps,nbins,massgrid,gijsub_in,tabflux_coag)
         dtCFLsub = coeff_CFL*dtCFLsub

      enddo

      !last timestep to reach dthydro
      dtlast = dthydro-dtsub

      nsub = nsub + 1

      call time_solver_k0(eps,nbins,massgrid,gijsub_in,tabflux_coag,dtlast,gijnew)

   !when coagulation CFL > hydro timstep
   else

      ndt = ndt + 1

      call time_solver_k0(eps,nbins,massgrid,gij,tabflux_coag,dt,gijnew)

   endif

end subroutine compute_coag_k0




!> @brief Compute coagulation solver (physical collision kernel) for 1 hydro time-step 
!!
!! DG scheme k=0, piecewise constant approximation
!!
!! @param[in]    eps           minimum value for mass distribution approximation gij
!! @param[in]    coeff_CFL     timestep coefficient for stability of the SSPRK order 3 scheme
!! @param[in]    nbins         number of dust bins
!! @param[in]    massgrid      grid of masses, borders value of mass bins
!! @param[in]    tabflux_coag  3D array to evaluate coagulation flux
!! @param[in]    dv             2D array of the differential velocity between grains in bins lp and l
!! @param[in]    dthydro       hydro timestep
!! @param[in]    gij           components of g on the polynomial basis
!! @param[out]   gijnew        evolved components of g on the polynomial basis after dthydro
!! @param[out]   nsub          number of subcycling coagulation timestep to reach dthydro
!! @param[out]   ndt           number of hydro timestep, when coagulation CFL > dthydro
subroutine compute_coag_k0_kdv(eps,coeff_CFL,nbins,massgrid,tabflux_coag,dv,dthydro,gij,gijnew,nsub,ndt)
   implicit none
   integer,  intent(in)   :: nbins
   real(wp), intent(in)   :: massgrid(nbins+1)
   real(wp), intent(in)   :: tabflux_coag(nbins,nbins,nbins),dv(nbins,nbins)
   real(wp), intent(in)   :: dthydro,coeff_CFL,eps
   real(wp), intent(in)   :: gij(nbins)
   real(wp), intent(out)  :: gijnew(nbins)
   integer,  intent(out)  :: nsub,ndt

   real(wp) :: dtCFLsub,dtsub,dtlast,dt,gijsub_in(nbins),gijsub_out(nbins)

   nsub = 0
   ndt = 0

   !evaluate coagulation CFL
   dtCFLsub = compute_CFL_k0_kdv(eps,nbins,massgrid,gij,tabflux_coag,dv)
   dtCFLsub = coeff_CFL*dtCFLsub

   !compare hydro timestep and coagulation CFL
   dt = min(dtCFLsub,dthydro)

   !coagulation subcycling timesteps
   if (dt<dthydro) then
      gijsub_in = gij
      dtsub = 0._wp
      do while ( dtsub<dthydro .and. dthydro-dtsub>dtCFLsub)
         dtsub=dtsub+dtCFLsub

         nsub = nsub + 1

         call time_solver_k0_kdv(eps,nbins,massgrid,gijsub_in,tabflux_coag,dv,dtCFLsub,gijsub_out)
         gijsub_in = gijsub_out

         dtCFLsub = compute_CFL_k0_kdv(eps,nbins,massgrid,gijsub_in,tabflux_coag,dv)
         dtCFLsub = coeff_CFL*dtCFLsub

      enddo

      !last timestep to reach dthydro
      dtlast = dthydro-dtsub

      nsub = nsub + 1

      call time_solver_k0_kdv(eps,nbins,massgrid,gijsub_in,tabflux_coag,dv,dtlast,gijnew)

   !when coagulation CFL > hydro timestep
   else

      ndt = ndt + 1

      call time_solver_k0_kdv(eps,nbins,massgrid,gij,tabflux_coag,dv,dt,gijnew)

   endif

end subroutine compute_coag_k0_kdv



!> @brief Compute coagulation solver (simple collision kernels) for 1 hydro time-step 
!!
!! DG scheme k>0, piecewise polynomial approximation
!!
!! @param[in]    eps               minimum value for mass distribution approximation gij
!! @param[in]    coeff_CFL         timestep coefficient for stability of the SSPRK order 3 scheme
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    massbins          arithmetic mean value of massgrid for each mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[in]    tabflux_coag      5D array to evaluate coagulation flux
!! @param[in]    tabintflux_coag   6D array to evaluate the term including the integral of coagulation flux
!! @param[in]    dthydro           hydro timestep
!! @param[in]    gij               components of g on the polynomial basis
!! @param[out]   gijnew            evolved components of g on the polynomial basis after dthydro
!! @param[out]   nsub              number of subcycling coagulation timestep to reach dthydro
!! @param[out]   ndt               number of hydro timestep, when coagulation CFL > dthydro
subroutine compute_coag(eps,coeff_CFL,nbins,kpol,massgrid,massbins,mat_coeffs_leg,&
                              tabflux_coag,tabintflux_coag,&
                              dthydro,gij,gijnew,nsub,ndt)
   implicit none
   integer,  intent(in)  :: nbins,kpol
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins),mat_coeffs_leg(kpol+1,kpol+1)
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   real(wp), intent(in)  :: tabintflux_coag(nbins,kpol+1,nbins,nbins,kpol+1,kpol+1)
   real(wp), intent(in)  :: dthydro,coeff_CFL,eps
   real(wp), intent(in)  :: gij(nbins,kpol+1)
   integer,  intent(out) :: nsub,ndt
   real(wp), intent(out) :: gijnew(nbins,kpol+1)


   real(wp) :: dtCFLsub,dtsub,dtlast,dt,gijsub_in(nbins,kpol+1),gijsub_out(nbins,kpol+1)

   nsub = 0
   ndt = 0

   !evaluate coagulation CFL
   dtCFLsub = compute_CFL(eps,nbins,kpol,massgrid,gij,tabflux_coag)
   dtCFLsub = coeff_CFL*dtCFLsub

   !compare hydro timestep and coagulation CFL
   dt = min(dtCFLsub,dthydro)

   !coagulation subcycling timesteps
   if (dt<dthydro) then
      gijsub_in = gij
      dtsub = 0._wp
      do while ( dtsub<dthydro .and. dthydro-dtsub>dtCFLsub)
         dtsub=dtsub+dtCFLsub

         nsub = nsub + 1

         call time_solver(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,&
                          gijsub_in,tabflux_coag,tabintflux_coag,dtCFLsub,gijsub_out)
         gijsub_in = gijsub_out

         dtCFLsub = compute_CFL(eps,nbins,kpol,massgrid,gijsub_in,tabflux_coag)
         dtCFLsub = coeff_CFL*dtCFLsub

      enddo

      !last timestep to reach dthydro
      dtlast = dthydro-dtsub

      nsub = nsub + 1

      call time_solver(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gijsub_in,tabflux_coag,tabintflux_coag,dtlast,gijnew)

   !when coagulation CFL > hydro timestep
   else

      ndt = ndt + 1

      call time_solver(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tabflux_coag,tabintflux_coag,dt,gijnew)

   endif


end subroutine compute_coag





! 


!> @brief Compute coagulation solver (physical kernel) for 1 hydro time-step 
!!
!! DG scheme k>0, piecewise polynomial approximation
!!
!! @param[in]    eps               minimum value for mass distribution approximation gij
!! @param[in]    coeff_CFL         timestep coefficient for stability of the SSPRK order 3 scheme
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    massbins          arithmetic mean value of massgrid for each mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[in]    tabflux_coag      5D array to evaluate coagulation flux
!! @param[in]    tabintflux_coag   6D array to evaluate the term including the integral of coagulation flux
!! @param[in]    dv                2D array of the differential velocity between grains in bins lp and l
!! @param[in]    dthydro           hydro timestep
!! @param[in]    gij               components of g on the polynomial basis
!! @param[out]   gijnew            evolved components of g on the polynomial basis after dthydro
!! @param[out]   nsub              number of subcycling coagulation timestep to reach dthydro
!! @param[out]   ndt               number of hydro timestep, when coagulation CFL > dthydro
subroutine compute_coag_kdv(eps,coeff_CFL,nbins,kpol,massgrid,massbins,mat_coeffs_leg,&
                              tabflux_coag,tabintflux_coag,dv,&
                              dthydro,gij,gijnew,nsub,ndt)
   implicit none
   integer,  intent(in)  :: nbins,kpol
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins),mat_coeffs_leg(kpol+1,kpol+1)
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   real(wp), intent(in)  :: tabintflux_coag(nbins,kpol+1,nbins,nbins,kpol+1,kpol+1)
   real(wp), intent(in)  :: dv(nbins,nbins)
   real(wp), intent(in)  :: dthydro,coeff_CFL,eps
   real(wp), intent(in)  :: gij(nbins,kpol+1)
   integer,  intent(out) :: nsub,ndt
   real(wp), intent(out) :: gijnew(nbins,kpol+1)

   real(wp) :: dtCFLsub,dtsub,dtlast,dt,gijsub_in(nbins,kpol+1),gijsub_out(nbins,kpol+1)

   nsub = 0
   ndt = 0

   !evaluate coagulation CFL
   dtCFLsub = compute_CFL_kdv(eps,nbins,kpol,massgrid,gij,tabflux_coag,dv)
   dtCFLsub = coeff_CFL*dtCFLsub

   !compare hydro timestep and coagulation CFL
   dt = min(dtCFLsub,dthydro)

   !coagulation subcycling timesteps
   if (dt<dthydro) then
      gijsub_in = gij
      dtsub = 0._wp
      do while ( dtsub<dthydro .and. dthydro-dtsub>dtCFLsub)
         dtsub=dtsub+dtCFLsub

         nsub = nsub + 1

         call time_solver_kdv(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gijsub_in,&
                           tabflux_coag,tabintflux_coag,dv,dtCFLsub,gijsub_out)
         gijsub_in = gijsub_out

         dtCFLsub = compute_CFL_kdv(eps,nbins,kpol,massgrid,gijsub_in,tabflux_coag,dv)
         dtCFLsub = coeff_CFL*dtCFLsub

      enddo

      !last timestep to reach dthydro
      dtlast = dthydro-dtsub

      nsub = nsub + 1

      call time_solver_kdv(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gijsub_in,tabflux_coag,tabintflux_coag,dv,dtlast,gijnew)

   !when coagulation CFL > hydro timestep
   else

      ndt = ndt + 1

      call time_solver_kdv(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tabflux_coag,tabintflux_coag,dv,dt,gijnew)

   endif


end subroutine compute_coag_kdv








!-----------------------------------------------------
! Optimised version to compute functions for DG scheme
!-----------------------------------------------------

!> @brief Compute coagulation solver (simple kernels) for 1 hydro time-step
!!
!! with indices optimisation (only with gij > eps)
!!
!! DG scheme k=0, piecewise constant approximation
!!
!! @param[in]    eps           minimum value for mass distribution approximation gij
!! @param[in]    coeff_CFL     timestep coefficient for stability of the SSPRK order 3 scheme
!! @param[in]    nbins         number of dust bins
!! @param[in]    massgrid      grid of masses, borders value of mass bins
!! @param[in]    massbins      arithmetic mean value of massgrid for each mass bins
!! @param[in]    tabflux_coag  3D array to evaluate coagulation flux
!! @param[in]    dthydro       hydro timestep
!! @param[in]    gij           components of g on the polynomial basis
!! @param[out]   gijnew        evolved components of g on the polynomial basis after dthydro
!! @param[out]   nsub          number of subcycling coagulation timestep to reach dthydro
!! @param[out]   ndt           number of hydro timestep, when coagulation CFL > dthydro
subroutine compute_coag_k0_opt(eps,coeff_CFL,nbins,massgrid,massbins,tabflux_coag,dthydro,gij,gijnew,nsub,ndt)
   implicit none
   integer,  intent(in)   :: nbins
   real(wp), intent(in)   :: massgrid(nbins+1),massbins(nbins)
   real(wp), intent(in)   :: tabflux_coag(nbins,nbins,nbins)
   real(wp), intent(in)   :: dthydro,coeff_CFL,eps
   real(wp), intent(in)   :: gij(nbins)
   real(wp), intent(out)  :: gijnew(nbins)
   integer,  intent(out)  :: nsub,ndt

   real(wp) :: dtCFLsub,dtsub,dtlast,dt,gijsub_in(nbins),gijsub_out(nbins)

   nsub = 0
   ndt = 0

   !evaluate coagulation CFL
   dtCFLsub = compute_CFL_k0_opt(eps,nbins,massgrid,massbins,gij,tabflux_coag)
   dtCFLsub = coeff_CFL*dtCFLsub
   
   !compare hydro timestep and coagulation CFL
   dt = min(dtCFLsub,dthydro)

   !coagulation subcycling timesteps
   if (dt<dthydro) then
      gijsub_in = gij
      dtsub = 0._wp
      do while ( dtsub<dthydro .and. dthydro-dtsub>dtCFLsub)
         dtsub=dtsub+dtCFLsub

         nsub = nsub + 1

         call time_solver_k0_opt(eps,nbins,massgrid,massbins,gijsub_in,tabflux_coag,dtCFLsub,gijsub_out)
         gijsub_in = gijsub_out

         dtCFLsub = compute_CFL_k0_opt(eps,nbins,massgrid,massbins,gijsub_in,tabflux_coag)
         dtCFLsub = coeff_CFL*dtCFLsub

      enddo

      !last timestep to reach dthydro
      dtlast = dthydro-dtsub

      nsub = nsub + 1

      call time_solver_k0_opt(eps,nbins,massgrid,massbins,gijsub_in,tabflux_coag,dtlast,gijnew)

   !when coagulation CFL > hydro timstep
   else

      ndt = ndt + 1

      call time_solver_k0_opt(eps,nbins,massgrid,massbins,gij,tabflux_coag,dt,gijnew)

   endif

end subroutine compute_coag_k0_opt


!> @brief Compute coagulation solver (physical collision kernel) for 1 hydro time-step
!! 
!! with indices optimisation (only with gij > eps)
!!
!! DG scheme k=0, piecewise constant approximation
!!
!! @param[in]    eps           minimum value for mass distribution approximation gij
!! @param[in]    coeff_CFL     timestep coefficient for stability of the SSPRK order 3 scheme
!! @param[in]    nbins         number of dust bins
!! @param[in]    massgrid      grid of masses, borders value of mass bins
!! @param[in]    gij           components of g on the polynomial basis
!! @param[in]    tabflux_coag  3D array to evaluate coagulation flux
!! @param[in]    dv             2D array of the differential velocity between grains in bins lp and l
!! @param[in]    dthydro       hydro timestep
!! @param[out]   gijnew        evolved components of g on the polynomial basis after dthydro
!! @param[out]   nsub          number of subcycling coagulation timestep to reach dthydro
!! @param[out]   ndt           number of hydro timestep, when coagulation CFL > dthydro
subroutine compute_coag_k0_kdv_opt(eps,coeff_CFL,nbins,massgrid,massbins,&
                                       tabflux_coag,dv,dthydro,gij,gijnew,nsub,ndt)
   implicit none
   integer,  intent(in)   :: nbins
   real(wp), intent(in)   :: massgrid(nbins+1),massbins(nbins)
   real(wp), intent(in)   :: tabflux_coag(nbins,nbins,nbins),dv(nbins,nbins)
   real(wp), intent(in)   :: dthydro,coeff_CFL,eps
   real(wp), intent(in)   :: gij(nbins)
   real(wp), intent(out)  :: gijnew(nbins)
   integer,  intent(out)  :: nsub,ndt

   real(wp) :: dtCFLsub,dtsub,dtlast,dt,gijsub_in(nbins),gijsub_out(nbins)

   nsub = 0
   ndt = 0

   !evaluate coagulation CFL
   dtCFLsub = compute_CFL_k0_kdv_opt(eps,nbins,massgrid,massbins,gij,tabflux_coag,dv)
   dtCFLsub = coeff_CFL*dtCFLsub
   ! print*,"dtCFLsub=",dtCFLsub

   !compare hydro timestep and coagulation CFL
   dt = min(dtCFLsub,dthydro)

   !coagulation subcycling timesteps
   if (dt<dthydro) then
      gijsub_in = gij
      dtsub = 0._wp
      do while ( dtsub<dthydro .and. dthydro-dtsub>dtCFLsub)
         dtsub=dtsub+dtCFLsub

         nsub = nsub + 1
         call time_solver_k0_kdv_opt(eps,nbins,massgrid,massbins,gijsub_in,tabflux_coag,dv,dtCFLsub,gijsub_out)
         gijsub_in = gijsub_out

         dtCFLsub = compute_CFL_k0_kdv_opt(eps,nbins,massgrid,massbins,gijsub_in,tabflux_coag,dv)
         dtCFLsub = coeff_CFL*dtCFLsub
         ! print*,"dtCFLsub=",dtCFLsub

      enddo

      !last timestep to reach dthydro
      dtlast = dthydro-dtsub

      nsub = nsub + 1

      call time_solver_k0_kdv_opt(eps,nbins,massgrid,massbins,gijsub_in,tabflux_coag,dv,dtlast,gijnew)

   !when coagulation CFL > hydro timestep
   else

      ndt = ndt + 1

      call time_solver_k0_kdv_opt(eps,nbins,massgrid,massbins,gij,tabflux_coag,dv,dt,gijnew)

   endif

end subroutine compute_coag_k0_kdv_opt



!> @brief Compute coagulation solver (simple collision kernels) for 1 hydro time-step
!! 
!! with indices optimisation (only with gij > eps)
!!
!! DG scheme k>0, piecewise polynomial approximation
!!
!! @param[in]    eps               minimum value for mass distribution approximation gij
!! @param[in]    coeff_CFL         timestep coefficient for stability of the SSPRK order 3 scheme
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    massbins          arithmetic mean value of massgrid for each mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[in]    tabflux_coag      5D array to evaluate coagulation flux
!! @param[in]    tabintflux_coag   6D array to evaluate the term including the integral of coagulation flux
!! @param[in]    dthydro           hydro timestep
!! @param[in]    gij               components of g on the polynomial basis
!! @param[out]   gijnew            evolved components of g on the polynomial basis after dthydro
!! @param[out]   nsub              number of subcycling coagulation timestep to reach dthydro
!! @param[out]   ndt               number of hydro timestep, when coagulation CFL > dthydro
subroutine compute_coag_opt(eps,coeff_CFL,nbins,kpol,massgrid,massbins,mat_coeffs_leg,&
                              tabflux_coag,tabintflux_coag,&
                              dthydro,gij,gijnew,nsub,ndt)
   implicit none
   integer,  intent(in)  :: nbins,kpol
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins),mat_coeffs_leg(kpol+1,kpol+1)
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   real(wp), intent(in)  :: tabintflux_coag(nbins,kpol+1,nbins,nbins,kpol+1,kpol+1)
   real(wp), intent(in)  :: dthydro,coeff_CFL,eps
   real(wp), intent(in)  :: gij(nbins,kpol+1)
   integer,  intent(out) :: nsub,ndt
   real(wp), intent(out) :: gijnew(nbins,kpol+1)

   real(wp) :: dtCFLsub,dtsub,dtlast,dt,gijsub_in(nbins,kpol+1),gijsub_out(nbins,kpol+1)

   nsub = 0
   ndt = 0

   !evaluate coagulation CFL
   dtCFLsub = compute_CFL_opt(eps,nbins,kpol,massgrid,massbins,gij,tabflux_coag)
   dtCFLsub = coeff_CFL*dtCFLsub


   !compare hydro timestep and coagulation CFL
   dt = min(dtCFLsub,dthydro)

   !coagulation subcycling timesteps
   if (dt<dthydro) then
      gijsub_in = gij
      dtsub = 0._wp
      do while ( dtsub<dthydro .and. dthydro-dtsub>dtCFLsub)
         dtsub=dtsub+dtCFLsub

         nsub = nsub + 1

         call time_solver_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gijsub_in,&
                              tabflux_coag,tabintflux_coag,dtCFLsub,gijsub_out)
         gijsub_in = gijsub_out

         dtCFLsub = compute_CFL_opt(eps,nbins,kpol,massgrid,massbins,gijsub_in,tabflux_coag)
         dtCFLsub = coeff_CFL*dtCFLsub


      enddo

      !last timestep to reach dthydro
      dtlast = dthydro-dtsub

      nsub = nsub + 1

      call time_solver_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gijsub_in,tabflux_coag,tabintflux_coag,dtlast,gijnew)

   !when coagulation CFL > hydro timestep
   else

      ndt = ndt + 1

      call time_solver_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tabflux_coag,tabintflux_coag,dt,gijnew)

   endif


end subroutine compute_coag_opt


!> @brief Compute coagulation solver (physical collision kernel) for 1 hydro time-step
!!
!! with indices optimisation (only with gij > eps)
!!
!! DG scheme k>0, piecewise polynomial approximation
!!
!! @param[in]    eps               minimum value for mass distribution approximation gij
!! @param[in]    coeff_CFL         timestep coefficient for stability of the SSPRK order 3 scheme
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    massbins          arithmetic mean value of massgrid for each mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[in]    tabflux_coag      5D array to evaluate coagulation flux
!! @param[in]    tabintflux_coag   6D array to evaluate the term including the integral of coagulation flux
!! @param[in]    dv                2D array of the differential velocity between grains in bins lp and l
!! @param[in]    dthydro           hydro timestep
!! @param[in]    gij               components of g on the polynomial basis
!! @param[out]   gijnew            evolved components of g on the polynomial basis after dthydro
!! @param[out]   nsub              number of subcycling coagulation timestep to reach dthydro
!! @param[out]   ndt               number of hydro timestep, when coagulation CFL > dthydro
subroutine compute_coag_kdv_opt(eps,coeff_CFL,nbins,kpol,massgrid,massbins,mat_coeffs_leg,&
                                    tabflux_coag,tabintflux_coag,dv,&
                                    dthydro,gij,gijnew,nsub,ndt)
   implicit none
   integer,  intent(in)  :: nbins,kpol
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins),mat_coeffs_leg(kpol+1,kpol+1)
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   real(wp), intent(in)  :: tabintflux_coag(nbins,kpol+1,nbins,nbins,kpol+1,kpol+1)
   real(wp), intent(in)  :: dv(nbins,nbins)
   real(wp), intent(in)  :: dthydro,coeff_CFL,eps
   real(wp), intent(in)  :: gij(nbins,kpol+1)
   integer,  intent(out) :: nsub,ndt
   real(wp), intent(out) :: gijnew(nbins,kpol+1)

   real(wp) :: dtCFLsub,dtsub,dtlast,dt,gijsub_in(nbins,kpol+1),gijsub_out(nbins,kpol+1)

   nsub = 0
   ndt = 0

   !evaluate coagulation CFL
   dtCFLsub = compute_CFL_kdv_opt(eps,nbins,kpol,massgrid,massbins,gij,tabflux_coag,dv)
   dtCFLsub = coeff_CFL*dtCFLsub

   !compare hydro timestep and coagulation CFL
   dt = min(dtCFLsub,dthydro)

   !coagulation subcycling timesteps
   if (dt<dthydro) then
      gijsub_in = gij
      dtsub = 0._wp
      do while ( dtsub<dthydro .and. dthydro-dtsub>dtCFLsub)
         dtsub=dtsub+dtCFLsub

         nsub = nsub + 1

         call time_solver_kdv_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gijsub_in,&
                                 tabflux_coag,tabintflux_coag,dv,dtCFLsub,gijsub_out)
         gijsub_in = gijsub_out

         dtCFLsub = compute_CFL_kdv_opt(eps,nbins,kpol,massgrid,massbins,gijsub_in,tabflux_coag,dv)
         dtCFLsub = coeff_CFL*dtCFLsub


      enddo

      !last timestep to reach dthydro
      dtlast = dthydro-dtsub

      nsub = nsub + 1

      call time_solver_kdv_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gijsub_in,&
                              tabflux_coag,tabintflux_coag,dv,dtlast,gijnew)

   !when coagulation CFL > hydro timestep
   else

      ndt = ndt + 1

      call time_solver_kdv_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tabflux_coag,tabintflux_coag,dv,dt,gijnew)

   endif


end subroutine compute_coag_kdv_opt

end module compute_solver_coag





