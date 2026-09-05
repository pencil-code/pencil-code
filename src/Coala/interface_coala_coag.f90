!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


!-------------------------------------------
! MODULE: Function to interface COALA with other code
!-------------------------------------------

module interface_coala_coag
   use precision
   use initial_condition_L2proj_legbasis_GQ
   use limiter
   use compute_solver_coag

   implicit none
contains


!> @brief Compute coagulation solver in hydro code for 1 hydro time-step 
!!
!! DG scheme k=0, piecewise constant approximation
!! Function for interface with hydro code
!!
!! @param[in]    nbins         number of dust bins
!! @param[in]    massgrid      grid of masses, borders value of mass bins
!! @param[in]    tabflux_coag  3D array to evaluate coagulation flux
!! @param[in]    rhodust       1D array dust density for each grain size
!! @param[in]    eps_rhodust   minimum value for dust density
!! @param[in]    dv            2D array of the differential velocity between grains in bins lp and l
!! @param[in]    dthydro       hydro timestep
!! @param[out]   new_rhodust   evolved 1D array dust density for each grain size
subroutine coala_coag_k0(nbins,massgrid,tabflux_coag,rhodust,eps_rhodust,dv,dthydro,new_rhodust)
   implicit none
   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: massgrid(nbins+1)
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins),dv(nbins,nbins)
   real(wp), intent(in)  :: dthydro,eps_rhodust
   real(wp), intent(in)  :: rhodust(nbins)
   real(wp), intent(out) :: new_rhodust(nbins)

   integer  :: j,iprogress,nsub,ndt,tot_nsub,tot_ndt
   real(wp) :: gij(nbins),gijnew(nbins),coeff_CFL,eps_gij
   real(wp) :: massbins(nbins)

   massbins = 0.5_wp * (massgrid(2:)+massgrid(1:nbins))

   eps_gij =  eps_rhodust/massgrid(nbins+1)
   ! print*,"eps_gij=",eps_gij

   !rhodust -> gij
   gij = eps_gij
   do j=1,nbins
      if (rhodust(j) > eps_rhodust) then
         gij(j) = rhodust(j)/(massgrid(j+1)-massgrid(j))
      endif
   enddo

  
   iprogress=1
   tot_nsub = 0
   tot_ndt = 0

   coeff_CFL = 3e-1_wp

   !Solve coagulation for 1 dthydro
   ! optimised version using only values rhodust > eps_rhodust
   ! call compute_coag_k0_kdv_opt(eps_gij,coeff_CFL,nbins,massgrid,massbins,&
   !                                  tabflux_coag,dv,dthydro,gij,gijnew,nsub,ndt)
   
   !for debug use this non optimised function
   call compute_coag_k0_kdv(eps_gij,coeff_CFL,nbins,massgrid,tabflux_coag,dv,dthydro,gij,gijnew,nsub,ndt)
   
   tot_nsub = tot_nsub + nsub
   tot_ndt = tot_ndt + ndt


   !gijnew -> new_rhodust
   new_rhodust = 0._wp
   do j=1,nbins
      new_rhodust(j) = max(eps_rhodust,gijnew(j)*(massgrid(j+1)-massgrid(j)))
   enddo

end subroutine coala_coag_k0



!> @brief Compute coagulation solver in hydro code for 1 hydro time-step 
!!
!! DG scheme k>0, piecewise polynomial approximation
!! Function for interface with hydro code
!!
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    massbins          arithmetic mean value of massgrid for each mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[in]    Q                 number of points for Gauss-Legendre quadrature
!! @param[in]    vecnodes          nodes of the Legendre polynomials
!! @param[in]    vecweights        weights coefficients for the Gauss-Legendre polynomials
!! @param[in]    tabflux_coag      5D array to evaluate coagulation flux
!! @param[in]    tabintflux_coag   6D array to evaluate the term including the integral of coagulation flux
!! @param[in]    rhodust           1D array dust density for each grain size
!! @param[in]    eps_rhodust       minimum value for dust density
!! @param[in]    dv                2D array of the differential velocity between grains in bins lp and l
!! @param[in]    dthydro           hydro timestep
!! @param[out]   new_rhodust       evolved 1D array dust density for each grain size
subroutine coala_coag(nbins,kpol,massgrid,massbins,mat_coeffs_leg,Q,vecnodes,vecweights,tabflux_coag,tabintflux_coag,&
                        rhodust,eps_rhodust,dv,dthydro,new_rhodust)
   implicit none
   integer,  intent(in)  :: nbins,kpol,Q
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins),mat_coeffs_leg(kpol+1,kpol+1)
   real(wp), intent(in)  :: vecnodes(Q),vecweights(Q)
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   real(wp), intent(in)  :: tabintflux_coag(nbins,kpol+1,nbins,nbins,kpol+1,kpol+1)
   real(wp), intent(in)  :: dv(nbins,nbins)
   real(wp), intent(in)  :: dthydro,eps_rhodust
   real(wp), intent(in)  :: rhodust(nbins)
   real(wp), intent(out) :: new_rhodust(nbins)

   integer  :: j,k,iprogress,nsub,ndt,tot_nsub,tot_ndt
   real(wp) :: gij(nbins,kpol+1),gijnew(nbins,kpol+1),tab_gamma(nbins),coeff_CFL,eps_gij

   eps_gij = eps_rhodust/massgrid(nbins+1)

   !rhodust -> gij + enforcing mass conservation after interpolation
   gij = 0._wp
   call L2proj_gij_GQ(nbins,kpol,massgrid,massbins,mat_coeffs_leg,Q,vecnodes,vecweights,eps_rhodust,rhodust,gij)

   !positivity of the polynomials approx
   !apply limiter coefficient
   call gammafunction(eps_gij,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tab_gamma)
   do k=1,kpol
      gij(:,k+1) = tab_gamma(:)*gij(:,k+1)
   enddo

   !limit to eps value
   do j=1,nbins
      if (gij(j,1) < 0._wp) then
         print*,"j=",j,", gij =",gij(j,1)
         stop
      else if ( gij(j,1) <= eps_gij  ) then
         gij(j,1) = eps_gij
         gij(j,2:) = 0._wp

      endif
   enddo


   !variable for analysis
   iprogress=1
   tot_nsub = 0
   tot_ndt = 0

   !CFL coefficient for SSPRK 3
   coeff_CFL = 3e-1_wp

   !Solve coagulation for 1 dthydro
   ! optimised version using only values rhodust > eps_rhodust
   call compute_coag_kdv_opt(eps_gij,coeff_CFL,nbins,kpol,massgrid,massbins,mat_coeffs_leg,&
                              tabflux_coag,tabintflux_coag,dv,&
                              dthydro,gij,gijnew,nsub,ndt)

   !for debug use this non optimised function
   ! call compute_coag_kdv(eps_gij,coeff_CFL,nbins,kpol,massgrid,massbins,mat_coeffs_leg,&
   !                            tabflux_coag,tabintflux_coag,dv,&
   !                            dthydro,gij,gijnew,nsub,ndt)

   
   tot_nsub = tot_nsub + nsub
   tot_ndt = tot_ndt + ndt

   

   !gij -> rhodust
   new_rhodust = 0._wp
   do j=1,nbins
      new_rhodust(j) = max(eps_rhodust,gij(j,1)*(massgrid(j+1)-massgrid(j)))
   enddo

   
end subroutine coala_coag


end module interface_coala_coag





