!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


!-----------------------------------------------
! MODULE: Limiter for DG scheme (Zhang,Shu 2010)
!-----------------------------------------------

module limiter
   use precision
   use polynomials_legendre
   use reconstruction_g
   use functions_flux_intflux
   implicit none

contains




!> @brief Minimum value of the approximation of g in each bin
!!
!! @param[in]    nbins                number of dust bins
!! @param[in]    kpol                 degree of polynomials for approximation
!! @param[in]    massgrid             grid of masses, borders value of mass bins
!! @param[in]    massbins             arithmetic mean value of massgrid for each mass bins
!! @param[in]    mat_coeffs_leg       array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[in]    gij                  components of g on the polynomial basis
!! @param[in]    j                    indice of bin
!! @param[in]    x                    value to evalute the polynomial
!! @param[out]   tabminval_approx_g   reconstruction of g in bin j evaluated in x
subroutine minval_approx_g(nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tabminval_approx_g)
   implicit none
   integer,  intent(in)  :: kpol,nbins
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins),mat_coeffs_leg(kpol+1,kpol+1)
   real(wp), intent(in)  :: gij(nbins,kpol+1)
   real(wp), intent(out) :: tabminval_approx_g(nbins)

   real(wp) :: xjgridl,xjgridr,dx,xval,g_left,g_right,res
   real(wp), allocatable :: func_pol(:)
   integer  :: j,npoints,i


   tabminval_approx_g = 0._wp

   if (kpol == 1) then
   
      do j=1,nbins
         call recons_g(nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,j,massgrid(j),g_left)
         call recons_g(nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,j,massgrid(j+1),g_right)
         tabminval_approx_g(j) = min(g_left,g_right)
      enddo

   else if (kpol > 1) then


      npoints = 200
      allocate(func_pol(npoints))
      func_pol = 0._wp
      do j=1,nbins
         xjgridl = massgrid(j)
         xjgridr = massgrid(j+1)
         dx = (xjgridr-xjgridl)/real(npoints-1,wp)

         
         do i=1,npoints
            xval = xjgridl + (i-1)*dx
            call recons_g(nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,j,xval,res)
            func_pol(i)  = res
            
         enddo

         tabminval_approx_g(j) = minval(func_pol)

      enddo

      deallocate(func_pol)

   else
      stop "Need correct kpol for minvalpol function in minval_approx_g subroutine"
   endif   

end subroutine minval_approx_g



!> @brief Limiter coefficient to ensure positivity of the numerical solution (Zhang and Shu 2010)
!!
!! @param[in]    eps               minimum value for mass distribution approximation gij
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    massbins          arithmetic mean value of massgrid for each mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[in]    gij               components of g on the polynomial basis
!! @param[out]   tab_gamma         limiter coefficient in each bin
subroutine gammafunction(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tab_gamma)
   implicit none
   integer,  intent(in)  :: kpol,nbins
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins),eps,mat_coeffs_leg(kpol+1,kpol+1)
   real(wp), intent(in)  :: gij(nbins,kpol+1)
   real(wp), intent(out) :: tab_gamma(nbins)

   real(wp) :: tabminval_approx_g(nbins)
   integer  :: j

   tab_gamma = 0._wp

   !Liu et al. 2019
   if(kpol==0) then 
      tab_gamma = 1._wp
   else
      
      call minval_approx_g(nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tabminval_approx_g)

      
      do j=1,nbins
        if (gij(j,1) == tabminval_approx_g(j)) then    
            tab_gamma(j) = 1._wp
        else
            tab_gamma(j) = min(1._wp,abs((gij(j,1)-eps)/(gij(j,1)-tabminval_approx_g(j))))
        endif

           
      enddo
   endif

end subroutine gammafunction


!> @brief Limiter coefficient to ensure positivity of the numerical solution (Zhang and Shu 2010), with indices optimisation (only with gij > eps) 
!!
!! @param[in]    eps               minimum value for mass distribution approximation gij
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    massbins          arithmetic mean value of massgrid for each mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[in]    gij               components of g on the polynomial basis
!! @param[out]   tab_gamma         limiter coefficient in each bin
subroutine gammafunction_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tab_gamma)
   implicit none
   integer,  intent(in)  :: kpol,nbins
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins),eps,mat_coeffs_leg(kpol+1,kpol+1)
   real(wp), intent(in)  :: gij(nbins,kpol+1)
   real(wp), intent(out) :: tab_gamma(nbins)

   real(wp) :: tabminval_approx_g(nbins)
   integer  :: j,ind_min,ind_max


   call compute_min_max_indices(eps,nbins,massbins,gij(:,1),ind_min,ind_max)

   tab_gamma = 0._wp
      
   call minval_approx_g(nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tabminval_approx_g)

   do j=ind_min,min(ind_max,nbins)
        
      tab_gamma(j) = min(1._wp,abs((gij(j,1)-eps)/(gij(j,1)-tabminval_approx_g(j))))
       
   enddo

end subroutine gammafunction_opt

end module limiter
