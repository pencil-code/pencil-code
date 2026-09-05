!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************

!-------------------------------------------
! MODULE: compute intflux and flux 
!-------------------------------------------

module functions_flux_intflux
   use precision
   implicit none

contains

!----------------------------------------------
! Original version to compute flux and intflux
!----------------------------------------------


!> @brief 2D array with element g_lp x g_l, for DG scheme k=0 for simple kernels
!!
!! @param[in]    nbins     number of dust bins
!! @param[in]    gij       components of g on the polynomial basis
!! @param[out]   mat_gij   matrix with element g_lp x g_l
subroutine compute_arr_gij_k0(nbins,gij,arr_gij)
   implicit none
   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: gij(nbins)
   real(wp), intent(out) :: arr_gij(nbins,nbins) 

   integer :: lp,l


   do lp=1,nbins
      do l=1,nbins
         arr_gij(lp,l) = gij(lp)*gij(l)
      enddo
   enddo

end subroutine compute_arr_gij_k0


!> @brief 2D array with element g_lp x g_l x dv(lp,l), for DG scheme k=0 with physical kernel
!!
!! with indices optimisation (only with gij > eps)
!!
!! @param[in]    nbins     number of dust bins
!! @param[in]    gij       components of g on the polynomial basis
!! @param[in]    dv        2D array of the differential velocity between grains in bins lp and l
!! @param[out]   arr_gij   matrix with element g_lp x g_l x dv(lp,l)
subroutine compute_arr_gij_k0_kdv(nbins,gij,dv,arr_gij_dv)
   implicit none
   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: gij(nbins),dv(nbins,nbins)
   real(wp), intent(out) :: arr_gij_dv(nbins,nbins) 

   integer :: lp,l

   do lp=1,nbins
      do l=1,nbins
         arr_gij_dv(lp,l) = gij(lp)*gij(l)*dv(lp,l)
      enddo
   enddo

end subroutine compute_arr_gij_k0_kdv


!> @brief compute the approximation of the coagulation flux with DG scheme k=0
!!
!! Coagulation flux is defined at the right boundary of mass bins, i.e. flux(j) ~ F(m_{j+1/2})
!!
!! @param[in]    nbins          number of dust bins
!! @param[in]    arr_gij        matrix with element g_lp x g_l
!! @param[in]    tabflux_coag   3D array to evaluate coagulation flux
!! @param[out]   flux           approximation of the coagulation flux in each bin
subroutine compute_flux_k0(nbins,arr_gij,tabflux_coag,flux)
   implicit none
   integer,  intent(in)   :: nbins
   real(wp), intent(in)   :: arr_gij(nbins,nbins)
   real(wp), intent(in)   :: tabflux_coag(nbins,nbins,nbins)
   real(wp), intent(out)  :: flux(nbins)

   integer  :: j

   flux = 0._wp

   do j=1,nbins-1
      flux(j) = sum(tabflux_coag(j,:,:)*arr_gij(:,:))

   enddo

end subroutine compute_flux_k0



!> @brief 4D array with element g_{lp}^{ip} x g_{l}^{i}, for DG scheme k>0, for simple kernels
!!
!! @param[in]    nbins     number of dust bins
!! @param[in]    kpol      degree of polynomials for approximation
!! @param[in]    gij       components of g on the polynomial basis
!! @param[out]   arr_gij   matrix with element g_{lp}^{ip} x g_{l}^{i}
subroutine compute_arr_gij(nbins,kpol,gij,arr_gij)
   implicit none
   integer,  intent(in)  :: nbins,kpol
   real(wp), intent(in)  :: gij(nbins,kpol+1)
   real(wp), intent(out) :: arr_gij(nbins,nbins,kpol+1,kpol+1) 

   integer :: lp,l,ip,i

   do ip=0,kpol
      do i=0,kpol
         do lp=1,nbins
            do l=1,nbins
               arr_gij(lp,l,ip+1,i+1) = gij(lp,ip+1)*gij(l,i+1)
            enddo
         enddo
      enddo
   enddo

end subroutine compute_arr_gij





!> @brief 4D array with element g_{lp}^{ip} x g_{l}^{i}, for DG scheme k>0, for physical kernel
!!
!! @param[in]    nbins     number of dust bins
!! @param[in]    kpol      degree of polynomials for approximation
!! @param[in]    gij       components of g on the polynomial basis
!! @param[in]    dv        2D array of the differential velocity between grains in bins lp and l
!! @param[out]   arr_gij   4D array with element g_{lp}^{ip} x g_{l}^{i}
subroutine compute_arr_gij_kdv(nbins,kpol,gij,dv,arr_gij_dv)

   implicit none
   integer,  intent(in)  :: nbins,kpol
   real(wp), intent(in)  :: gij(nbins,kpol+1),dv(nbins,nbins)
   real(wp), intent(out) :: arr_gij_dv(nbins,nbins,kpol+1,kpol+1) 

   integer :: lp,l,ip,i


   do ip=0,kpol
      do i=0,kpol
         do lp=1,nbins
            do l=1,nbins
               arr_gij_dv(lp,l,ip+1,i+1) = gij(lp,ip+1)*gij(l,i+1)*dv(lp,l)
            enddo
         enddo
      enddo
   enddo

end subroutine compute_arr_gij_kdv


!> @brief compute the approximation of the coagulation flux with DG scheme k>0
!!
!! Coagulation flux is defined at the right boundary of mass bins, i.e. flux[j] ~ F(m_{j+1/2})
!!
!! @param[in]    nbins          number of dust bins
!! @param[in]    kpol           degree of polynomials for approximation
!! @param[in]    arr_gij        4D array with element g_{lp}^{ip} x g_{l}^{i}
!! @param[in]    tabflux_coag   3D array to evaluate coagulation flux
!! @param[out]   flux           approximation of the coagulation flux in each bin
subroutine compute_flux(nbins,kpol,arr_gij,tabflux_coag,flux)
   use omp_lib
   implicit none
   integer,  intent(in)  :: nbins,kpol
   real(wp), intent(in)  :: arr_gij(nbins,nbins,kpol+1,kpol+1)
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   real(wp), intent(out) :: flux(nbins)

   integer :: j

   flux = 0._wp

   do j=1,nbins-1
      flux(j) =  sum(tabflux_coag(j,:,:,:,:)*arr_gij(:,:,:,:))

   enddo

end subroutine compute_flux





!> @brief compute the approximation of the coagulation flux and the term including the integral of the flux with DG scheme k>0
!!
!! Coagulation flux is defined at the right boundary of mass bins, i.e. flux[j] ~ F(m_{j+1/2})
!!
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    arr_gij           4D array with element g_{lp}^{ip} x g_{l}^{i}
!! @param[in]    tabflux_coag      3D array to evaluate coagulation flux
!! @param[in]    tabintflux_coag   5D array to evaluate the term including the integral of the coagulation flux
!! @param[out]   flux              approximation of the coagulation flux in each bin
!! @param[out]   intflux           approximation of the term including the integral of coagulation flux at each bin
subroutine compute_flux_intflux(nbins,kpol,arr_gij,tabflux_coag,tabintflux_coag,flux,intflux)
   implicit none
   integer, intent(in)  :: nbins,kpol
   real(wp),intent(in)  :: arr_gij(nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(in)  :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(in)  :: tabintflux_coag(nbins,kpol+1,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(out) :: flux(nbins)
   real(wp),intent(out) :: intflux(nbins,kpol+1)

   integer :: j,k

   flux = 0._wp
   intflux = 0._wp

   do j=1,nbins
      if (j<nbins) then
         flux(j) = sum(tabflux_coag(j,:,:,:,:)*arr_gij(:,:,:,:))

      endif
      do k=1,kpol
         intflux(j,k+1) = sum(tabintflux_coag(j,k+1,:,:,:,:)*arr_gij(:,:,:,:)) 
      enddo
   enddo

end subroutine compute_flux_intflux






!----------------------------------------------
! Optimised version to compute flux and intflux
!----------------------------------------------

!> @brief Minimum and maximum indices of bins defining mass distribution, where gij > eps (minimum value)
!!
!! When gij < eps, grains are not considered
!!
!! @param[in]    eps       minimum value for mass distribution approximation gij
!! @param[in]    nbins     number of dust bins
!! @param[in]    massbins  arithmetic mean value of massgrid for each mass bins
!! @param[in]    gij       components of g on the polynomial basis
!! @param[out]   ind_min   minimum index for bin where gij > eps
!! @param[out]   ind_max   maximum index for bin where gij > eps
subroutine compute_min_max_indices(eps,nbins,massbins,gij,ind_min,ind_max)
   implicit none
   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: eps,massbins(nbins),gij(nbins)
   integer,  intent(out) :: ind_min,ind_max


   !minloc and maxloc need sorted values array to find first index 
   !corresponding min index and max index for which gij > eps
   ind_min = minloc(massbins, mask = gij > eps, dim = 1)
   ind_max = maxloc(massbins, mask = gij > eps, dim = 1)

end subroutine compute_min_max_indices

!> @brief 2D array with element g_lp x g_l, for DG scheme k=0 for simple kernels
!!
!! with indices optimisation (only with gij > eps)
!!
!! @param[in]    nbins     number of dust bins
!! @param[in]    gij       components of g on the polynomial basis
!! @param[in]    ind_min   smallest index for bins where gij > eps
!! @param[in]    ind_max   largest index for bins where gij > eps
!! @param[in]    gij       components of g on the polynomial basis
!! @param[out]   arr_gij   array with element g_lp x g_l
subroutine compute_arr_gij_k0_opt(nbins,gij,ind_min,ind_max,arr_gij)
   implicit none
   integer,  intent(in)  :: nbins,ind_min,ind_max
   real(wp), intent(in)  :: gij(nbins)
   real(wp), intent(out) :: arr_gij(nbins,nbins) 

   integer :: lp,l

   arr_gij = 0._wp
   do lp=ind_min,ind_max
      do l=ind_min,ind_max
         arr_gij(lp,l) = gij(lp)*gij(l)
      enddo
   enddo

end subroutine compute_arr_gij_k0_opt


!> @brief 2D array with element g_lp x g_l x dv(lp,l), for DG scheme k=0 with physical kernel
!!
!! @param[in]    nbins     number of dust bins
!! @param[in]    gij       components of g on the polynomial basis
!! @param[in]    ind_min   smallest index for bins where gij > eps
!! @param[in]    ind_max   largest index for bins where gij > eps
!! @param[in]    dv        2D array of the differential velocity between grains in bins lp and l
!! @param[out]   arr_gij   matrix with element g_lp x g_l x dv(lp,l)
subroutine compute_arr_gij_k0_kdv_opt(nbins,gij,ind_min,ind_max,dv,arr_gij_dv)
   implicit none
   integer,  intent(in)  :: nbins,ind_min,ind_max
   real(wp), intent(in)  :: gij(nbins),dv(nbins,nbins)
   real(wp), intent(out) :: arr_gij_dv(nbins,nbins) 

   integer :: lp,l

   arr_gij_dv = 0._wp
   do lp=ind_min,ind_max
      do l=ind_min,ind_max
         arr_gij_dv(lp,l) = gij(lp)*gij(l)*dv(lp,l)
      enddo
   enddo

end subroutine compute_arr_gij_k0_kdv_opt


!> @brief compute the approximation of the coagulation flux with DG scheme k=0 (simple kernels)
!!
!! with indices optimisation (only with gij > eps)
!!
!! Coagulation flux is defined at the right boundary of mass bins, i.e. flux[j] ~ F(m_{j+1/2})
!!
!! @param[in]    nbins          number of dust bins
!! @param[in]    ind_min        smallest index for bins where gij > eps
!! @param[in]    arr_gij        array with element g_lp x g_l
!! @param[in]    tabflux_coag   3D array to evaluate coagulation flux
!! @param[out]   flux           approximation of the coagulation flux in each bin
subroutine compute_flux_k0_opt(nbins,ind_min,arr_gij,tabflux_coag,flux)
   implicit none
   integer,  intent(in)   :: nbins,ind_min
   real(wp), intent(in)   :: arr_gij(nbins,nbins)
   real(wp), intent(in)   :: tabflux_coag(nbins,nbins,nbins)
   real(wp), intent(out)  :: flux(nbins)

   integer  :: j

   flux = 0._wp

   do j=ind_min,nbins-1
   
      flux(j) = sum(tabflux_coag(j,:,:)*arr_gij(:,:))

   enddo

end subroutine compute_flux_k0_opt


!> @brief 4D array with element g_{lp}^{ip} x g_{l}^{i}, for DG scheme k>0 (simple kernels), with indices optimisation (where gij > eps) 
!!
!! @param[in]    nbins     number of dust bins
!! @param[in]    kpol      degree of polynomials for approximation
!! @param[in]    gij       components of g on the polynomial basis
!! @param[in]    ind_min   smallest index for bins where gij(:,1) > eps
!! @param[in]    ind_max   largest index for bins where gij(:,1) > eps
!! @param[out]   arr_gij   matrix with element g_{lp}^{ip} x g_{l}^{i}
subroutine compute_arr_gij_opt(nbins,kpol,gij,ind_min,ind_max,arr_gij)

   implicit none
   integer,  intent(in)  :: nbins,kpol,ind_min,ind_max
   real(wp), intent(in)  :: gij(nbins,kpol+1)
   real(wp), intent(out) :: arr_gij(nbins,nbins,kpol+1,kpol+1) 

   integer :: lp,l,ip,i


   arr_gij = 0._wp

   do ip=0,kpol
      do i=0,kpol
         do lp=ind_min,ind_max
            do l=ind_min,ind_max
               arr_gij(lp,l,ip+1,i+1) = gij(lp,ip+1)*gij(l,i+1)
            enddo
         enddo
      enddo
   enddo

end subroutine compute_arr_gij_opt


!> @brief 4D array with element g_{lp}^{ip} x g_{l}^{i}, for DG scheme k>0, for physical kernel
!! with indices optimisation (where gij > eps)
!!
!! @param[in]    nbins     number of dust bins
!! @param[in]    kpol      degree of polynomials for approximation
!! @param[in]    gij       components of g on the polynomial basis
!! @param[in]    ind_min   smallest index for bins where gij(:,1) > eps
!! @param[in]    ind_max   largest index for bins where gij(:,1) > eps
!! @param[in]    dv        2D array of the differential velocity between grains in bins lp and l
!! @param[out]   arr_gij   4D array with element g_{lp}^{ip} x g_{l}^{i}
subroutine compute_arr_gij_kdv_opt(nbins,kpol,gij,ind_min,ind_max,dv,arr_gij_dv)

   implicit none
   integer,  intent(in)  :: nbins,kpol,ind_min,ind_max
   real(wp), intent(in)  :: gij(nbins,kpol+1),dv(nbins,nbins)
   real(wp), intent(out) :: arr_gij_dv(nbins,nbins,kpol+1,kpol+1) 

   integer :: lp,l,ip,i

   do ip=0,kpol
      do i=0,kpol
         do lp=ind_min,ind_max
            do l=ind_min,ind_max
               arr_gij_dv(lp,l,ip+1,i+1) = gij(lp,ip+1)*gij(l,i+1)*dv(lp,l)
            enddo
         enddo
      enddo
   enddo

end subroutine compute_arr_gij_kdv_opt


!> @brief compute the approximation of the coagulation flux with DG scheme k>0
!!        
!! with indices optimisation (where gij > eps)
!!
!! Coagulation flux is defined at the right boundary of mass bins, i.e. flux[j] ~ F(m_{j+1/2})
!!
!! @param[in]    nbins          number of dust bins
!! @param[in]    kpol           degree of polynomials for approximation
!! @param[in]    ind_min        smallest index for bins where gij > eps
!! @param[in]    arr_gij   matrix with element g_{lp}^{ip} x g_{l}^{i}
!! @param[in]    tabflux_coag   3D array to evaluate coagulation flux
!! @param[out]   flux           approximation of the coagulation flux in each bin
subroutine compute_flux_opt(nbins,kpol,ind_min,arr_gij,tabflux_coag,flux)
   implicit none
   integer,  intent(in)  :: nbins,kpol,ind_min
   real(wp), intent(in)  :: arr_gij(nbins,nbins,kpol+1,kpol+1) 
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   real(wp), intent(out) :: flux(nbins)

   integer  :: j

   flux = 0._wp


   do j=ind_min,nbins-1
      flux(j) =  sum(tabflux_coag(j,:,:,:,:)*arr_gij(:,:,:,:))
   enddo

end subroutine compute_flux_opt


!> @brief compute the approximation of the coagulation flux and the term including the integral of the flux with DG scheme k>0 (simple kernels),
!!
!! with indices optimisation (only with gij > eps)
!!
!! Coagulation flux is defined at the right boundary of mass bins, i.e. flux[j] ~ F(m_{j+1/2})
!!
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    ind_min           smallest index for bins where gij > eps
!! @param[in]    gij               components of g on the polynomial basis
!! @param[in]    tabflux_coag      3D array to evaluate coagulation flux
!! @param[in]    tabintflux_coag   5D array to evaluate the term including the integral of the coagulation flux
!! @param[out]   flux              approximation of the coagulation flux in each bin
!! @param[out]   intflux           approximation of the term including the integral of coagulation flux at each bin
subroutine compute_flux_intflux_opt(nbins,kpol,ind_min,arr_gij,tabflux_coag,tabintflux_coag,flux,intflux)
   implicit none
   integer, intent(in)  :: nbins,kpol,ind_min
   real(wp),intent(in)  :: arr_gij(nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(in)  :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(in)  :: tabintflux_coag(nbins,kpol+1,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(out) :: flux(nbins)
   real(wp),intent(out) :: intflux(nbins,kpol+1)

   integer  :: j,k

   flux = 0._wp
   intflux = 0._wp


   do j=ind_min,nbins
      if (j<nbins) then
         flux(j) = sum(tabflux_coag(j,:,:,:,:)*arr_gij(:,:,:,:))

      endif
      do k=1,kpol
         intflux(j,k+1) = sum(tabintflux_coag(j,k+1,:,:,:,:)*arr_gij(:,:,:,:)) 
      enddo
   enddo


end subroutine compute_flux_intflux_opt

end module functions_flux_intflux