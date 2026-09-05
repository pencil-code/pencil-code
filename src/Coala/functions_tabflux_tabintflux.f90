!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


!-------------------------------------------
! MODULE: compute tabflux and tabintflux
!-------------------------------------------

module functions_tabflux_tabintflux
   use precision
   use progress_bar
   use coagflux_functions_GQ
   use coagintflux_functions_GQ

   implicit none

contains

!> @brief Precompute array depending only on massgrid to evaluate the coagulation flux
!!
!! DG scheme with piecewise constant approximation
!!
!! @param[in]    kernel            select the collisional kernel function
!! @param[in]    K0                constant value of the kernel function (used to adapt to code unit)
!! @param[in]    Q                 number of points for Gauss-Legendre quadrature
!! @param[in]    vecnodes          nodes of the Legendre polynomials
!! @param[in]    vecweights        weights coefficients for the Gauss-Legendre polynomials
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[out]   tabflux           3D array to evaluate coagulation flux
subroutine compute_coagtabflux_GQ_k0(kernel,K0,Q,vecnodes,vecweights,nbins,kpol,massgrid,mat_coeffs_leg,tabflux)
   implicit none
   character(len=*), intent(in)  :: kernel
   integer,          intent(in)  :: Q,nbins,kpol
   real(wp),         intent(in)  :: vecnodes(Q),vecweights(Q),K0,mat_coeffs_leg(kpol+1,kpol+1)
   real(wp),         intent(in)  :: massgrid(nbins+1)
   real(wp),         intent(out) :: tabflux(nbins,nbins,nbins)

   integer  :: j,lp,l,iprogress
   real(wp) :: res

   tabflux = 0._wp

   iprogress = 1
   do j = 1,nbins
      do lp = 1,j
         do l = 1,nbins
            res = coagfluxfunction_GQ(kernel,K0,Q,vecnodes,vecweights,nbins,kpol,massgrid,mat_coeffs_leg,j,lp,l,0,0)
   
            if (res /= res) then
               print*,"NAN in coagfluxfunction_GQ for k=",kpol; stop
            endif
            tabflux(j,lp,l) = res

         enddo
      enddo
      call display_progress_bar(nbins,iprogress)
      iprogress = iprogress+1
   enddo
   print*,""


end subroutine compute_coagtabflux_GQ_k0


!> @brief Precompute array depending only on massgrid to evaluate the coagulation flux
!!
!! DG scheme with piecewise polynomial approximation
!!
!! @param[in]    kernel            select the collisional kernel function
!! @param[in]    K0                constant value of the kernel function (used to adapt to code unit)
!! @param[in]    Q                 number of points for Gauss-Legendre quadrature
!! @param[in]    vecnodes          nodes of the Legendre polynomials
!! @param[in]    vecweights        weights coefficients for the Gauss-Legendre polynomials
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[in]    tabflux           5D array to evaluate coagulation flux
subroutine compute_coagtabflux_GQ(kernel,K0,Q,vecnodes,vecweights,nbins,kpol,massgrid,mat_coeffs_leg,tabflux)
   implicit none
   character(len=*), intent(in)  :: kernel
   integer,          intent(in)  :: Q,nbins,kpol
   real(wp),         intent(in)  :: vecnodes(Q),vecweights(Q),K0,mat_coeffs_leg(kpol+1,kpol+1)
   real(wp),         intent(in)  :: massgrid(nbins+1)
   real(wp),         intent(out) :: tabflux(nbins,nbins,nbins,kpol+1,kpol+1)

   integer  :: j,lp,l,ip,i,iprogress
   real(wp) :: res

   tabflux = 0._wp

   iprogress = 1
   do j = 1,nbins
      do lp = 1,j
         do l = 1,nbins
            do ip= 0,kpol
               do i= 0,kpol
                  res = coagfluxfunction_GQ(kernel,K0,Q,vecnodes,vecweights,nbins,kpol,massgrid,mat_coeffs_leg,j,lp,l,ip,i)
         
                  if (res /= res) then
                     print*,"NAN in coagfluxfunction_GQ for k=",kpol; stop
                  endif
                  tabflux(j,lp,l,ip+1,i+1) = res

               enddo
            enddo
         enddo
      enddo
      call display_progress_bar(nbins,iprogress)
      iprogress = iprogress+1
   enddo
   print*,""

end subroutine compute_coagtabflux_GQ


!> @brief Precompute array depending only on massgrid to evaluate the term including integral of the coagulation flux
!!
!! DG scheme with piecewise polynomial approximation
!!
!! @param[in]    kernel            select the collisional kernel function
!! @param[in]    K0                constant value of the kernel function (used to adapt to code unit)
!! @param[in]    Q                 number of points for Gauss-Legendre quadrature
!! @param[in]    vecnodes          nodes of the Legendre polynomials
!! @param[in]    vecweights        weights coefficients for the Gauss-Legendre polynomials
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[out]   tabintflux        6D array to evaluate coagulation flux
subroutine compute_coagtabintflux_GQ(kernel,K0,Q,vecnodes,vecweights,nbins,kpol,massgrid,mat_coeffs_leg,tabintflux)
   implicit none
   character(len=*), intent(in)  :: kernel
   integer,          intent(in)  :: Q,nbins,kpol
   real(wp),         intent(in)  :: vecnodes(Q),vecweights(Q),K0,mat_coeffs_leg(kpol+1,kpol+1)
   real(wp),         intent(in)  :: massgrid(nbins+1)
   real(wp),         intent(out) :: tabintflux(nbins,kpol+1,nbins,nbins,kpol+1,kpol+1)

   integer  :: j,k,lp,l,ip,i,iprogress
   real(wp) :: res

   tabintflux = 0._wp

   iprogress = 1


   do j = 1,nbins
      do k=1,kpol
         do lp = 1,j
            do l = 1,nbins
               do ip=0,kpol
                  do i=0,kpol
                     res = coagintfluxfunction_GQ(kernel,K0,Q,vecnodes,vecweights,nbins,kpol,massgrid,mat_coeffs_leg,j,k,lp,l,ip,i)
                     if (res /= res) then
                        print*,"NAN in coagintfluxfunction_GQ for k=",kpol; stop
                     endif
                     tabintflux(j,k+1,lp,l,ip+1,i+1) = res
                  enddo
               enddo
            enddo
         enddo
      enddo
      call display_progress_bar(nbins,iprogress)
      iprogress = iprogress+1

   enddo

   print*,""

end subroutine compute_coagtabintflux_GQ


end module functions_tabflux_tabintflux


