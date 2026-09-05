!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


!---------------------------------------------------------
! MODULE: Reconstruction polynomials of approximation of g
!---------------------------------------------------------
module reconstruction_g
   use precision
   use polynomials_legendre

contains



!> @brief Reconstruct the approximation of g using the components gij and the Legendre polynomial basis
!!
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    massbins          arithmetic mean value of massgrid for each mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[in]    gij               components of g on the polynomial basis
!! @param[in]    j                 indice of bin
!! @param[in]    x                 value to evalute the polynomial
!! @param[out]   res               reconstruction of g in bin j evaluated in x
subroutine recons_g(nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,j,x,res)
   implicit none
   integer,  intent(in)  :: kpol,nbins,j
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins),mat_coeffs_leg(kpol+1,kpol+1)
   real(wp), intent(in)  :: x
   real(wp), intent(in)  :: gij(nbins,kpol+1)
   real(wp), intent(out) :: res


   integer  :: i
   real(wp) :: xij,hj,xj
   real(wp) :: vecLegP(kpol+1)
   real(wp), allocatable :: ai(:)

   vecLegP = 0._wp
   hj = massgrid(j+1)-massgrid(j)
   xj = massbins(j)
   xij = 2._wp*(x-xj)/hj

   do i=0,kpol
      allocate(ai(i+1))
      ai = mat_coeffs_leg(i+1,:i+1)
      vecLegP(i+1)  = phi_pol(i,ai,xij)
      deallocate(ai)
   enddo

   res = dot_product(gij(j,:),vecLegP)


end subroutine recons_g

end module reconstruction_g
