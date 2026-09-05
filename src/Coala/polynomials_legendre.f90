!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


!-------------------------------------------
! MODULE: functions for Legendre polynomials
!-------------------------------------------

module polynomials_legendre
   use precision
   implicit none

contains

!> @brief Normalisation coefficient of Legendre polynomials
!!
!! @param[in]    i   order of polynomials
!! @param[out]   c   normalisation coefficient
subroutine coeffnorm(i,c)
   implicit none
   integer, intent(in)  :: i
   real(wp),intent(out) :: c

   c = 2._wp/(2._wp*i+1._wp)
end subroutine coeffnorm


!> @brief Coefficients of Legendre polynomial of order k 
!!
!! @param[in]    k        order of polynomials
!! @param[out]       coeffs   polynomials coefficients, sorted from low to high order 
subroutine coeff_Leg(k,coeffs)
   implicit none
   integer,  intent(in)  :: k
   real(wp), intent(out) :: coeffs(k+1)

   select case (k)
      case (0)
         coeffs = [1._wp]
      case (1)
         coeffs = [0._wp,1._wp]
      case (2)
         coeffs = [-0.5_wp,0._wp,1.5_wp]
      case (3)
         coeffs = [0._wp,-1.5_wp,0._wp,2.5_wp]
      case (4)
         coeffs = [3._wp/8._wp,0._wp,-30._wp/8._wp,0._wp,35._wp/8._wp]
      case (5)
         coeffs = [0._wp,15._wp/8._wp,0._wp,-70._wp/8._wp,0._wp,63._wp/8._wp]
      case (6)
         coeffs = [-5._wp/16._wp,0._wp,105._wp/16._wp,0._wp,-315._wp/16._wp,0._wp,231._wp/16._wp]
      case (7)
         coeffs = [0._wp,-35._wp/16._wp,0._wp,315._wp/16._wp,0._wp,-693._wp/16._wp,0._wp,429._wp/16._wp]
      case (8)
         coeffs = [35._wp/128._wp,0._wp,-1260._wp/128._wp,0._wp,6930._wp/128._wp,0._wp,-12012._wp/128._wp,&
                   0._wp,6435._wp/128._wp]
      case (9)
         coeffs = [0._wp,315._wp/128._wp,0._wp,-4620._wp/128._wp,0._wp,18018._wp/128._wp,0._wp,&
                   -25740._wp/128._wp,0._wp,12155._wp/128._wp]
      case(10)
         coeffs = [-63._wp/256._wp,0._wp,3465._wp/256._wp,0._wp,-30030._wp/256._wp,0._wp,&
                   90090._wp/256._wp,0._wp,-109395._wp/256._wp,0._wp,46189._wp/256._wp]
      case default
         stop 'polynomials_legendre.f90 -> subr coeff_Leg, Wrong order, need kpol <= 10'
   end select


end subroutine coeff_Leg


!> @brief Matrix of coefficients of Legendre polynomial of order up to k 
!!
!! @param[in]    k            order of polynomials
!! @param[out]   mat_coeffs   polynomials coefficients, sorted from low to high order 
subroutine compute_mat_coeffs(k,mat_coeffs)
   implicit none
   integer, intent(in)  :: k
   real(wp),intent(out) :: mat_coeffs(k+1,k+1)

   integer               :: i
   real(wp), allocatable :: coeffs(:)

   mat_coeffs = 0._wp

   do i=0,k
      allocate(coeffs(i+1))
      call coeff_Leg(i,coeffs)
      mat_coeffs(i+1,:i+1) = coeffs

      deallocate(coeffs)
   enddo

end subroutine compute_mat_coeffs


!> @brief Evaluate polynomial sum_{i=0}^{k} a_i x^i by Horner's method
!!
!! @param[in]    i    order of polynomials
!! @param[in]    ai   polynomials coefficients, sorted from low to high order
!! @return            evaluation of polynomial of order i at x
real(wp) function phi_pol(i,ai,x) result(res)
   implicit none 
   integer,  intent(in) :: i
   real(wp), intent(in) :: ai(0:i),x

   integer  :: j


   res = 0._wp
   ! Horner method need coefficients from high to low ordre => reverse do loop
   do j=i,0,-1
      res = res * x + ai(j)
   enddo
   
end function


!> @brief Coefficients of the derivative of Legendre polynomial of order k
!!
!! @param[in]    k             order of polynomials
!! @param[in]    pol_coeffs    polynomials coefficients, sorted from low to high order
!! @param[out]   dpol_coeffs   coefficients of the derivative of polynomials, sorted from low to high order
subroutine polynomial_derivative_coeffs(k, pol_coeffs,dpol_coeffs)
   implicit none
   integer, intent(in)  :: k
   real(wp),intent(in)  :: pol_coeffs(0:k)
   real(wp),intent(out) :: dpol_coeffs(k)

   integer :: i

   if (k==0) then
      dpol_coeffs = [0._wp]
   else
      do i=1,k
         dpol_coeffs(i) = i * pol_coeffs(i)
      enddo
   endif
end subroutine polynomial_derivative_coeffs


!> @brief Derivative of P_k(xij) with respect to x, where xij = 2/hj*(x-xj)
!!
!! @param[in]    k    order of polynomials
!! @param[in]    ak   polynomials coefficients, sorted from low to high order
!! @param[in]    hj   width of bin j
!! @param[in]    xij  variable mapping the mass bin j in [-1,1], needed for Legendre polynomials
!! @return            evaluation at xij of the derivative of P_k(xij) with respect to x
real(wp) function dphi_pol_k(k,ak,hj,xij) result(res)
   implicit none 
   integer,  intent(in) :: k
   real(wp), intent(in) :: ak(0:k),hj,xij

   real(wp) :: dpol_coeffs(k)

   if (k==0) then
      res = 0._wp
   else
      call polynomial_derivative_coeffs(k,ak,dpol_coeffs)
      res = phi_pol(k-1,dpol_coeffs,xij)*2._wp/hj
   endif

   
end function



end module polynomials_legendre