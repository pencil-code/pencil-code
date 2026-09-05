!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


!-------------------------------------------
! MODULE: collision kernel functions
!-------------------------------------------

module collision_kernels
   use precision
   implicit none

contains

!> @brief Compute constant kernel 
!!
!! @param[in]    K0    constant value of the kernel function (used to adapt to code unit)
!! @param[in]    u     mass variable (colliding grain of mass u)
!! @param[in]    v     mass variable (colliding grain of mass v)
!! @return       res   evaluate constant kernel
real(wp) function kconst(K0,u,v) result(res)
   implicit none 
   real(wp), intent(in)  :: K0,u,v

   res = K0
   
end function

!> @brief Compute additive kernel 
!!
!! @param[in]    K0    constant value of the kernel function (used to adapt to code unit)
!! @param[in]    u     mass variable (colliding grain of mass u)
!! @param[in]    v     mass variable (colliding grain of mass v)
!! @return       res   evaluate additive kernel at u and v
real(wp) function kadd(K0,u,v) result(res)
   implicit none 
   real(wp), intent(in)  :: K0,u,v

   res = K0*(u+v)
   
end function

!> @brief Compute multiplicative kernel 
!!
!! @param[in]    K0    constant value of the kernel function (used to adapt to code unit)
!! @param[in]    u     mass variable (colliding grain of mass u)
!! @param[in]    v     mass variable (colliding grain of mass v)
!! @return       res   evaluate multiplicative kernel at u and v
real(wp) function kmul(K0,u,v) result(res)
   implicit none 
   real(wp), intent(in)  :: K0,u,v

   res = K0*u*v
   
end function

!> @brief Compute the cross-section term in the ballistic kernel K = cross-section * dv
!!
!! @param[in]    K0    constant value of the kernel function (used to adapt to code unit)
!! @param[in]    u     mass variable (colliding grain of mass u)
!! @param[in]    v     mass variable (colliding grain of mass v)
!! @return       res   evaluate cross-section term of the balistic kernel at u and v
real(wp) function k_cross_section(K0,u,v) result(res)
   implicit none 
   real(wp), intent(in)  :: K0,u,v

   res = K0 * (u**(1._wp/3._wp) + v**(1._wp/3._wp))**2

   
end function

!> @brief Compute collision kernel from Brownian motion, K = sigma * dv
!!
!! @param[in]    K0    constant value of the kernel function (used to adapt to code unit)
!! @param[in]    u     mass variable (colliding grain of mass u)
!! @param[in]    v     mass variable (colliding grain of mass v)
!! @return       res   Brownian motion collision kernel at u and v
real(wp) function k_Br(K0,u,v) result(res)
   implicit none 
   real(wp), intent(in)  :: K0,u,v

   res = K0 * (u**(1._wp/3._wp) + v**(1._wp/3._wp))**2 *sqrt(1._wp/u+1._wp/v)

end function


!> @brief Compute kernels at u and v
!!
!! @param[in]    kernel    select the collisional kernel function
!! @param[in]    K0        constant value of the kernel function (used to adapt to code unit)
!! @param[in]    u         mass variable (colliding grain of mass u)
!! @param[in]    v         mass variable (colliding grain of mass v)
!! @return       res       evaluate kernel at u and v
real(wp) function func_kernel(kernel,K0,u,v) result(res)
   implicit none 
   character(len=*),  intent(in) :: kernel
   real(wp),          intent(in) :: K0,u,v

   select case (kernel)
   case ("kconst")
      res = kconst(K0,u,v)
   case ("kadd")
      res = kadd(K0,u,v)
   case ("k_cross_section")
      res = k_cross_section(K0,u,v)
   case ("k_brownian")
      res = k_Br(K0,u,v)
   case default
      print*,kernel
      stop "Need to choose available collision kernel function in func_kernel"
   end select
   
end function


end module collision_kernels