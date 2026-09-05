!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************

!-------------------------------------------
! MODULE: coag flux function GQ
!-------------------------------------------
module coagflux_functions_GQ
   use precision
   use collision_kernels
   use polynomials_legendre
   implicit none

   public :: coagfluxfunction_GQ

   private

contains

!> @brief Evaluate the integrand of the coagulation flux, i.e. integrand of the double integral
!!
!! @param[in]    kernel    select the collisional kernel function
!! @param[in]    K0        constant value of the kernel function (used to adapt to code unit)
!! @param[in]    i         order of polynomials
!! @param[in]    ip        order of polynomials
!! @param[in]    ai        coefficients of polynomial of degree i, sorted from low to high order
!! @param[in]    aip       coefficients of polynomial of degree ip, sorted from low to high order
!! @param[in]    u         mass variable (colliding grain of mass u)
!! @param[in]    v         mass variable (colliding grain of mass v)
!! @param[in]    xilp      variable mapping the mass bin lp in [-1,1], needed for Legendre polynomials
!! @param[in]    xil       variable mapping the mass bin l in [-1,1], needed for Legendre polynomials
!! @return                 integrand of cogulation flux evaluated at u,v,xilp,xil
real(wp) function func_coag_flux(kernel,K0,i,ip,ai,aip,u,v,xilp,xil) result(res)
   implicit none 
   character(len=*), intent(in) :: kernel
   integer,          intent(in) :: i,ip
   real(wp),         intent(in) :: ai(0:i),aip(0:ip),u,v,xilp,xil,K0

   res = func_kernel(kernel,K0,u,v)*phi_pol(ip,aip,xilp)*phi_pol(i,ai,xil)/v

end function


!> @brief Evaluate the double integral for coagulation flux depending only masses with Gauss-Legendre quadrature method.
!!
!! This function is used to calculate the array for the coagulation flux as precomputation.
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
!! @param[in]    j                 index corresponding to the mass of the new formed grain
!! @param[in]    lp                index corresponding to the mass of one colliding grain
!! @param[in]    l                 index corresponding to the mass of the second colliding grain
!! @param[in]    ip                degree of polynomials in Legendre basis for approximation in bin lp
!! @param[in]    i                 degree of polynomials in Legendre basis for approximation in bin l
!! @return                         double integral for the coagulation flux evaluated at j,lp,l,ip,i
real(wp) function coagfluxfunction_GQ(kernel,K0,Q,vecnodes,vecweights,nbins,kpol,massgrid,&
                                       mat_coeffs_leg,j,lp,l,ip,i) result(res)
   implicit none

   character(len=*), intent(in)  :: kernel
   integer,          intent(in)  :: Q,nbins,kpol
   real(wp),         intent(in)  :: K0,vecnodes(Q),vecweights(Q),mat_coeffs_leg(kpol+1,kpol+1)
   integer,          intent(in)  :: j,lp,l,ip,i
   real(wp),         intent(in)  :: massgrid(nbins+1)

   real(wp) :: xlgridl,xlgridr,xlpgridl,xlpgridr,xjgridr,xmin,xmax
   real(wp) :: xlp,hlp,xl,hl
   real(wp) :: ulp_alpha,vl_alpha,xilp,xil,a_vl,b_vl
   integer  :: alpha_u,alpha_v
   real(wp) :: aip(0:ip),ai(0:i)


   xlgridr  = massgrid(l + 1)
   xlgridl  = massgrid(l)
   xlpgridl = massgrid(lp)
   xlpgridr = massgrid(lp + 1)

   hlp = xlpgridr-xlpgridl
   xlp = 0.5_wp*(xlpgridr+xlpgridl)

   hl = xlgridr-xlgridl
   xl = 0.5_wp*(xlgridr+xlgridl)

   xjgridr  = massgrid(j+1)
   xmin     = massgrid(1)
   xmax     = massgrid(nbins+1)


   aip = mat_coeffs_leg(ip+1,:ip+1)
   ai  = mat_coeffs_leg(i+1,:i+1)

   res = 0._wp

   do alpha_u=1,Q
      ulp_alpha = xlp + 0.5_wp*hlp*vecnodes(alpha_u)
      xilp = vecnodes(alpha_u)

      do alpha_v=1,Q
         a_vl = max(xjgridr - ulp_alpha + xmin, xlgridl)
         b_vl = min(xmax - ulp_alpha + xmin, xlgridr)

         vl_alpha = 0.5_wp*(b_vl + a_vl) + 0.5_wp*(b_vl - a_vl)*vecnodes(alpha_v)
         xil = 2._wp*(vl_alpha-xl)/hl

         if (xmax - ulp_alpha + xmin > xlgridl .and. xlgridr > xjgridr - ulp_alpha + xmin) then

            res = res + 0.25_wp*hlp*(b_vl - a_vl) &
                        *vecweights(alpha_u)*vecweights(alpha_v)&
                        *func_coag_flux(kernel,K0,i,ip,ai,aip,ulp_alpha,vl_alpha,xilp,xil)


         endif
                    
      enddo

   enddo

end function

end module coagflux_functions_GQ
