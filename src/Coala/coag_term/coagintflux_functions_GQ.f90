!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************



!-------------------------------------------
! MODULE: coag intflux function GQ
!-------------------------------------------
module coagintflux_functions_GQ
   use precision
   use collision_kernels
   use polynomials_legendre
   implicit none

   public :: coagintfluxfunction_GQ

   private

contains


!> @brief Evaluate the integrand of the term including the integral of the coagulation flux, i.e. integrand of the triple integral
!!
!! @param[in]    kernel    select the collisional kernel function
!! @param[in]    K0        constant value of the kernel function (used to adapt to code unit)
!! @param[in]    k         order of polynomials
!! @param[in]    i         order of polynomials
!! @param[in]    ip        order of polynomials
!! @param[in]    ak        coefficients of polynomial of degree k, sorted from low to high order
!! @param[in]    ai        coefficients of polynomial of degree i, sorted from low to high order
!! @param[in]    aip       coefficients of polynomial of degree ip, sorted from low to high order
!! @param[in]    hj        size of mass bin j
!! @param[in]    u         mass variable (colliding grain of mass u)
!! @param[in]    v         mass variable (colliding grain of mass v)
!! @param[in]    xij       variable mapping the mass bin j in [-1,1], needed for Legendre polynomials
!! @param[in]    xilp      variable mapping the mass bin lp in [-1,1], needed for Legendre polynomials
!! @param[in]    xil       variable mapping the mass bin l in [-1,1], needed for Legendre polynomials
!! @return                 integrand of the term including the integral of the coagulation flux
real(wp) function func_coag_intflux(kernel,K0,k,i,ip,ak,ai,aip,hj,u,v,xij,xilp,xil) result(res)
   implicit none 
   character(len=*), intent(in) :: kernel
   integer,          intent(in) :: k,i,ip
   real(wp),         intent(in) :: ak(0:k),ai(0:i),aip(0:ip),u,v,xij,xilp,xil,hj,K0


   res = dphi_pol_k(k,ak,hj,xij)*func_kernel(kernel,K0,u,v)*phi_pol(ip,aip,xilp)*phi_pol(i,ai,xil)/v

   
end function


!> @brief Evaluate the triple integral for coagulation (term in DG scheme) depending only masses with Gauss-Legendre quadrature method.
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
!! @param[in]    k                 degree of polynomials in Legendre basis for approximation in bin j
!! @param[in]    lp                index corresponding to the mass of one colliding grain
!! @param[in]    l                 index corresponding to the mass of the second colliding grain
!! @param[in]    ip                degree of polynomials in Legendre basis for approximation in bin lp
!! @param[in]    i                 degree of polynomials in Legendre basis for approximation in bin l
!! @return                         triple integral for the coagulation flux evaluated at j,lp,l,ip,i
real(wp) function coagintfluxfunction_GQ(kernel,K0,Q,vecnodes,vecweights,nbins,kpol,massgrid,&
                                          mat_coeffs_leg,j,k,lp,l,ip,i) result(res)
   implicit none

   character(len=*), intent(in)  :: kernel
   integer,          intent(in)  :: Q,nbins,kpol
   real(wp),         intent(in)  :: vecnodes(Q),vecweights(Q),K0,mat_coeffs_leg(kpol+1,kpol+1)
   integer,          intent(in)  :: j,k,lp,l,ip,i
   real(wp),         intent(in)  :: massgrid(nbins+1)

   real(wp) :: xlgridl,xlgridr,xlpgridl,xlpgridr,xjgridl,xjgridr,xmin,xmax
   real(wp) :: xlp,hlp,xl,hl,xj,hj
   real(wp) :: xj_alpha,ulp_alpha,vl_alpha
   real(wp) :: xilp,xil,xij
   real(wp) :: a_vl,b_vl,a_ulp,b_ulp
   integer  :: alpha_u,alpha_v,alpha_x
   real(wp) :: aip(0:ip),ai(0:i),ak(0:k)

   xlgridr  = massgrid(l + 1)
   xlgridl  = massgrid(l)
   xlpgridl = massgrid(lp)
   xlpgridr = massgrid(lp + 1)

   xjgridr  = massgrid(j + 1)
   xjgridl  = massgrid(j)
   xmin     = massgrid(1)
   xmax     = massgrid(nbins+1)

   hlp = xlpgridr-xlpgridl
   xlp = 0.5_wp*(xlpgridr+xlpgridl)

   hl = xlgridr-xlgridl
   xl = 0.5_wp*(xlgridr+xlgridl)

   hj = xjgridr-xjgridl
   xj = 0.5_wp*(xjgridr+xjgridl)

   aip = mat_coeffs_leg(ip+1,:ip+1)
   ai  = mat_coeffs_leg(i+1,:i+1)
   ak  = mat_coeffs_leg(k+1,:k+1)

   res = 0._wp

   do alpha_x = 1,Q
      xj_alpha = xj + 0.5_wp*hj*vecnodes(alpha_x)
      xij = vecnodes(alpha_x)


      do alpha_u=1,Q
         a_ulp = max(xmin, xlpgridl)
         b_ulp = min(xj_alpha, xlpgridr)
         ulp_alpha = 0.5_wp*(b_ulp + a_ulp) + 0.5_wp*(b_ulp - a_ulp)*vecnodes(alpha_u)
         xilp = 2._wp*(ulp_alpha-xlp)/hlp

         do alpha_v=1,Q
            a_vl = max(xj_alpha - ulp_alpha + xmin, xlgridl)
            b_vl = min(xmax - ulp_alpha + xmin, xlgridr)

            vl_alpha = 0.5_wp*(b_vl + a_vl) + 0.5_wp*(b_vl - a_vl)*vecnodes(alpha_v)
            xil = 2._wp*(vl_alpha-xl)/hl

            if (xmax - ulp_alpha + xmin > xlgridl .and. xlgridr > xj_alpha - ulp_alpha + xmin) then

               res = res + 0.125_wp*hj*(b_ulp - a_ulp)*(b_vl - a_vl) &
                           *vecweights(alpha_x)*vecweights(alpha_u)*vecweights(alpha_v)&
                           *func_coag_intflux(kernel,K0,k,i,ip,ak,ai,aip,hj,ulp_alpha,vl_alpha,xij,xilp,xil)


            endif
                       
         enddo

      enddo
   enddo

   
end function


end module coagintflux_functions_GQ
