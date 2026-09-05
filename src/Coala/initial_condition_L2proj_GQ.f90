!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


!--------------------------------------------------------------------
! MODULE: L2 projection on Legendre polynomials basis with GQ method
!         For initial condition of g
!--------------------------------------------------------------------
module initial_condition_L2proj_legbasis_GQ
   use functions_flux_intflux
   implicit none

contains

!> @brief Compute the initial gij coefficient from the dimensionless initial function g(x)=x*exp(-x)
!!
!! DG scheme k=0, piecewise constant approximation
!!
!! @param[in]    eps           minimum value for mass distribution approximation gij
!! @param[in]    nbins         number of dust bins
!! @param[in]    massgrid      grid of masses, borders value of mass bins
!! @param[in]    massbins      arithmetic mean value of massgrid for each mass bins
!! @param[in]    Q             number of points for Gauss-Legendre quadrature
!! @param[in]    vecnodes      nodes of the Legendre polynomials
!! @param[in]    vecweights    weights coefficients for the Gauss-Legendre polynomials
!! @param[out]   gij           initial components of g on the polynomial basis
subroutine L2proj_GQ_k0(eps,nbins,massgrid,massbins,Q,vecnodes,vecweights,gij)
   use precision
   use polynomials_legendre
   implicit none
   integer,  intent(in)  :: nbins,Q
   real(wp), intent(in)  :: eps
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins)
   real(wp), intent(in)  :: vecnodes(Q),vecweights(Q)
   real(wp), intent(out) :: gij(nbins)

   real(wp) :: c,xj,hj,term_sum,xjalpha,node,init_func
   integer  :: alpha,j

   do j=1,nbins
      
      xj = massbins(j)
      hj = massgrid(j+1)-massgrid(j)

      call coeffnorm(0,c)
      term_sum = 0._wp

      do alpha=1,Q
         node = vecnodes(alpha)
         xjalpha = xj + hj*node/2._wp

         if (xjalpha > 1e2_wp) then
            init_func = 0._wp
         else
            init_func = xjalpha*exp(-xjalpha)
         endif

         term_sum = term_sum + vecweights(alpha)*init_func
      enddo

      if ( (abs(term_sum) < eps .and. abs(term_sum) >= 0._wp)  ) then
         gij(j) = sign(eps,term_sum)
      else
         gij(j) = term_sum/c
      endif
   enddo

end subroutine L2proj_GQ_k0


!> @brief Compute the initial gij coefficient from the dimensionless initial function g(x)=x*exp(-x)
!!
!! DG scheme k>0, piecewise polynomial approximation
!!
!! @param[in]    eps               minimum value for mass distribution approximation gij
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    massbins          arithmetic mean value of massgrid for each mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[in]    Q                 number of points for Gauss-Legendre quadrature
!! @param[in]    vecnodes          nodes of the Legendre polynomials
!! @param[in]    vecweights        weights coefficients for the Gauss-Legendre polynomials
!! @param[out]   gij               initial components of g on the polynomial basis
subroutine L2proj_GQ(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,Q,vecnodes,vecweights,gij)
   use precision
   use polynomials_legendre
   implicit none
   integer,  intent(in)  :: nbins,kpol,Q
   real(wp), intent(in)  :: eps
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins),mat_coeffs_leg(kpol+1,kpol+1)
   real(wp), intent(in)  :: vecnodes(Q),vecweights(Q)
   real(wp), intent(out) :: gij(nbins,kpol+1)

   real(wp) :: c,xj,hj,term_sum,xjalpha,node,init_func
   integer  :: alpha,k,j
   real(wp), allocatable :: ak(:)

   do j=1,nbins
      xj = massbins(j)
      hj = massgrid(j+1)-massgrid(j)

      do k=0,kpol
         allocate(ak(k+1))
         ak = mat_coeffs_leg(k+1,:k+1)

         call coeffnorm(k,c)
         term_sum = 0._wp

         do alpha=1,Q
            node = vecnodes(alpha)
            xjalpha = xj + hj*node/2._wp
            
            if (xjalpha > 1e2_wp) then
               init_func = 0._wp
            else
               init_func = xjalpha*exp(-xjalpha)
            endif

            term_sum = term_sum + vecweights(alpha)*init_func*phi_pol(k,ak,vecnodes(alpha))
         
            ! print*,"phi_pol(k,ak,vecnodes(alpha))=",phi_pol(k,ak,vecnodes(alpha))
         enddo

         if ( (abs(term_sum) < eps .and. abs(term_sum) >= 0._wp)  ) then
            gij(j,k+1) = sign(eps,term_sum)
         else
            gij(j,k+1) = term_sum/c
         endif

         deallocate(ak)

      enddo

      !check negative values for mass density distribution
      if (gij(j,1) < 0._wp) then
         print*,"j=",j,", gij =",gij(j,1)
         stop

      else if ( gij(j,1) <= eps  ) then
         gij(j,1) = eps
         gij(j,2:) = 0._wp
      endif
   enddo

end subroutine L2proj_GQ


!> @brief Compute the initial gij coefficient from interpolation of cumulative function of the dust density
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
!! @param[in]    eps_rhodust       minimum value for dust density
!! @param[in]    rhodust           1D array dust density for each grain size
!! @param[out]   gij               initial components of g on the polynomial basis
subroutine L2proj_gij_GQ(nbins,kpol,massgrid,massbins,mat_coeffs_leg,Q,vecnodes,vecweights,eps_rhodust,rhodust,gij)
   use precision
   use polynomials_legendre
   use tension_mod
   implicit none
   integer,  intent(in)  :: nbins,kpol,Q
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins),mat_coeffs_leg(kpol+1,kpol+1)
   real(wp), intent(in)  :: vecnodes(Q),vecweights(Q)
   real(wp), intent(in)  :: eps_rhodust
   real(wp), intent(in)  :: rhodust(nbins)
   real(wp), intent(out) :: gij(nbins,kpol+1)

   real(wp) :: c,xj,hj,term_sum,xjalpha,node,init_func
   real(wp) :: eps_gij
   integer  :: alpha,k,j
   integer  :: IER,SigErr
   real(wp) :: Vfunc(nbins+1),YP(nbins+1),SIGMA(nbins+1)
   real(wp), allocatable :: ak(:)


   !value of the cumulative function Vfunc
   Vfunc = 0._wp
   do j=2,nbins+1
      Vfunc(j) = sum(rhodust(1:j-1))
   enddo

   ! print*,"Vfunc=",Vfunc
   ! stop

   YP = 0._wp
   SIGMA = 0._wp

   !interpolation of the cumulative function Vfunc to get gij
   
   !interpolation Vfunc(x)
   ! call TSPSI (nbins+1,massgrid, Vfunc,YP, SIGMA,IER,SigErr)

   !interpolation Vfunc(log(x))
   call TSPSI (nbins+1,log(massgrid), Vfunc,YP, SIGMA,IER,SigErr)

   ! print*,"YP=",YP

   gij = 0._wp
   eps_gij = eps_rhodust/massgrid(nbins+1)

   do j=1,nbins
      if (rhodust(j) > eps_rhodust) then
         xj = massbins(j)
         hj = massgrid(j+1)-massgrid(j)

         do k=0,kpol
            allocate(ak(k+1))
            ak = mat_coeffs_leg(k+1,:k+1)

            call coeffnorm(k,c)
            term_sum = 0._wp

            do alpha=1,Q
               node = vecnodes(alpha)
               xjalpha = xj + hj*node/2._wp

               !interpolation V(x)
               ! init_func = HPVAL(xjalpha,nbins+1,massgrid, Vfunc,YP,SIGMA, IER)
               
               !interpolation V(log(x))
               init_func = HPVAL(log(xjalpha),nbins+1,log(massgrid), Vfunc,YP,SIGMA, IER)/xjalpha
               

               term_sum = term_sum + vecweights(alpha)*init_func*phi_pol(k,ak,vecnodes(alpha))
            enddo

            gij(j,k+1) = term_sum/c

            deallocate(ak)

         enddo

         !enforce mass conservation after interpolation + L2 projection
         gij(j,1) = rhodust(j)/(massgrid(j+1)-massgrid(j))

         !check negative values for mass density distribution
         if (gij(j,1) < 0._wp) then
            print*,"j=",j,", gij =",gij(j,1)
            stop

         else if ( gij(j,1) <= eps_gij  ) then
            gij(j,1) = eps_gij
            gij(j,2:) = 0._wp
         endif

         

      endif
   enddo


end subroutine L2proj_gij_GQ




!> @brief Compute the initial gij coefficient from interpolation of cumulative function of the dust density, with indices optimisation (only with gij > eps)
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
!! @param[in]    eps_rhodust       minimum value for dust density
!! @param[in]    rhodust           1D array dust density for each grain size
!! @param[out]   gij               initial components of g on the polynomial basis
subroutine L2proj_gij_GQ_opt(nbins,kpol,massgrid,massbins,mat_coeffs_leg,Q,vecnodes,vecweights,eps_rhodust,rhodust,gij)
   use precision
   use polynomials_legendre
   use tension_mod
   implicit none

   integer,  intent(in)  :: nbins,kpol,Q
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins),mat_coeffs_leg(kpol+1,kpol+1)
   real(wp), intent(in)  :: vecnodes(Q),vecweights(Q)
   real(wp), intent(in)  :: eps_rhodust
   real(wp), intent(in)  :: rhodust(nbins)
   real(wp), intent(out) :: gij(nbins,kpol+1)
   

   real(wp) :: c,xj,hj,term_sum,xjalpha,node,init_func,M1
   real(wp) :: eps_gij
   integer  :: alpha,k,j
   integer  :: IER,SigErr
   integer  :: ind_min,ind_max,nbins_loc

   real(wp), allocatable :: Vfunc(:),YP(:),SIGMA(:)
   real(wp), allocatable :: ak(:)


   !range rhodust > eps_rhodust
   call compute_min_max_indices(eps_rhodust,nbins,massbins,rhodust,ind_min,ind_max)
   ! print*,"ind_min=",ind_min
   ! print*,"ind_max=",ind_max
   ! stop

   nbins_loc = ind_max-ind_min+1
   ! print*,"nbins_loc=",nbins_loc

   ! cumulative function Vfunc
   allocate(Vfunc(nbins_loc+1))
   Vfunc = 0._wp
   do j=ind_min+1,ind_max+1
      ! print*,"j=",j
      Vfunc(j-ind_min+1) = sum(rhodust(ind_min:j-1))
   enddo

   ! print*,"Vfunc=",Vfunc
   ! print*,"rhodust(ind_min:nbins_loc)=",rhodust(ind_min:ind_max)
   ! stop

   allocate(YP(nbins_loc+1))
   allocate(SIGMA(nbins_loc+1))
   YP = 0._wp
   SIGMA = 0._wp

   !interpolation V(x)
   ! call TSPSI (nbins_loc+1,massgrid(ind_min:ind_max+1), Vfunc,YP, SIGMA,IER,SigErr)

   !interpolation V(log(x))
   call TSPSI (nbins_loc+1,log(massgrid(ind_min:ind_max+1)), Vfunc,YP, SIGMA,IER,SigErr)

   gij = 0._wp
   eps_gij = eps_rhodust/massgrid(nbins+1)

   !interpolation to get gij
   do j=ind_min,ind_max
      xj = massbins(j)
      hj = massgrid(j+1)-massgrid(j)

      do k=0,kpol
         allocate(ak(k+1))
         ak = mat_coeffs_leg(k+1,:k+1)

         call coeffnorm(k,c)
         term_sum = 0._wp

         do alpha=1,Q
            node = vecnodes(alpha)
            xjalpha = xj + hj*node/2._wp

            !interpolation V(x)
            ! init_func = HPVAL(xjalpha,nbins_loc+1,massgrid(ind_min:ind_max+1), Vfunc,YP,SIGMA, IER)
            
            !interpolation V(log(x))
            init_func = HPVAL(log(xjalpha),nbins_loc+1,log(massgrid(ind_min:ind_max+1)), Vfunc,YP,SIGMA, IER)/xjalpha
            

            term_sum = term_sum + vecweights(alpha)*init_func*phi_pol(k,ak,vecnodes(alpha))
         enddo

         gij(j,k+1) = term_sum/c

         deallocate(ak)

      enddo

      !enforce mass conservation after interpolation + L2 projection
      if ( gij(j,1) > eps_gij) then
         gij(j,1) = rhodust(j)/(massgrid(j+1)-massgrid(j))
      endif

      !check negative values for mass density distribution
      if (gij(j,1) < 0._wp) then
         print*,"j=",j,", gij =",gij(j,1)
         stop

      else if ( gij(j,1) <= eps_gij  ) then
         gij(j,1) = eps_gij
         gij(j,2:) = 0._wp
      endif
   enddo

   ! do j=1,nbins
   !    print*,"j=",j,", gij =",reshape(gij(j,:),(/kpol+1/))
   ! enddo

   ! print*,"total mass init =",sum(rhodust)

   M1 = 0._wp
   do j=1,nbins
      M1 = M1 + gij(j,1)*(massgrid(j+1)-massgrid(j))
   enddo
   ! print*,"total mass after interp =",M1
   ! print*,"abs error =", abs(sum(rhodust) - M1)/sum(rhodust)

   
   deallocate(Vfunc)
   deallocate(YP)
   deallocate(SIGMA)


end subroutine L2proj_gij_GQ_opt


end module initial_condition_L2proj_legbasis_GQ
