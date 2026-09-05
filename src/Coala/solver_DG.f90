!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


!--------------------
! MODULE: solver DG
!--------------------

module solver_DG
   use precision
   use functions_flux_intflux
   use polynomials_legendre
   use limiter
   implicit none

contains

!-----------------------------------------------------
! Original version to compute functions for DG scheme
!-----------------------------------------------------

!> @brief Coagulation CFL for DG scheme k=0 piecewise constant approximation (simple kernels)
!!
!! CFL formulation from Filbet & Laurencot 2004, dt <= mean_g * dm/dF
!!
!! @param[in]    eps            minimum value for mass distribution approximation gij
!! @param[in]    nbins          number of dust bins
!! @param[in]    massgrid       grid of masses, borders value of mass bins
!! @param[in]    gij            components of g on the polynomial basis
!! @param[in]    tabflux_coag   3D array to evaluate coagulation flux
!! @return                      CFL restriction for coagulation, DG k=0
real(wp) function compute_CFL_k0(eps,nbins,massgrid,gij,tabflux_coag) result(CFL)
   implicit none
   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: massgrid(nbins+1),eps
   real(wp), intent(in)  :: gij(nbins)
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins)
   
   real(wp) :: flux(nbins),tabdtCFL(nbins)
   real(wp) :: arr_gij(nbins,nbins) 
   real(wp) :: hj
   integer  :: j

   call compute_arr_gij_k0(nbins,gij,arr_gij)
   call compute_flux_k0(nbins,arr_gij,tabflux_coag,flux)

   ! for debug
   ! print*,"display flux from compute_CFL_k0"
   ! do j=1,nbins
   !    print*,"j=",j,", gij = ",gij(j),",flux =",flux(j)
   ! enddo

   tabdtCFL = 0._wp

   !loop in mass domain [massmin,massmax]
   do j=1,nbins
      hj = massgrid(j+1)-massgrid(j)
      
      if (gij(j) > eps) then
         if (j==1) then
            !bin 1, where flux(0) = 0, no mass entering the mass domain [massmin,massmax]
            tabdtCFL(j) = abs(gij(j)*hj/flux(j))
         else
            tabdtCFL(j) = abs(gij(j)*hj/(flux(j)-flux(j-1)))
         endif
      endif

   enddo
      
   CFL = minval(tabdtCFL,tabdtCFL > 0._wp)

   ! for debug
   ! do j=1,nbins
   !    print*,"j=",j,", tabCFL(j) = ",tabdtCFL(j)
   ! enddo
   ! print*,"CFL = ",CFL
   ! stop

   ! check CFL value
   if ( CFL == 0._wp) then
      print*,"j,gij,hj,flux,tabdtCFL"
      do j=1,nbins
         print*,j,gij(j),massgrid(j+1)-massgrid(j),flux(j),tabdtCFL(j)
      enddo
      print*,"CFL = ",CFL
      print*,"issue in CFL coagulation k=0"
      stop
   endif


end function




!> @brief Coagulation CFL for DG scheme k=0 piecewise constant approximation (physical kernel)
!!
!! CFL formulation from Filbet & Laurencot 2004, dt <= mean_g * dm/dF
!!
!! @param[in]    eps            minimum value for mass distribution approximation gij
!! @param[in]    nbins          number of dust bins
!! @param[in]    massgrid       grid of masses, borders value of mass bins
!! @param[in]    gij            components of g on the polynomial basis
!! @param[in]    tabflux_coag   3D array to evaluate coagulation flux
!! @param[in]    dv             2D array of the differential velocity between grains in bins
!! @return                      CFL restriction for coagulation, DG k=0
real(wp) function compute_CFL_k0_kdv(eps,nbins,massgrid,gij,tabflux_coag,dv) result(CFL)
   implicit none
   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: massgrid(nbins+1),eps
   real(wp), intent(in)  :: gij(nbins)
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins),dv(nbins,nbins)
   
   real(wp) :: flux(nbins),tabdtCFL(nbins)
   real(wp) :: arr_gij_dv(nbins,nbins) 
   real(wp) :: hj
   integer  :: j


   call compute_arr_gij_k0_kdv(nbins,gij,dv,arr_gij_dv)
   call compute_flux_k0(nbins,arr_gij_dv,tabflux_coag,flux)

   ! for debug
   ! do j=1,nbins
   !    print*,"j=",j,", gij = ",gij(j),",flux =",flux(j)
   ! enddo

   tabdtCFL = 0._wp

   !loop in mass domain [massmin,massmax]
   do j=1,nbins
      hj = massgrid(j+1)-massgrid(j)
      
      if (gij(j) > eps) then
         if (j==1) then
            !bin 1, where flux(0) = 0, no mass entering the mass domain [massmin,massmax]
            tabdtCFL(j) = abs(gij(j)*hj/flux(j))
         else
            tabdtCFL(j) = abs(gij(j)*hj/(flux(j)-flux(j-1)))
         endif
      endif

   enddo
      
   CFL = minval(tabdtCFL,tabdtCFL > 0._wp)

   ! for debug
   ! do j=1,nbins
   !    print*,"j=",j,", tabCFL(j) = ",tabdtCFL(j)
   ! enddo
   ! print*,"CFL = ",CFL
   ! stop

   if ( CFL ==0._wp) then
      print*,"j,gij,hj,flux,tabdtCFL"
      do j=1,nbins
         print*,j,gij(j),massgrid(j+1)-massgrid(j),flux(j),tabdtCFL(j)
      enddo
      print*,"CFL = ",CFL
      print*,"issue in CFL coagulation k=0, physical kernel cross-section * dv"
      stop
   endif

end function






!> @brief compute the DG operator L for piecewise constant approximation and simple kernels (see Lombart et al., 2021) 
!!
!! It is used for the time solver
!!
!! @param[in]    nbins          number of dust bins
!! @param[in]    massgrid       grid of masses, borders value of mass bins
!! @param[in]    gij            components of g on the polynomial basis
!! @param[in]    tabflux_coag   3D array to evaluate coagulation flux
!! @param[out]   L_k0            DG operator for piecewise constant approximation in each bin
subroutine operator_L_k0(nbins,massgrid,gij,tabflux_coag,L_k0)
   implicit none
   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: massgrid(nbins+1)
   real(wp), intent(in)  :: gij(nbins)
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins)
   real(wp), intent(out) :: L_k0(nbins)

   real(wp) :: flux(nbins)
   real(wp) :: arr_gij(nbins,nbins)
   real(wp) :: hj
   integer  :: j


   call compute_arr_gij_k0(nbins,gij,arr_gij)
   call compute_flux_k0(nbins,arr_gij,tabflux_coag,flux)

   ! for debug
   ! do j=1,nbins
   !    print*,"j=",j,",flux =",flux(j)
   ! enddo

   L_k0 = 0._wp
   !bin 1 with flux(0) = 0
   L_k0(1) = -flux(1)/(massgrid(2)-massgrid(1))
   do j=2,nbins
      hj = massgrid(j+1)-massgrid(j)

      ! L_k0(j) = -(flux(j) - flux(j-1))/hj

      !normalisation to reduce errors from substraction operation
      L_k0(j) = -(flux(j)/hj - flux(j-1)/hj)

   enddo

end subroutine operator_L_k0





!> @brief compute the DG operator L for piecewise constant approximation and physical kernel 
!!
!! It is used for the time solver
!!
!! @param[in]    nbins          number of dust bins
!! @param[in]    massgrid       grid of masses, borders value of mass bins
!! @param[in]    gij            components of g on the polynomial basis
!! @param[in]    tabflux_coag   3D array to evaluate coagulation flux
!! @param[in]    dv             2D array of the differential velocity between grains in bins
!! @param[out]   L_k0            DG operator for piecewise constant approximation in each bin
subroutine operator_L_k0_kdv(nbins,massgrid,gij,tabflux_coag,dv,L_k0)
   implicit none
   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: massgrid(nbins+1)
   real(wp), intent(in)  :: gij(nbins)
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins),dv(nbins,nbins)
   real(wp), intent(out) :: L_k0(nbins)

   real(wp) :: flux(nbins)
   real(wp) :: arr_gij_dv(nbins,nbins)
   real(wp) :: hj
   integer  :: j

   
   call compute_arr_gij_k0_kdv(nbins,gij,dv,arr_gij_dv)
   call compute_flux_k0(nbins,arr_gij_dv,tabflux_coag,flux)

   ! for debug
   ! do j=1,nbins
   !    print*,"j=",j,",flux =",flux(j)
   ! enddo
   ! stop
   
   L_k0 = 0._wp

   !bin 1 with flux(0) = 0
   L_k0(1) = -flux(1)/(massgrid(2)-massgrid(1))

   do j=2,nbins
      hj = massgrid(j+1)-massgrid(j)

      !original expression
      ! L_k0(j) = -(flux(j) - flux(j-1))/hj

      !normalisation to reduce errors from substraction
      L_k0(j) = -(flux(j)/hj - flux(j-1)/hj)
   enddo


end subroutine operator_L_k0_kdv






!> @brief Function to compute SSPRK order 3 time solver with piecewise constant approximation for simple kernels
!!
!! See Zhang & Shu 2010 and Lombart et al., 2021
!!
!! @param[in]     eps            minimum value for mass distribution approximation gij
!! @param[in]     nbins          number of dust bins
!! @param[in]     massgrid       grid of masses, borders value of mass bins
!! @param[in]     gij            components of g on the polynomial basis
!! @param[in]     tabflux_coag   3D array to evaluate coagulation flux
!! @param[in]     dt             timestep
!! @param[out]    gijnew         evolved components of g on the polynomial basis
subroutine time_solver_k0(eps,nbins,massgrid,gij,tabflux_coag,dt,gijnew)

   implicit none

   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: massgrid(nbins+1)
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins)
   real(wp), intent(in)  :: dt,eps
   real(wp), intent(in)  :: gij(nbins)
   real(wp), intent(out) :: gijnew(nbins)


   real(wp) :: gij_1(nbins),gij_2(nbins)
   real(wp) :: L_k0(nbins),L_k0_1(nbins),L_k0_2(nbins)

   integer :: j,i


   !SSPRK3 algorithm (Zhang & Shu 2010)
   !step 1
   call operator_L_k0(nbins,massgrid,gij,tabflux_coag,L_k0)

   gij_1 = gij + dt*L_k0


   !check gij values and limit to eps 
   do j=1,nbins
      if ( gij_1(j) < 0._wp  ) then
         print*,"error in calculating gij_1"
         print*,"i,gij,dt*Lk0"
         do i=1,nbins
            print*,i,gij(i),dt*L_k0(i)
         enddo
         print*,"j=",j,", gij_1 =",gij_1(j)
         stop
      else if ( gij_1(j) <= eps  ) then
         gij_1(j) = eps

      endif
   enddo


   !step 2
   call operator_L_k0(nbins,massgrid,gij_1,tabflux_coag,L_k0_1)

   gij_2 = 3._wp*gij/4._wp + (gij_1 + dt*L_k0_1)/4._wp


   !check gij values and limit to eps 
   do j=1,nbins
      if ( gij_2(j) < 0._wp  ) then
         print*,"error in calculating gij_2"
         print*,"i,gij1,dt*Lk0_2"
         do i=1,nbins
            print*,i,gij_1(i),dt*L_k0_1(i)
         enddo
         print*,"j=",j,", gij_2 =",gij_2(j)
         stop

      else if ( gij_2(j) <= eps  ) then
         gij_2(j) = eps
      endif
   enddo


   !step 3
   call operator_L_k0(nbins,massgrid,gij_2,tabflux_coag,L_k0_2)

   gijnew = gij/3._wp + 2._wp*(gij_2 + dt*L_k0_2)/3._wp


   !check gij values and limit to eps 
   do j=1,nbins
      if ( gijnew(j) < 0._wp ) then
         print*,"error in calculating gijnew"
         print*,"i,gij_2,dt*Lk0_2"
         do i=1,nbins
            print*,i,gij_2(i),dt*L_k0_2(i)
         enddo
         print*,"j=",j,", gijnew =",gijnew(j)
         stop
      else if ( gijnew(j) <= eps  ) then
         gijnew(j) = eps
      endif
   enddo

end subroutine time_solver_k0




!> @brief Function to compute SSPRK order 3 time solver with piecewise constant approximation for physical kernel
!!
!! See Zhang & Shu 2010 and Lombart et al., 2021
!!
!! @param[in]     eps            minimum value for mass distribution approximation gij
!! @param[in]     nbins          number of dust bins
!! @param[in]     massgrid       grid of masses, borders value of mass bins
!! @param[in]     gij            components of g on the polynomial basis
!! @param[in]     tabflux_coag   3D array to evaluate coagulation flux
!! @param[in]     dv             2D array of the differential velocity between grains in bins lp and l
!! @param[in]     dt             timestep
!! @param[out]    gijnew         evolved components of g on the polynomial basis
subroutine time_solver_k0_kdv(eps,nbins,massgrid,gij,tabflux_coag,dv,dt,gijnew)
   implicit none

   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: massgrid(nbins+1)
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins),dv(nbins,nbins)
   real(wp), intent(in)  :: dt,eps
   real(wp), intent(in)  :: gij(nbins)
   real(wp), intent(out) :: gijnew(nbins)


   real(wp) :: gij_1(nbins),gij_2(nbins)
   real(wp) :: L_k0(nbins),L_k0_1(nbins),L_k0_2(nbins)

   integer :: j,i


   !SSPRK3 algo (Zhang & Shu 2010)
   !step 1
   call operator_L_k0_kdv(nbins,massgrid,gij,tabflux_coag,dv,L_k0)

   gij_1 = gij + dt*L_k0

   !check gij values and limit to eps
   do j=1,nbins
      if ( gij_1(j) < 0._wp  ) then
         print*,"error in calculating gij_1"
         print*,"i,gij,dt*L_k0"
         do i=1,nbins
            print*,i,gij(i),dt*L_k0(i)
         enddo
         print*,"j=",j,", gij_1 =",gij_1(j)
         stop
      else if ( gij_1(j) <= eps  ) then
         gij_1(j) = eps
      endif
   enddo


   !step 2
   call operator_L_k0_kdv(nbins,massgrid,gij_1,tabflux_coag,dv,L_k0_1)

   gij_2 = 3._wp*gij/4._wp + (gij_1 + dt*L_k0_1)/4._wp

   !check gij values and limit to eps
   do j=1,nbins
      if ( gij_2(j) < 0._wp  ) then
         print*,"error in calculating gij_2"
         print*,"i,gij_1,dt*L_k0_2"
         do i=1,nbins
            print*,i,gij_1(i),dt*L_k0_1(i)
         enddo
         print*,"j=",j,", gij_2 =",gij_2(j)
         stop

      else if ( gij_2(j) <= eps  ) then
         gij_2(j) = eps
      endif
   enddo


   !step 3
   call operator_L_k0_kdv(nbins,massgrid,gij_2,tabflux_coag,dv,L_k0_2)

   gijnew = gij/3._wp + 2._wp*(gij_2 + dt*L_k0_2)/3._wp


   !check gij values and limit to eps
   do j=1,nbins
      if ( gijnew(j) < 0._wp ) then
         print*,"error in calculating gijnew"
         print*,"i,gij_2,dt*L_k0_2"
         do i=1,nbins
            print*,i,gij_2(i),dt*L_k0_2(i)
         enddo
         print*,"j=",j,", gij =",gijnew(j)
         stop
      else if ( gijnew(j) <= eps  ) then
         gijnew(j) = eps
      endif
   enddo


end subroutine time_solver_k0_kdv




!> @brief Coagulation CFL for DG scheme k>0 piecewise polynomial approximation (simple kernels)
!!
!! CFL formulation from Filbet & Laurencot 2004, dt <= mean_g * dm/dF
!!
!! @param[in]    eps            minimum value for mass distribution approximation gij
!! @param[in]    nbins          number of dust bins
!! @param[in]    kpol           degree of polynomials for approximation
!! @param[in]    massgrid       grid of masses, borders value of mass bins
!! @param[in]    gij            components of g on the polynomial basis
!! @param[in]    tabflux_coag   5D array to evaluate coagulation flux
!! @return                      CFL restriction for coagulation, DG k>0
real(wp) function compute_CFL(eps,nbins,kpol,massgrid,gij,tabflux_coag) result(CFL)
   implicit none
   integer,            intent(in)  :: nbins,kpol
   real(wp),           intent(in)  :: massgrid(nbins+1),eps
   real(wp),           intent(in)  :: gij(nbins,kpol+1)
   real(wp),           intent(in)  :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   
   real(wp) :: flux(nbins),tabdtCFL(nbins)
   real(wp) :: arr_gij(nbins,nbins,kpol+1,kpol+1) 
   real(wp) :: hj
   integer  :: j
   
   call compute_arr_gij(nbins,kpol,gij,arr_gij)
   call compute_flux(nbins,kpol,arr_gij,tabflux_coag,flux)

   ! for debug
   ! do j=1,nbins
   !    print*,"j=",j,",gij =",reshape(gij(j,:),(/kpol+1/)),",flux=",flux(j)
   ! enddo

   ! stop

   tabdtCFL = 0._wp

   !loop in mass domain [massmin,massmax]
   do j=1,nbins
      hj = massgrid(j+1)-massgrid(j)
      
      if (gij(j,1) > eps) then
         if (j==1) then
            !bin 1, where flux(0) = 0, no mass entering the mass domain [massmin,massmax]
            tabdtCFL(j) = abs(gij(j,1)*hj/flux(j))
         else
            tabdtCFL(j) = abs(gij(j,1)*hj/(flux(j)-flux(j-1)))
         endif
      endif

   enddo
      
   CFL = minval(tabdtCFL,tabdtCFL > 0._wp)

   ! check CFL value
   if ( CFL == 0._wp) then
      print*,"j,gij,hj,flux,tabdtCFL"
      do j=1,nbins
         print*,j,reshape(gij(j,:),(/kpol+1/)),massgrid(j+1)-massgrid(j),flux(j),tabdtCFL(j)
      enddo
      print*,"CFL = ",CFL
      print*,"issue in CFL coagulation k > 0, simple kernels"
      stop
   endif


end function






!> @brief Coagulation CFL for DG scheme k>0 piecewise polynomial approximation (physical kernel)
!!
!! CFL formulation from Filbet & Laurencot 2004, dt <= mean_g * dm/dF
!!
!! @param[in]    eps            minimum value for mass distribution approximation gij
!! @param[in]    nbins          number of dust bins
!! @param[in]    kpol           degree of polynomials for approximation
!! @param[in]    massgrid       grid of masses, borders value of mass bins
!! @param[in]    gij            components of g on the polynomial basis
!! @param[in]    tabflux_coag   5D array to evaluate coagulation flux
!! @param[in]    dv             2D array of the differential velocity between grains in bins
!! @return                      CFL restriction for coagulation, DG k>0
real(wp) function compute_CFL_kdv(eps,nbins,kpol,massgrid,gij,tabflux_coag,dv) result(CFL)
   implicit none
   integer,            intent(in)  :: nbins,kpol
   real(wp),           intent(in)  :: massgrid(nbins+1),eps
   real(wp),           intent(in)  :: gij(nbins,kpol+1)
   real(wp),           intent(in)  :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   real(wp),           intent(in)  :: dv(nbins,nbins)
   
   real(wp) :: flux(nbins),tabdtCFL(nbins)
   real(wp) :: arr_gij_dv(nbins,nbins,kpol+1,kpol+1) 
   real(wp) :: hj
   integer  :: j
   
   call compute_arr_gij_kdv(nbins,kpol,gij,dv,arr_gij_dv)
   call compute_flux(nbins,kpol,arr_gij_dv,tabflux_coag,flux)


   tabdtCFL = 0._wp

   !loop in mass domain [massmin,massmax]
   do j=1,nbins
      hj = massgrid(j+1)-massgrid(j)
      
      if (gij(j,1) > eps) then
         if (j==1) then
            !bin 1, where flux(0) = 0, no mass entering the mass domain [massmin,massmax]
            tabdtCFL(j) = abs(gij(j,1)*hj/flux(j))
         else
            tabdtCFL(j) = abs(gij(j,1)*hj/(flux(j)-flux(j-1)))
         endif
      endif

   enddo
      
   CFL = minval(tabdtCFL,tabdtCFL > 0._wp)

   ! check CFL value
   if ( CFL == 0._wp) then
      print*,"j,gij,hj,flux,tabdtCFL"
      do j=1,nbins
         print*,j,reshape(gij(j,:),(/kpol+1/)),massgrid(j+1)-massgrid(j),flux(j),tabdtCFL(j)
      enddo
      print*,"CFL = ",CFL
      print*,"issue in CFL coagulation k > 0, physical kernel cross-section * dv"
      stop
   endif


end function







!> @brief compute the DG operator L for piecewise polynomial approximation and simple kernels (see Lombart et al., 2021) 
!!
!! It is used for the time solver
!!
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[in]    gij               components of g on the polynomial basis
!! @param[in]    tabflux_coag      5D array to evaluate coagulation flux
!! @param[in]    tabintflux_coag   6D array to evaluate the term including the integral of coagulation flux
!! @param[out]   L_k                DG operator for piecewise polynomial approximation in each bin
subroutine operator_L(nbins,kpol,massgrid,mat_coeffs_leg,gij,tabflux_coag,tabintflux_coag,L_k)
   
   implicit none
   integer, intent(in)  :: nbins,kpol
   real(wp),intent(in)  :: massgrid(nbins+1),mat_coeffs_leg(kpol+1,kpol+1)
   real(wp),intent(in)  :: gij(nbins,kpol+1)
   real(wp),intent(in)  :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(in)  :: tabintflux_coag(nbins,kpol+1,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(out) :: L_k(nbins,kpol+1)

   real(wp) :: flux(nbins),intflux(nbins,kpol+1) 
   real(wp) :: arr_gij(nbins,nbins,kpol+1,kpol+1) 
   real(wp) :: c,LegPleft,LegPright,hj
   integer  :: j,k
   real(wp), allocatable :: ak(:)

   call compute_arr_gij(nbins,kpol,gij,arr_gij)
   call compute_flux_intflux(nbins,kpol,arr_gij,tabflux_coag,tabintflux_coag,flux,intflux)

   ! for debug
   ! do j=1,nbins
   !    print*,"j=",j,"flux = ",flux(j)
   ! enddo

   ! do j=1,nbins
   !    print*,"j=",j,"flux = ",flux(j),",intflux =",reshape(intflux(j,:),(/kpol+1/))

   ! enddo
   ! stop
   
   L_k = 0._wp
   do k=0,kpol
      call coeffnorm(k,c)

      allocate(ak(k+1))
      ak = mat_coeffs_leg(k+1,:k+1)

      LegPleft  = phi_pol(k,ak,-1._wp)
      LegPright = phi_pol(k,ak,1._wp)

      do j=1,nbins
         hj = massgrid(j+1)-massgrid(j)
      
         if (j==1) then
            L_k(j,k+1) = (2._wp*(intflux(j,k+1)- flux(j)*LegPright))/(c*hj)
         else
            L_k(j,k+1) = (2._wp*(intflux(j,k+1)-(flux(j)*LegPright - flux(j-1)*LegPleft)))/(c*hj)
         endif
      enddo

      deallocate(ak)

   enddo


end subroutine operator_L




!> @brief compute the DG operator L for piecewise polynomial approximation and physical kernel (see Lombart et al., 2021) 
!!
!! It is used for the time solver
!!
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[in]    gij               components of g on the polynomial basis
!! @param[in]    tabflux_coag      5D array to evaluate coagulation flux
!! @param[in]    tabintflux_coag   6D array to evaluate the term including the integral of coagulation flux
!! @param[in]    dv                2D array of the differential velocity between grains in bins lp and l
!! @param[out]   L_k                DG operator for piecewise polynomial approximation in each bin
subroutine operator_L_kdv(nbins,kpol,massgrid,mat_coeffs_leg,gij,tabflux_coag,tabintflux_coag,dv,L_k)
   
   implicit none
   integer, intent(in)  :: nbins,kpol
   real(wp),intent(in)  :: massgrid(nbins+1),mat_coeffs_leg(kpol+1,kpol+1)
   real(wp),intent(in)  :: gij(nbins,kpol+1)
   real(wp),intent(in)  :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(in)  :: tabintflux_coag(nbins,kpol+1,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(in)  :: dv(nbins,nbins)
   real(wp),intent(out) :: L_k(nbins,kpol+1)

   real(wp) :: flux(nbins),intflux(nbins,kpol+1) 
   real(wp) :: arr_gij_dv(nbins,nbins,kpol+1,kpol+1) 
   real(wp) :: c,LegPleft,LegPright,hj
   integer  :: j,k
   real(wp), allocatable :: ak(:)


   call compute_arr_gij_kdv(nbins,kpol,gij,dv,arr_gij_dv)
   call compute_flux_intflux(nbins,kpol,arr_gij_dv,tabflux_coag,tabintflux_coag,flux,intflux)

   ! for debug
   ! do j=1,nbins
   !    print*,"j=",j,"flux = ",flux(j)
   ! enddo


   ! do j=1,nbins
   !    print*,"j=",j,",intflux =",reshape(intflux(j,:),(/kpol+1/))

   ! enddo
   ! stop
   
   
   L_k = 0._wp
   do k=0,kpol
      call coeffnorm(k,c)

      allocate(ak(k+1))
      ak = mat_coeffs_leg(k+1,:k+1)

      LegPleft  = phi_pol(k,ak,-1._wp)
      LegPright = phi_pol(k,ak,1._wp)

      do j=1,nbins
         hj = massgrid(j+1)-massgrid(j)
      
         if (j==1) then
            L_k(j,k+1) = (2._wp*(intflux(j,k+1)- flux(j)*LegPright))/(c*hj)
         else
            L_k(j,k+1) = (2._wp*(intflux(j,k+1)-(flux(j)*LegPright - flux(j-1)*LegPleft)))/(c*hj)
         endif
      enddo

      deallocate(ak)

   enddo


end subroutine operator_L_kdv





!> @brief compute SSPRK order 3 time solver with piecewise polynomial approximation for simple kernels 
!!
!! See Zhang & Shu 2010 and Lombart et al., 2021
!!
!! @param[in]    eps               minimum value for mass distribution approximation gij
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    massbins          arithmetic mean value of massgrid for each mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[in]    gij               components of g on the polynomial basis
!! @param[in]    tabflux_coag      5D array to evaluate coagulation flux
!! @param[in]    tabintflux_coag   6D array to evaluate the term including the integral of coagulation flux
!! @param[in]    dt                timestep
!! @param[out]   gijnew            evolved components of g on the polynomial basis
subroutine time_solver(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tabflux_coag,tabintflux_coag,dt,gijnew)
   implicit none
   integer, intent(in)    :: nbins,kpol
   real(wp),intent(in)    :: massgrid(nbins+1),massbins(nbins),mat_coeffs_leg(kpol+1,kpol+1)
   real(wp),intent(in)    :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(in)    :: tabintflux_coag(nbins,kpol+1,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(in)    :: dt,eps
   real(wp),intent(in)    :: gij(nbins,kpol+1)
   real(wp),intent(out)   :: gijnew(nbins,kpol+1)


   real(wp) :: gij_1(nbins,kpol+1),gij_2(nbins,kpol+1)
   real(wp) :: L_k(nbins,kpol+1),L_k_1(nbins,kpol+1),L_k_2(nbins,kpol+1)
   real(wp) :: tab_gamma(nbins)
   integer  :: j,k



   !SSPRK3 algo (Zhang & Shu 2010)
   !step 1

   call operator_L(nbins,kpol,massgrid,mat_coeffs_leg,gij,tabflux_coag,tabintflux_coag,L_k)

   gij_1 = gij + dt*L_k


   !apply limiter coefficient to ensure positivity
   call gammafunction(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij_1,tab_gamma)
   do k=1,kpol
      gij_1(:,k+1) = tab_gamma(:)*gij_1(:,k+1)
   enddo


   !check gij values and limit to eps
   do j=1,nbins
      if (gij_1(j,1) < 0._wp) then
         print*,"error in calculating gij_1"
         print*,"j=",j,", gij_1 =",gij_1(j,1)
         stop

      else if ( gij_1(j,1) < eps  ) then
         gij_1(j,1) = eps
         gij_1(j,2:) = 0._wp
      endif
   enddo



   !step 2
   call operator_L(nbins,kpol,massgrid,mat_coeffs_leg,gij_1,tabflux_coag,tabintflux_coag,L_k_1)

   gij_2 = 3._wp*gij/4._wp + (gij_1 + dt*L_k_1)/4._wp

   !apply limiter coefficient to ensure positivity
   call gammafunction(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij_2,tab_gamma)
   do k=1,kpol
      gij_2(:,k+1) = tab_gamma(:)*gij_2(:,k+1)
   enddo

   !check gij values and limit to eps
   do j=1,nbins
      if (gij_2(j,1) < 0._wp) then
         print*,"error in calculating gij_2"
         print*,"j=",j,", gij_2 =",gij_2(j,1)
         stop

      else if ( gij_2(j,1) <= eps  ) then
         gij_2(j,1) = eps
         gij_2(j,2:) = 0._wp

      endif
   enddo


   !step 3
   call operator_L(nbins,kpol,massgrid,mat_coeffs_leg,gij_2,tabflux_coag,tabintflux_coag,L_k_2)
   gijnew = gij/3._wp + 2._wp*(gij_2 + dt*L_k_2)/3._wp

   !apply scale limiter to ensure positivity
   call gammafunction(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gijnew,tab_gamma)
   do k=1,kpol
      gijnew(:,k+1) = tab_gamma(:)*gijnew(:,k+1)
   enddo

   !limit to eps value
   do j=1,nbins
      if (gijnew(j,1) < 0._wp) then
         print*,"error in calculating gijnew"
         print*,"j=",j,", gijnew =",gijnew(j,1)
         stop

      else if ( gijnew(j,1) <= eps  ) then
         gijnew(j,1) = eps
         gijnew(j,2:) = 0._wp

      endif
   enddo


end subroutine time_solver









!> @brief compute SSPRK order 3 time solver with piecewise polynomial approximation for physical kernel
!!
!! See Zhang & Shu 2010 and Lombart et al., 2021
!!
!! @param[in]    eps               minimum value for mass distribution approximation gij
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    massbins          arithmetic mean value of massgrid for each mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[in]    gij               components of g on the polynomial basis
!! @param[in]    tabflux_coag      5D array to evaluate coagulation flux
!! @param[in]    tabintflux_coag   6D array to evaluate the term including the integral of coagulation flux
!! @param[in]    dv                2D array of the differential velocity between grains in bins lp and l
!! @param[in]    dt                timestep
!! @param[out]   gijnew            evolved components of g on the polynomial basis
subroutine time_solver_kdv(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tabflux_coag,tabintflux_coag,dv,dt,gijnew)
   implicit none
   integer, intent(in)  :: nbins,kpol
   real(wp),intent(in)  :: massgrid(nbins+1),massbins(nbins),mat_coeffs_leg(kpol+1,kpol+1)
   real(wp),intent(in)  :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(in)  :: tabintflux_coag(nbins,kpol+1,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(in)  :: dv(nbins,nbins)
   real(wp),intent(in)  :: dt,eps
   real(wp),intent(in)  :: gij(nbins,kpol+1)
   real(wp),intent(out) :: gijnew(nbins,kpol+1)


   real(wp) :: gij_1(nbins,kpol+1),gij_2(nbins,kpol+1)
   real(wp) :: L_k(nbins,kpol+1),L_k_1(nbins,kpol+1),L_k_2(nbins,kpol+1)
   real(wp) :: tab_gamma(nbins)
   integer  :: j,k


   ! call cpu_time(start)

   !SSPRK3 algo (Zhang & Shu 2010)
   !step 1
   call operator_L_kdv(nbins,kpol,massgrid,mat_coeffs_leg,gij,tabflux_coag,tabintflux_coag,dv,L_k)
   gij_1 = gij + dt*L_k


   !apply limiter coefficient to ensure positivity
   call gammafunction(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij_1,tab_gamma)
   do k=1,kpol
      gij_1(:,k+1) = tab_gamma(:)*gij_1(:,k+1)
   enddo


   !check gij values and limit to eps
   do j=1,nbins
      if (gij_1(j,1) < 0._wp) then
         print*,"error in calculating gij_1"
         print*,"j=",j,", gij_1 =",gij_1(j,1)
         stop

      else if ( gij_1(j,1) <= eps  ) then
         gij_1(j,1) = eps
         gij_1(j,2:) = 0._wp
      endif
   enddo



   !step 2
   call operator_L_kdv(nbins,kpol,massgrid,mat_coeffs_leg,gij_1,tabflux_coag,tabintflux_coag,dv,L_k_1)

   gij_2 = 3._wp*gij/4._wp + (gij_1 + dt*L_k_1)/4._wp

   !apply limiter coefficient to ensure positivity
   call gammafunction(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij_2,tab_gamma)

   do k=1,kpol
      gij_2(:,k+1) = tab_gamma(:)*gij_2(:,k+1)
   enddo

   !check gij values and limit to eps
   do j=1,nbins
      if (gij_2(j,1) < 0._wp) then
         print*,"error in calculating gij_2"
         print*,"j=",j,", gij_2 =",gij_2(j,1)
         stop

      else if ( gij_2(j,1) <= eps  ) then
         gij_2(j,1) = eps
         gij_2(j,2:) = 0._wp
      endif
   enddo


   !step 3
   call operator_L_kdv(nbins,kpol,massgrid,mat_coeffs_leg,gij_2,tabflux_coag,tabintflux_coag,dv,L_k_2)
   gijnew = gij/3._wp + 2._wp*(gij_2 + dt*L_k_2)/3._wp

   !apply limiter coefficient to ensure positivity
   call gammafunction(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gijnew,tab_gamma)
   do k=1,kpol
      gijnew(:,k+1) = tab_gamma(:)*gijnew(:,k+1)
   enddo

   !check gij values and limit to eps
   do j=1,nbins
      if (gijnew(j,1) < 0._wp) then
         print*,"error in calculating gijnew"
         print*,"j=",j,", gijnew =",gijnew(j,1)
         stop

      else if ( gijnew(j,1) <= eps  ) then
         gijnew(j,1) = eps
         gijnew(j,2:) = 0._wp

      endif
   enddo


end subroutine time_solver_kdv






!-----------------------------------------------------
! Optimised version to compute functions for DG scheme
!-----------------------------------------------------



!> @brief Coagulation CFL for DG scheme k=0 piecewise constant approximation (simple kernels), with indices optimisation
!!
!! with indices optimisation (where gij > eps)
!!
!! CFL formulation from Filbet & Laurencot 2004, dt <= mean_g * dm/dF
!!
!! @param[in]    eps            minimum value for mass distribution approximation gij
!! @param[in]    nbins          number of dust bins
!! @param[in]    massgrid       grid of masses, borders value of mass bins
!! @param[in]    massbins       arithmetic mean value of massgrid for each mass bins
!! @param[in]    gij            components of g on the polynomial basis
!! @param[in]    tabflux_coag   3D array to evaluate coagulation flux
!! @return                      CFL restriction for coagulation, DG k=0
real(wp) function compute_CFL_k0_opt(eps,nbins,massgrid,massbins,gij,tabflux_coag) result(CFL)
   implicit none
   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins),eps
   real(wp), intent(in)  :: gij(nbins)
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins)
   
   real(wp) :: flux(nbins),tabdtCFL(nbins)
   real(wp) :: hj,arr_gij(nbins,nbins)
   integer  :: j,ind_min,ind_max


   call compute_min_max_indices(eps,nbins,massbins,gij,ind_min,ind_max)
   call compute_arr_gij_k0_opt(nbins,gij,ind_min,ind_max,arr_gij)
   call compute_flux_k0_opt(nbins,ind_min,arr_gij,tabflux_coag,flux)

   ! for debug
   ! !compute min and max indices where flux /= 0
   ! ind_min_flux = minloc(massbins, mask = flux .ne. 0._wp, dim = 1)
   ! ind_max_flux = maxloc(massbins, mask = flux .ne. 0._wp, dim = 1)

   ! do j=1,nbins
   !    print*,"j=",j,",flux =",flux(j)
   ! enddo

   ! do j=1,nbins
   !    print*,"j=",j,", gij = ",gij(j),",flux =",flux(j)
   ! enddo

   ! print*,"ind_min_flux =",ind_min_flux,",ind_max_flux =",ind_max_flux
   ! stop

   tabdtCFL = 0._wp

   !loop in mass domain [massgrid(ind_min),massgrid(ind_max+1)]
   do j=ind_min,ind_max

      hj = massgrid(j+1)-massgrid(j)

      if (j==1) then
         !bin 1, where flux(0) = 0, no mass entering the mass domain
         tabdtCFL(j) = abs(gij(j)*hj/flux(j))
      else
         tabdtCFL(j) = abs(gij(j)*hj/(flux(j)-flux(j-1)))
      endif

   enddo
   
   !CFL calculation only mass domain [massgrid(ind_min),massgrid(ind_max+1)] where flux /= 0.
   CFL = minval(tabdtCFL(ind_min:ind_max))

   
   ! for debug
   ! do j=1,nbins
   !    print*,"j=",j,", tabCFL(j) = ",tabdtCFL(j)
   ! enddo
   ! print*,"CFL = ",CFL
   ! if (tabdtCFL(nbins-4) .ne. 0._wp) stop

   ! check CFL value
   if ( CFL == 0._wp) then
      print*,"j,gij,hj,flux,tabdtCFL"
      do j=1,nbins
         print*,j,gij(j),massgrid(j+1)-massgrid(j),flux(j),tabdtCFL(j)
      enddo
      print*,"CFL = ",CFL
      print*,"issue in CFL coagulation k=0, optimised version"
      stop
   endif

end function


!> @brief Coagulation CFL for DG scheme k=0 piecewise constant approximation (physical kernel), with indices optimisation
!!
!! with indices optimisation (where gij > eps)
!!
!! CFL formulation from Filbet & Laurencot 2004, dt <= mean_g * dm/dF
!!
!! @param[in]    eps            minimum value for mass distribution approximation gij
!! @param[in]    nbins          number of dust bins
!! @param[in]    massgrid       grid of masses, borders value of mass bins
!! @param[in]    massbins       arithmetic mean value of massgrid for each mass bins
!! @param[in]    gij            components of g on the polynomial basis
!! @param[in]    tabflux_coag   3D array to evaluate coagulation flux
!! @param[in]    dv             2D array of the differential velocity between grains in bins
!! @return                      CFL restriction for coagulation, DG k=0
real(wp) function compute_CFL_k0_kdv_opt(eps,nbins,massgrid,massbins,gij,tabflux_coag,dv) result(CFL)
   implicit none
   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins)
   real(wp), intent(in)  :: gij(nbins),eps
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins),dv(nbins,nbins)

   
   real(wp) :: flux(nbins),tabdtCFL(nbins)
   real(wp) :: arr_gij_dv(nbins,nbins)
   real(wp) :: hj
   integer  :: j,ind_min,ind_max


   call compute_min_max_indices(eps,nbins,massbins,gij,ind_min,ind_max)
   call compute_arr_gij_k0_kdv_opt(nbins,gij,ind_min,ind_max,dv,arr_gij_dv)
   call compute_flux_k0_opt(nbins,ind_min,arr_gij_dv,tabflux_coag,flux)

   ! for debug
   ! !compute min and max indices where flux /= 0
   ! ind_min_flux = minloc(massbins, mask = flux .ne. 0._wp, dim = 1)
   ! ind_max_flux = maxloc(massbins, mask = flux .ne. 0._wp, dim = 1)

   ! do j=1,nbins
   !    print*,"j=",j,",flux =",flux(j)
   ! enddo

   ! do j=1,nbins
   !    print*,"j=",j,", gij = ",gij(j),",flux =",flux(j)
   ! enddo

   ! print*,"ind_min =",ind_min,",ind_max =",ind_max
   ! print*,"ind_min_flux =",ind_min_flux,",ind_max_flux =",ind_max_flux
   ! stop

   tabdtCFL = 0._wp

   !loop in mass domain [massgrid(ind_min),massgrid(ind_max+1)]
   do j=ind_min,ind_max

      hj = massgrid(j+1)-massgrid(j)

      if (j==1) then
         !bin 1, where flux(0) = 0, no mass entering the mass domain
         tabdtCFL(j) = abs(gij(j)*hj/flux(j))
      else
         tabdtCFL(j) = abs(gij(j)*hj/(flux(j)-flux(j-1)))
      endif

   enddo
   
   !CFL calculation only mass domain [massgrid(ind_min),massgrid(ind_max+1)] where flux .ne. 0.
   CFL = minval(tabdtCFL(ind_min:ind_max))

   ! for debug
   ! do j=1,nbins
   !    print*,"j=",j,", tabCFL(j) = ",tabdtCFL(j)
   ! enddo
   ! print*,"CFL = ",CFL
   ! stop

   ! check CFL value
   if ( CFL == 0._wp) then
      print*,"j,gij,hj,flux,tabdtCFL"
      do j=1,nbins
         print*,j,gij(j),massgrid(j+1)-massgrid(j),flux(j),tabdtCFL(j)
      enddo
      print*,"CFL = ",CFL
      print*,"issue in CFL coagulation k=0, physical kernel cross-section * dv,  optimised version"
      stop
   endif



end function


!> @brief compute the DG operator L for piecewise constant approximation and simple kernels (see Lombart et al., 2021)
!!
!! with indices optimisation (where gij > eps)
!!
!! It is used for the time solver
!!
!! @param[in]    eps            minimum value for mass distribution approximation gij
!! @param[in]    nbins          number of dust bins
!! @param[in]    massgrid       grid of masses, borders value of mass bins
!! @param[in]    massbins       arithmetic mean value of massgrid for each mass bins
!! @param[in]    gij            components of g on the polynomial basis
!! @param[in]    tabflux_coag   3D array to evaluate coagulation flux
!! @param[out]   L_k0            DG operator for piecewise constant approximation in each bin
subroutine operator_L_k0_opt(eps,nbins,massgrid,massbins,gij,tabflux_coag,L_k0)
   implicit none
   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins),eps
   real(wp), intent(in)  :: gij(nbins)
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins)
   real(wp), intent(out) :: L_k0(nbins)

   real(wp) :: flux(nbins)
   real(wp) :: arr_gij(nbins,nbins)
   real(wp) :: hj
   integer  :: j
   integer  :: ind_min,ind_max
   integer  :: ind_min_flux,ind_max_flux

   
   call compute_min_max_indices(eps,nbins,massbins,gij,ind_min,ind_max)
   call compute_arr_gij_k0_opt(nbins,gij,ind_min,ind_max,arr_gij)
   call compute_flux_k0_opt(nbins,ind_min,arr_gij,tabflux_coag,flux)

   ! for debug
   ! do j=1,nbins
   !    print*,"j=",j,",flux =",flux(j)
   ! enddo

   ! do j=1,nbins
   !    print*,"j=",j,",gij(j) =",gij(j)
   ! enddo

   ! do j=1,nbins
   !    print*,"j=",j,",gij(j) =",gij(j),", flux =",flux(j)
   ! enddo


   !compute min and max indices where flux /= 0
   ind_min_flux = minloc(massbins, mask = flux .ne. 0._wp, dim = 1)
   ind_max_flux = maxloc(massbins, mask = flux .ne. 0._wp, dim = 1)

   ! print*,"ind_min_flux=",ind_min_flux,",ind_max_flux=",ind_max_flux

   L_k0 = 0._wp
   !loop on indices where flux /= 0
   !ind_max_flux+1 needed to calcul dF when flux(ind_max_flux+1) = 0
   do j=ind_min_flux,ind_max_flux+1
   
      hj = massgrid(j+1)-massgrid(j)
      if (j == 1) then
         !bin 1 with flux(0) = 0, no mass entering the mass domain [massmin,massmax]
         L_k0(j) = -flux(j)/hj
      else

         ! L_k0(j) = -(flux(j) - flux(j-1))/hj

         !normalisation to reduce errors from substraction operation
         L_k0(j) = -(flux(j)/hj - flux(j-1)/hj)
      endif


   enddo


end subroutine operator_L_k0_opt


!> @brief compute the DG operator L for piecewise constant approximation and physical kernel
!!
!! with indices optimisation (only with gij > eps)
!!
!! It is used for the time solver
!!
!! @param[in]    eps            minimum value for mass distribution approximation gij
!! @param[in]    nbins          number of dust bins
!! @param[in]    massgrid       grid of masses, borders value of mass bins
!! @param[in]    massbins       arithmetic mean value of massgrid for each mass bins
!! @param[in]    gij            components of g on the polynomial basis
!! @param[in]    tabflux_coag   3D array to evaluate coagulation flux
!! @param[in]    dv             2D array of the differential velocity between grains in bins
!! @param[out]   L_k0            DG operator for piecewise constant approximation in each bin
subroutine operator_L_k0_kdv_opt(eps,nbins,massgrid,massbins,gij,tabflux_coag,dv,L_k0)
   implicit none
   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins)
   real(wp), intent(in)  :: gij(nbins),eps
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins),dv(nbins,nbins)
   real(wp), intent(out) :: L_k0(nbins)

   real(wp) :: flux(nbins)
   real(wp) :: arr_gij_dv(nbins,nbins)
   real(wp) :: hj
   integer  :: j
   integer  :: ind_min,ind_max
   integer  :: ind_min_flux,ind_max_flux
   
   call compute_min_max_indices(eps,nbins,massbins,gij,ind_min,ind_max)
   call compute_arr_gij_k0_kdv_opt(nbins,gij,ind_min,ind_max,dv,arr_gij_dv)
   call compute_flux_k0_opt(nbins,ind_min,arr_gij_dv,tabflux_coag,flux)
   
   ! for debug
   ! do j=1,nbins
   !    print*,"j=",j,",flux =",flux(j)
   ! enddo

   ! do j=1,nbins
   !    print*,"j=",j,",gij(j) =",gij(j)
   ! enddo

   ! do j=1,nbins
   !    print*,"j=",j,",gij(j) =",gij(j),", flux =",flux(j)
   ! enddo
   ! stop
   
   !compute min and max indices where flux /= 0
   ind_min_flux = minloc(massbins, mask = flux .ne. 0._wp, dim = 1)
   ind_max_flux = maxloc(massbins, mask = flux .ne. 0._wp, dim = 1)

   ! print*,"ind_min_flux=",ind_min_flux,",ind_max_flux=",ind_max_flux
   ! stop
   
   L_k0 = 0._wp
   !loop on indices where flux /= 0
   !ind_max_flux+1 needed to calcul dF when flux(ind_max_flux+1) = 0
   do j=ind_min_flux,ind_max_flux+1
   
      hj = massgrid(j+1)-massgrid(j)
      if (j == 1) then
         !bin 1 with flux(0) = 0, no mass entering the mass domain [massmin,massmax]
         L_k0(j) = -flux(j)/hj
      else
         !original expression
         ! L_k0(j) = -(flux(j) - flux(j-1))/hj

         !normalisation to reduce errors from substraction
         L_k0(j) = -(flux(j)/hj - flux(j-1)/hj)
      endif
   enddo

end subroutine operator_L_k0_kdv_opt


!> @brief Function to compute SSPRK order 3 time solver with piecewise constant approximation for simple kernels
!!
!! with indices optimisation (only with gij > eps)
!!
!! @param[in]     eps            minimum value for mass distribution approximation gij
!! @param[in]     nbins          number of dust bins
!! @param[in]     massgrid       grid of masses, borders value of mass bins
!! @param[in]     massbins       arithmetic mean value of massgrid for each mass bins
!! @param[in]     gij            components of g on the polynomial basis
!! @param[in]     tabflux_coag   3D array to evaluate coagulation flux
!! @param[in]     dt             timestep
!! @param[out]    gijnew         evolved components of g on the polynomial basis
subroutine time_solver_k0_opt(eps,nbins,massgrid,massbins,gij,tabflux_coag,dt,gijnew)

   implicit none

   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins)
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins)
   real(wp), intent(in)  :: dt,eps
   real(wp), intent(in)  :: gij(nbins)
   real(wp), intent(out) :: gijnew(nbins)


   real(wp) :: gij_1(nbins),gij_2(nbins)
   real(wp) :: L_k0(nbins),L_k0_1(nbins),L_k0_2(nbins)

   integer :: j,i


   !SSPRK3 algorithm (Zhang & Shu 2010)
   !step 1

   call operator_L_k0_opt(eps,nbins,massgrid,massbins,gij,tabflux_coag,L_k0)

   gij_1 = gij + dt*L_k0

   !check gij values and limit to eps 
   do j=1,nbins
      if ( gij_1(j) < 0._wp  ) then
         print*,"error in calculating gij_1"
         print*,"i,gij,dt*Lk0"
         do i=1,nbins
            print*,i,gij(i),dt*L_k0(i)
         enddo
         print*,"j=",j,", gij_1 =",gij_1(j)
         stop
      else if ( gij_1(j) <= eps  ) then
         gij_1(j) = eps
      endif
   enddo



   !step 2
   call operator_L_k0_opt(eps,nbins,massgrid,massbins,gij_1,tabflux_coag,L_k0_1)

   gij_2 = 3._wp*gij/4._wp + (gij_1 + dt*L_k0_1)/4._wp


   !check gij values and limit to eps
   do j=1,nbins
      if ( gij_2(j) < 0._wp  ) then
         print*,"error in calculating gij_2"
         print*,"i,gij_1,dt*L_k0_2"
         do i=1,nbins
            print*,i,gij_1(i),dt*L_k0_1(i)
         enddo
         print*,"j=",j,", gij_2 =",gij_2(j)
         stop

      else if ( gij_2(j) <= eps  ) then
         gij_2(j) = eps
      endif
   enddo


   !step 3
   call operator_L_k0_opt(eps,nbins,massgrid,massbins,gij_2,tabflux_coag,L_k0_2)

   gijnew = gij/3._wp + 2._wp*(gij_2 + dt*L_k0_2)/3._wp

   !check gij values and limit to eps 
   do j=1,nbins
      if ( gijnew(j) < 0._wp ) then
         print*,"error in calculating gijnew"
         print*,"i,gij_2,dt*Lk0_2"
         do i=1,nbins
            print*,i,gij_2(i),dt*L_k0_2(i)
         enddo
         print*,"j=",j,", gijnew =",gijnew(j)
         stop
      else if ( gijnew(j) <= eps  ) then
         gijnew(j) = eps
      endif
   enddo

end subroutine time_solver_k0_opt


!> @brief Function to compute SSPRK order 3 time solver with piecewise constant approximation for physical kernel
!!
!! with indices optimisation (only with gij > eps) 
!!
!! See Zhang & Shu 2010 and Lombart et al., 2021
!!
!! @param[in]     eps            minimum value for mass distribution approximation gij
!! @param[in]     nbins          number of dust bins
!! @param[in]     massgrid       grid of masses, borders value of mass bins
!! @param[in]     massbins       arithmetic mean value of massgrid for each mass bins
!! @param[in]     gij            components of g on the polynomial basis
!! @param[in]     tabflux_coag   3D array to evaluate coagulation flux
!! @param[in]     dv             2D array of the differential velocity between grains in bins lp and l
!! @param[in]     dt             timestep
!! @param[out]    gijnew         evolved components of g on the polynomial basis
subroutine time_solver_k0_kdv_opt(eps,nbins,massgrid,massbins,gij,tabflux_coag,dv,dt,gijnew)

   implicit none

   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins)
   real(wp), intent(in)  :: tabflux_coag(nbins,nbins,nbins),dv(nbins,nbins)
   real(wp), intent(in)  :: dt,eps
   real(wp), intent(in)  :: gij(nbins)
   real(wp), intent(out) :: gijnew(nbins)


   real(wp) :: gij_1(nbins),gij_2(nbins)
   real(wp) :: L_k0(nbins),L_k0_1(nbins),L_k0_2(nbins)

   integer :: j,i


   !SSPRK3 algo (Zhang & Shu 2010)
   !step 1
   call operator_L_k0_kdv_opt(eps,nbins,massgrid,massbins,gij,tabflux_coag,dv,L_k0)

   gij_1 = gij + dt*L_k0


   !check gij values and limit to eps
   do j=1,nbins
      if ( gij_1(j) < 0._wp  ) then
         print*,"error in calculating gij_1"
         print*,"i,gij,dt*L_k0"
         do i=1,nbins
            print*,i,gij(i),dt*L_k0(i)
         enddo
         print*,"j=",j,", gij_1 =",gij_1(j)
         stop
      else if ( gij_1(j) <= eps  ) then
         gij_1(j) = eps

      endif
   enddo


   !step 2
   call operator_L_k0_kdv_opt(eps,nbins,massgrid,massbins,gij_1,tabflux_coag,dv,L_k0_1)

   gij_2 = 3._wp*gij/4._wp + (gij_1 + dt*L_k0_1)/4._wp

   !check gij values and limit to eps
   do j=1,nbins
      if ( gij_2(j) < 0._wp  ) then
         print*,"error in calculating gij_2"
         print*,"i,gij_1,dt*L_k0_2"
         do i=1,nbins
            print*,i,gij_1(i),dt*L_k0_1(i)
         enddo
         print*,"j=",j,", gij_2 =",gij_2(j)
         stop

      else if ( gij_2(j) <= eps  ) then
         gij_2(j) = eps
      endif
   enddo


   !step 3
   call operator_L_k0_kdv_opt(eps,nbins,massgrid,massbins,gij_2,tabflux_coag,dv,L_k0_2)

   gijnew = gij/3._wp + 2._wp*(gij_2 + dt*L_k0_2)/3._wp

   !check gij values and limit to eps
   do j=1,nbins
      if ( gijnew(j) < 0._wp ) then
         print*,"error in calculating gijnew"
         print*,"i,gij_2,dt*L_k0_2"
         do i=1,nbins
            print*,i,gij_2(i),dt*L_k0_2(i)
         enddo
         print*,"j=",j,", gij =",gijnew(j)
         stop
      else if ( gijnew(j) <= eps  ) then
         gijnew(j) = eps
      endif
   enddo


end subroutine time_solver_k0_kdv_opt


!> @brief Coagulation CFL for DG scheme k>0 piecewise polynomial approximation (simple kernels)
!!
!! with indices optimisation (only with gij > eps) 
!!
!! CFL formulation from Filbet & Laurencot 2004, dt <= mean_g * dm/dF
!!
!! @param[in]    eps            minimum value for mass distribution approximation gij
!! @param[in]    nbins          number of dust bins
!! @param[in]    kpol           degree of polynomials for approximation
!! @param[in]    massgrid       grid of masses, borders value of mass bins
!! @param[in]    massbins       arithmetic mean value of massgrid for each mass bins
!! @param[in]    gij            components of g on the polynomial basis
!! @param[in]    tabflux_coag   5D array to evaluate coagulation flux
!! @return                      CFL restriction for coagulation, DG k>0
real(wp) function compute_CFL_opt(eps,nbins,kpol,massgrid,massbins,gij,tabflux_coag) result(CFL)
   implicit none
   integer,            intent(in)  :: nbins,kpol
   real(wp),           intent(in)  :: massgrid(nbins+1),massbins(nbins),eps
   real(wp),           intent(in)  :: gij(nbins,kpol+1)
   real(wp),           intent(in)  :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   
   real(wp) :: flux(nbins),tabdtCFL(nbins)
   real(wp) :: arr_gij(nbins,nbins,kpol+1,kpol+1) 
   real(wp) :: hj
   integer  :: j,ind_min,ind_max

   call compute_min_max_indices(eps,nbins,massbins,gij(:,1),ind_min,ind_max)
   call compute_arr_gij_opt(nbins,kpol,gij,ind_min,ind_max,arr_gij)
   call compute_flux_opt(nbins,kpol,ind_min,arr_gij,tabflux_coag,flux)

   ! for debug
   ! !compute min and max indices where flux /= 0
   ! ind_min_flux = minloc(massbins, mask = flux .ne. 0._wp, dim = 1)
   ! ind_max_flux = maxloc(massbins, mask = flux .ne. 0._wp, dim = 1)


   ! do j=1,nbins
   !    print*,"j=",j,",gij =",reshape(gij(j,:),(/kpol+1/)),",flux=",flux(j)
   ! enddo

   ! print*,"ind_min =",ind_min,", ind_max = ",ind_max


   tabdtCFL = 0._wp

   !loop in mass domain [massgrid(ind_min),massgrid(ind_max+1)]
   do j=ind_min,ind_max

      hj = massgrid(j+1)-massgrid(j)

      if (j==1) then
         !bin 1, where flux(0) = 0, no mass entering the mass domain
         tabdtCFL(j) = abs(gij(j,1)*hj/flux(j))
      else
         tabdtCFL(j) = abs(gij(j,1)*hj/(flux(j)-flux(j-1)))
      endif
   enddo
   
   !CFL calculation only mass domain [massgrid(ind_min),massgrid(ind_max+1)] where gij > eps
   CFL = minval(tabdtCFL(ind_min:ind_max))

   ! check CFL value
   if ( CFL == 0._wp) then
      print*,"j,gij,hj,flux,tabdtCFL"
      do j=1,nbins
         print*,j,reshape(gij(j,:),(/kpol+1/)),massgrid(j+1)-massgrid(j),flux(j),tabdtCFL(j)
      enddo
      print*,"CFL = ",CFL
      print*,"issue in CFL coagulation k > 0, simple kernels, optimised version"
      stop
   endif

end function


!> @brief Coagulation CFL for DG scheme k>0 piecewise polynomial approximation (physical kernel)
!!
!! with indices optimisation (only with gij > eps) 
!!
!! CFL formulation from Filbet & Laurencot 2004, dt <= mean_g * dm/dF
!!
!! @param[in]    eps            minimum value for mass distribution approximation gij
!! @param[in]    nbins          number of dust bins
!! @param[in]    kpol           degree of polynomials for approximation
!! @param[in]    massgrid       grid of masses, borders value of mass bins
!! @param[in]    massbins       arithmetic mean value of massgrid for each mass bins
!! @param[in]    gij            components of g on the polynomial basis
!! @param[in]    tabflux_coag   5D array to evaluate coagulation flux
!! @param[in]    dv             2D array of the differential velocity between grains in bins
!! @return                      CFL restriction for coagulation, DG k>0
real(wp) function compute_CFL_kdv_opt(eps,nbins,kpol,massgrid,massbins,gij,tabflux_coag,dv) result(CFL)
   implicit none
   integer,            intent(in)  :: nbins,kpol
   real(wp),           intent(in)  :: massgrid(nbins+1),massbins(nbins)
   real(wp),           intent(in)  :: gij(nbins,kpol+1),eps
   real(wp),           intent(in)  :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   real(wp),           intent(in)  :: dv(nbins,nbins)
   
   real(wp) :: flux(nbins),tabdtCFL(nbins)
   real(wp) :: arr_gij_dv(nbins,nbins,kpol+1,kpol+1) 
   real(wp) :: hj
   integer  :: j,ind_min,ind_max
   
   call compute_min_max_indices(eps,nbins,massbins,gij(:,1),ind_min,ind_max)
   call compute_arr_gij_kdv_opt(nbins,kpol,gij,ind_min,ind_max,dv,arr_gij_dv)
   call compute_flux_opt(nbins,kpol,ind_min,arr_gij_dv,tabflux_coag,flux)
   

   tabdtCFL = 0._wp

   !loop in mass domain [massgrid(ind_min),massgrid(ind_max)] where gij > eps
   do j=ind_min,ind_max

      hj = massgrid(j+1)-massgrid(j)

      if (j==1) then
         !bin 1, where flux(0) = 0, no mass entering the mass domain
         tabdtCFL(j) = abs(gij(j,1)*hj/flux(j))
      else
         tabdtCFL(j) = abs(gij(j,1)*hj/(flux(j)-flux(j-1)))
      endif
   enddo
   
   !CFL calculation only mass domain [massgrid(ind_min),massgrid(ind_max+1)] where gij > eps
   CFL = minval(tabdtCFL(ind_min:ind_max))

   ! check CFL value
   if ( CFL == 0._wp) then
      print*,"j,gij,hj,flux,tabdtCFL"
      do j=1,nbins
         print*,j,reshape(gij(j,:),(/kpol+1/)),massgrid(j+1)-massgrid(j),flux(j),tabdtCFL(j)
      enddo
      print*,"CFL = ",CFL
      print*,"issue in CFL coagulation k > 0, physical kernel cross-section * dv, optimised version"
      stop
   endif

end function


!> @brief compute the DG operator L for piecewise polynomial approximation and simple kernels (see Lombart et al., 2021)
!!
!! with indices optimisation (only with gij > eps)  
!!
!! It is used for the time solver
!!
!! @param[in]    eps               minimum value for mass distribution approximation gij
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    massbins          arithmetic mean value of massgrid for each mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[in]    gij               components of g on the polynomial basis
!! @param[in]    tabflux_coag      5D array to evaluate coagulation flux
!! @param[in]    tabintflux_coag   6D array to evaluate the term including the integral of coagulation flux
!! @param[out]   L_k                DG operator for piecewise polynomial approximation in each bin
subroutine operator_L_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tabflux_coag,tabintflux_coag,L_k)
   
   implicit none
   integer, intent(in)  :: nbins,kpol
   real(wp),intent(in)  :: massgrid(nbins+1),massbins(nbins),mat_coeffs_leg(kpol+1,kpol+1),eps
   real(wp),intent(in)  :: gij(nbins,kpol+1)
   real(wp),intent(in)  :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(in)  :: tabintflux_coag(nbins,kpol+1,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(out) :: L_k(nbins,kpol+1)

   real(wp) :: flux(nbins),intflux(nbins,kpol+1)
   real(wp) :: c,LegPleft,LegPright,hj
   real(wp) :: arr_gij(nbins,nbins,kpol+1,kpol+1) 
   integer  :: j,k,ind_min,ind_max
   ! integer  :: ind_min_flux,ind_max_flux
   integer  :: ind_min_intflux,ind_max_intflux
   real(wp), allocatable :: ak(:)

   call compute_min_max_indices(eps,nbins,massbins,gij(:,1),ind_min,ind_max)
   call compute_arr_gij_opt(nbins,kpol,gij,ind_min,ind_max,arr_gij)
   call compute_flux_intflux_opt(nbins,kpol,ind_min,arr_gij,tabflux_coag,tabintflux_coag,flux,intflux)

   ! for debug
   ! do j=1,nbins
   !    print*,"j=",j,"gij =",reshape(gij(j,:),(/kpol+1/))
   ! enddo
   

   ! do j=1,nbins
   !    print*,"j=",j,"flux = ",flux(j),",intflux =",reshape(intflux(j,:),(/kpol+1/))
   ! enddo

   ! print*,"ind_min=",ind_min,", ind_max=",ind_max



   !compute min and max indices where flux /= 0
   ! ind_min_flux   = minloc(massbins, mask = flux .ne. 0._wp, dim = 1)
   ! ind_max_flux   = maxloc(massbins, mask = flux .ne. 0._wp, dim = 1)


   !compute min and max indices where intflux /= 0
   ind_min_intflux   = minloc(massbins, mask = intflux(:,kpol+1) .ne. 0._wp, dim = 1)
   ind_max_intflux   = maxloc(massbins, mask = intflux(:,kpol+1) .ne. 0._wp, dim = 1)

   ! for debug
   ! print*,"ind_min_flux=",ind_min_flux,", ind_max_flux=",ind_max_flux
   ! print*,"ind_min_intflux=",ind_min_intflux,", ind_max_intflux=",ind_max_intflux

   ! stop

   L_k = 0._wp
   do k=0,kpol
      call coeffnorm(k,c)
      
      allocate(ak(k+1))
      ak = mat_coeffs_leg(k+1,:k+1)

      LegPleft  = phi_pol(k,ak,-1._wp)
      LegPright = phi_pol(k,ak,1._wp)

      do j=ind_min_intflux,ind_max_intflux
         hj = massgrid(j+1)-massgrid(j)
         if (j==1)  then
            L_k(j,k+1) = (2._wp*(intflux(j,k+1)- flux(j)*LegPright))/(c*hj)
         else
            L_k(j,k+1) = (2._wp*(intflux(j,k+1)-(flux(j)*LegPright - flux(j-1)*LegPleft)))/(c*hj)
         endif
      enddo

      deallocate(ak)
   enddo

end subroutine operator_L_opt


!> @brief compute the DG operator L for piecewise polynomial approximation and physical kernel (see Lombart et al., 2021)
!!
!! with indices optimisation (only with gij > eps) 
!!
!! It is used for the time solver
!!
!! @param[in]    eps               minimum value for mass distribution approximation gij
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    massbins          arithmetic mean value of massgrid for each mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[in]    gij               components of g on the polynomial basis
!! @param[in]    tabflux_coag      5D array to evaluate coagulation flux
!! @param[in]    tabintflux_coag   6D array to evaluate the term including the integral of coagulation flux
!! @param[in]    dv                2D array of the differential velocity between grains in bins lp and l
!! @param[out]   L_k                DG operator for piecewise polynomial approximation in each bin
subroutine operator_L_kdv_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tabflux_coag,tabintflux_coag,dv,L_k)
   
   implicit none
   integer, intent(in)  :: nbins,kpol
   real(wp),intent(in)  :: massgrid(nbins+1),massbins(nbins),mat_coeffs_leg(kpol+1,kpol+1)
   real(wp),intent(in)  :: gij(nbins,kpol+1),eps
   real(wp),intent(in)  :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(in)  :: tabintflux_coag(nbins,kpol+1,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(in)  :: dv(nbins,nbins)
   real(wp),intent(out) :: L_k(nbins,kpol+1)

   real(wp) :: flux(nbins),intflux(nbins,kpol+1) 
   real(wp) :: arr_gij_dv(nbins,nbins,kpol+1,kpol+1) 
   real(wp) :: c,LegPleft,LegPright,hj
   integer  :: j,k,ind_min,ind_max
   ! integer  :: ind_min_flux,ind_max_flux
   integer  :: ind_min_intflux,ind_max_intflux
   real(wp), allocatable :: ak(:)

   call compute_min_max_indices(eps,nbins,massbins,gij(:,1),ind_min,ind_max)
   call compute_arr_gij_kdv_opt(nbins,kpol,gij,ind_min,ind_max,dv,arr_gij_dv)
   call compute_flux_intflux_opt(nbins,kpol,ind_min,arr_gij_dv,tabflux_coag,tabintflux_coag,flux,intflux)

   !for debug
   ! do j=1,nbins
   !    print*,"j=",j,"gij =",reshape(gij(j,:),(/kpol+1/))
   ! enddo
   

   ! do j=1,nbins
   !    print*,"j=",j,"flux = ",flux(j),",intflux =",reshape(intflux(j,:),(/kpol+1/))
   ! enddo

   ! print*,"ind_min=",ind_min,", ind_max=",ind_max



   !compute min and max indices where flux /= 0
   ! ind_min_flux   = minloc(massbins, mask = flux .ne. 0._wp, dim = 1)
   ! ind_max_flux   = maxloc(massbins, mask = flux .ne. 0._wp, dim = 1)

   

   !compute min and max indices where flux /= 0
   ind_min_intflux   = minloc(massbins, mask = intflux(:,kpol+1) .ne. 0._wp, dim = 1)
   ind_max_intflux   = maxloc(massbins, mask = intflux(:,kpol+1) .ne. 0._wp, dim = 1)

   ! for debug
   ! print*,"ind_min_flux=",ind_min_flux,", ind_max_flux=",ind_max_flux
   ! print*,"ind_min_intflux=",ind_min_intflux,", ind_max_intflux=",ind_max_intflux
   

   ! stop

   L_k = 0._wp
   do k=0,kpol
      call coeffnorm(k,c)
      
      allocate(ak(k+1))
      ak = mat_coeffs_leg(k+1,:k+1)

      LegPleft  = phi_pol(k,ak,-1._wp)
      LegPright = phi_pol(k,ak,1._wp)

      do j=ind_min_intflux,ind_max_intflux
         hj = massgrid(j+1)-massgrid(j)
         if (j==1)  then
            L_k(j,k+1) = (2._wp*(intflux(j,k+1)- flux(j)*LegPright))/(c*hj)
         else
            L_k(j,k+1) = (2._wp*(intflux(j,k+1)-(flux(j)*LegPright - flux(j-1)*LegPleft)))/(c*hj)
         endif
      enddo

      deallocate(ak)
   enddo

end subroutine operator_L_kdv_opt


!> @brief compute SSPRK order 3 time solver with piecewise polynomial approximation for simple kernels
!!
!! with indices optimisation (only with gij > eps) 
!!
!! See Zhang & Shu 2010 and Lombart et al., 2021
!!
!! @param[in]    eps               minimum value for mass distribution approximation gij
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    massbins          arithmetic mean value of massgrid for each mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[in]    gij               components of g on the polynomial basis
!! @param[in]    tabflux_coag      5D array to evaluate coagulation flux
!! @param[in]    tabintflux_coag   6D array to evaluate the term including the integral of coagulation flux
!! @param[in]    dt                timestep
!! @param[out]   gijnew            evolved components of g on the polynomial basis
subroutine time_solver_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tabflux_coag,tabintflux_coag,dt,gijnew)
   implicit none
   integer, intent(in)   :: nbins,kpol
   real(wp),intent(in)   :: massgrid(nbins+1),massbins(nbins),mat_coeffs_leg(kpol+1,kpol+1)
   real(wp),intent(in)   :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(in)   :: tabintflux_coag(nbins,kpol+1,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(in)   :: dt,eps
   real(wp),intent(in)   :: gij(nbins,kpol+1)
   real(wp),intent(out)  :: gijnew(nbins,kpol+1)


   real(wp) :: gij_1(nbins,kpol+1),gij_2(nbins,kpol+1)
   real(wp) :: L_k(nbins,kpol+1),L_k_1(nbins,kpol+1),L_k_2(nbins,kpol+1)
   real(wp) :: tab_gamma(nbins)
   integer  :: j,k


   !SSPRK3 algo (Zhang & Shu 2010)
   !step 1
   call operator_L_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tabflux_coag,tabintflux_coag,L_k)
   gij_1 = gij + dt*L_k


   !apply limiter coefficient to ensure positivity
   call gammafunction_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij_1,tab_gamma)
   do k=1,kpol
      gij_1(:,k+1) = tab_gamma(:)*gij_1(:,k+1)
   enddo

   !check gij values and limit to eps
   do j=1,nbins
      if (gij_1(j,1) < 0._wp) then
         print*,"error in calculating gij_1"
         print*,"j=",j,", gij_1 =",gij_1(j,1)
         stop

      else if ( gij_1(j,1) <= eps  ) then
         gij_1(j,1) = eps
         gij_1(j,2:) = 0._wp
      endif
   enddo

   

   !step 2
   call operator_L_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij_1,tabflux_coag,tabintflux_coag,L_k_1)

   gij_2 = 3._wp*gij/4._wp + (gij_1 + dt*L_k_1)/4._wp


   !apply limiter coefficient to ensure positivity
   call gammafunction_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij_2,tab_gamma)
   do k=1,kpol
      gij_2(:,k+1) = tab_gamma(:)*gij_2(:,k+1)
   enddo



   !check gij values and limit to eps
   do j=1,nbins
      if (gij_2(j,1) < 0._wp) then
         print*,"error in calculating gij_2"
         print*,"j=",j,", gij_2 =",gij_2(j,1)
         stop

      else if ( gij_2(j,1) <= eps  ) then
         gij_2(j,1) = eps
         gij_2(j,2:) = 0._wp
      endif
   enddo


   !step 3
   call operator_L_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij_2,tabflux_coag,tabintflux_coag,L_k_2)
   gijnew = gij/3._wp + 2._wp*(gij_2 + dt*L_k_2)/3._wp


   !apply limiter coefficient to ensure positivity
   call gammafunction_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gijnew,tab_gamma)

   do k=1,kpol
      gijnew(:,k+1) = tab_gamma(:)*gijnew(:,k+1)
   enddo


   !check gij values and limit to eps
   do j=1,nbins
      if (gijnew(j,1) < 0._wp) then
         print*,"error in calculating gijnew"
         print*,"j=",j,", gijnew =",gijnew(j,1)
         stop

      else if ( gijnew(j,1) <= eps  ) then
         gijnew(j,1) = eps
         gijnew(j,2:) = 0._wp
      endif
   enddo


end subroutine time_solver_opt

!> @brief compute SSPRK order 3 time solver with piecewise polynomial approximation for physical kernel
!! 
!! with indices optimisation (only with gij > eps) 
!!
!! See Zhang & Shu 2010 and Lombart et al., 2021
!!
!! @param[in]    eps               minimum value for mass distribution approximation gij
!! @param[in]    nbins             number of dust bins
!! @param[in]    kpol              degree of polynomials for approximation
!! @param[in]    massgrid          grid of masses, borders value of mass bins
!! @param[in]    massbins          arithmetic mean value of massgrid for each mass bins
!! @param[in]    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
!! @param[in]    gij               components of g on the polynomial basis
!! @param[in]    tabflux_coag      5D array to evaluate coagulation flux
!! @param[in]    tabintflux_coag   6D array to evaluate the term including the integral of coagulation flux
!! @param[in]    dv                2D array of the differential velocity between grains in bins lp and l
!! @param[in]    dt                timestep
!! @param[out]   gijnew            evolved components of g on the polynomial basis
subroutine time_solver_kdv_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tabflux_coag,tabintflux_coag,dv,dt,gijnew)
   implicit none
   integer, intent(in)  :: nbins,kpol
   real(wp),intent(in)  :: massgrid(nbins+1),massbins(nbins),mat_coeffs_leg(kpol+1,kpol+1)
   real(wp),intent(in)  :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(in)  :: tabintflux_coag(nbins,kpol+1,nbins,nbins,kpol+1,kpol+1)
   real(wp),intent(in)  :: dv(nbins,nbins)
   real(wp),intent(in)  :: dt,eps
   real(wp),intent(in)  :: gij(nbins,kpol+1)
   real(wp),intent(out) :: gijnew(nbins,kpol+1)


   real(wp) :: gij_1(nbins,kpol+1),gij_2(nbins,kpol+1)
   real(wp) :: L_k(nbins,kpol+1),L_k_1(nbins,kpol+1),L_k_2(nbins,kpol+1)
   real(wp) :: tab_gamma(nbins)
   integer  :: j,k


   ! call cpu_time(start)

   !SSPRK3 algo (Zhang & Shu 2010)
   !step 1
   call operator_L_kdv_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tabflux_coag,tabintflux_coag,dv,L_k)
   gij_1 = gij + dt*L_k


   !apply limiter coefficient to ensure positivity
   call gammafunction_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij_1,tab_gamma)
   do k=1,kpol
      gij_1(:,k+1) = tab_gamma(:)*gij_1(:,k+1)
   enddo


   !check gij values and limit to eps
   do j=1,nbins
      if (gij_1(j,1) < 0._wp) then
         print*,"error in calculating gij_1"
         print*,"j=",j,", gij_1 =",gij_1(j,1)
         stop

      else if ( gij_1(j,1) <= eps  ) then
         gij_1(j,1) = eps
         gij_1(j,2:) = 0._wp
      endif
   enddo



   !step 2
   call operator_L_kdv_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij_1,tabflux_coag,tabintflux_coag,dv,L_k_1)

   gij_2 = 3._wp*gij/4._wp + (gij_1 + dt*L_k_1)/4._wp

   !apply limiter coefficient to ensure positivity
   call gammafunction_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij_2,tab_gamma)
   do k=1,kpol
      gij_2(:,k+1) = tab_gamma(:)*gij_2(:,k+1)
   enddo

   !check gij values and limit to eps
   do j=1,nbins
      if (gij_2(j,1) < 0._wp) then
         print*,"error in calculating gij_2"
         print*,"j=",j,", gij_2 =",gij_2(j,1)
         stop

      else if ( gij_2(j,1) <= eps  ) then
         gij_2(j,1) = eps
         gij_2(j,2:) = 0._wp
      endif
   enddo


   !step 3
   call operator_L_kdv_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij_2,tabflux_coag,tabintflux_coag,dv,L_k_2)
   gijnew = gij/3._wp + 2._wp*(gij_2 + dt*L_k_2)/3._wp

   !apply limiter coefficient to ensure positivity
   call gammafunction_opt(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gijnew,tab_gamma)

   do k=1,kpol
      gijnew(:,k+1) = tab_gamma(:)*gijnew(:,k+1)
   enddo

   !check gij values and limit to eps
   do j=1,nbins
      if (gijnew(j,1) < 0._wp) then
         print*,"error in calculating gijnew"
         print*,"j=",j,", gijnew =",gijnew(j,1)
         stop

      else if ( gijnew(j,1) <= eps  ) then
         gijnew(j,1) = eps
         gijnew(j,2:) = 0._wp

      endif
   enddo


end subroutine time_solver_kdv_opt


end module solver_DG
