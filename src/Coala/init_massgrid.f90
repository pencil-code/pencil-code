!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


!--------------------------------------
! MODULE: initiate mass grid
!--------------------------------------
module init_massgrid 
   use precision
   
   implicit none

contains



!> @brief Generate massgrid and massbins in logscale from smin and smax
!!
!! @param[in]    nbins        number of dust bins
!! @param[in]    smax         maximum grain size
!! @param[in]    smin         minimum grain size
!! @param[in]    pi           pi value
!! @param[in]    rhograin     intrinsic grain density
!! @param[out]   sizegrid     grid of sizes, borders value of bins in size
!! @param[out]   sizemeanlog  geometric mean value of sizegrid for each bins
!! @param[out]   massgrid     grid of masses, borders value of mass bins
!! @param[out]   massbins     arithmetic mean value of massgrid for each mass bins
!! @param[out]   massmeanlog  geometric mean value of massgrid for each bins
subroutine compute_sizegrid_massgrid(nbins,smax,smin,pi,rhograin,sizegrid,sizemeanlog,massgrid,massbins,massmeanlog)
   implicit none
   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: smax,smin,pi,rhograin
   real(wp), intent(out) :: sizegrid(nbins+1),sizemeanlog(nbins)
   real(wp), intent(out) :: massgrid(nbins+1),massbins(nbins),massmeanlog(nbins)

   integer  :: j
   real(wp) :: r

   r=(smax/smin)**(1._wp/real(nbins,wp))
   
   sizegrid(1) = smin
   
   do j=2,nbins+1
      sizegrid(j) = r*sizegrid(j-1)
   enddo

   massgrid = 4._wp*pi*rhograin*sizegrid**3/3._wp
   
   do j=1,nbins
      sizemeanlog(j) = sqrt(sizegrid(j+1)*sizegrid(j))
      massbins(j)    = 0.5_wp*(massgrid(j+1)+massgrid(j))
      massmeanlog(j) = sqrt(massgrid(j+1)*massgrid(j))
   enddo
   print*,"rho minus plus: ",sizegrid(2),sizegrid(3),sizemeanlog(2)
      

end subroutine compute_sizegrid_massgrid


!> @brief Generate massgrid and massbins in logscale from massmin and massmax (for tests)
!!
!! @param[in]    nbins        number of dust bins
!! @param[in]    massmax      maximum grain mass
!! @param[in]    masssmin     minimum grain mass
!! @param[out]   massgrid     grid of masses, borders value of mass bins
!! @param[out]   massbins     arithmetic mean value of massgrid for each mass bins
!! @param[out]   massmeanlog  geometric mean value of massgrid for each mass bins
subroutine compute_massgrid(nbins,massmax,massmin,massgrid,massbins,massmeanlog)
   use precision
   implicit none
   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: massmax,massmin
   real(wp), intent(out) :: massgrid(nbins+1),massbins(nbins),massmeanlog(nbins)

   integer  :: j
   real(wp) :: r


   r=(massmax/massmin)**(1._wp/real(nbins,wp))
   
   !j==1
   massgrid(1) = massmin
   
   do j=2,nbins+1
      massgrid(j) = r*massgrid(j-1)
   enddo
   
   do j=1,nbins
      massbins(j) = 0.5_wp*(massgrid(j+1)+massgrid(j))
      massmeanlog(j) = sqrt(massgrid(j+1)*massgrid(j))
   enddo
      

end subroutine compute_massgrid


end module init_massgrid

