!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


!-------------------------------------------
! MODULE: display progress bar
!-------------------------------------------

module progress_bar
   use precision
   
contains


!> @brief Display progress bar when computing arrays for flux and integral of flux and for time solver
!!
!! @param[in]    nsteps      number of step for the outer loop
!! @param[in]    iprogress   index of step
subroutine display_progress_bar(nsteps,iprogress)
   implicit none
   integer,   intent(in) :: iprogress,nsteps
   
   character :: bar(50)
   
   integer :: ind

   bar(1) = "["
   bar(50) = "]"
   bar(2:49) = " "

   !char(13) is the character for return

   if (nsteps==1) then
      bar(2:49) = "#"
      print*,bar,"100%"
   else 
      ind = iprogress*(50-1)/(nsteps)
      bar(2:ind) = "#"

      write(*,'(1a1,50a,3x,i3,a)',advance='no') char(13),bar,100*iprogress/(nsteps),"%"

      ! print*,bar
   endif


end subroutine display_progress_bar


!> @brief Display progress bar for time solver with subcycling coagulation time step
!!
!! @param[in]    ndthydro    number of hydro timestep
!! @param[in]    dtCFLsub    coagulation CFL
!! @param[in]    iprogress   index of step
subroutine display_progress_bar_subcycling(ndthydro,dtCFLsub,iprogress)
   implicit none
   integer,  intent(in) :: iprogress,ndthydro
   real(wp), intent(in) :: dtCFLsub
   
   character :: bar(50)
   
   integer :: ind

   bar(1) = "["
   bar(50) = "]"
   bar(2:49) = " "

   !char(13) is the character for return

   if (ndthydro==1) then
      bar(2:49) = "#"
      print*,bar,"100%"
   else 
      ind = iprogress*(50-1)/(ndthydro)
      bar(2:ind) = "#"

      write(*,'(1a1,50a,3x,i3,a,2x,1a11,e30.16E3)',advance='no') char(13),bar,100*iprogress/(ndthydro),"%","dtCFLsub = ",dtCFLsub

      ! print*,bar
   endif


end subroutine display_progress_bar_subcycling


end module progress_bar