! $Id: planet.f90,v 1.71 2006-09-21 23:19:17 wlyra Exp $
!
!  This modules contains the routines for accretion disk and planet
!  building simulations. 
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lplanet = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
! This module takes care of (mostly) everything related to the 
! planet module
!
module Planet
!  
  use Cdata
  use Cparam
  use Messages
!  
  implicit none
!  
  public :: runtime_phiavg
!
  integer :: nr=10
! 
  contains
!
!***********************************************************************
    subroutine pencil_criteria_planet()
! 
!  All pencils that the Planet module depends on are specified here.
! 
!  06-nov-05/wlad: coded
!
      lpenc_requested(i_lnrho)=.true.
      lpenc_requested(i_rho)=.true.
      lpenc_requested(i_uu)=.true.
      lpenc_requested(i_u2)=.true.
      lpenc_requested(i_bb)=.true.
!
    endsubroutine pencil_criteria_planet
!***********************************************************************
    subroutine runtime_phiavg(p)
!
      type (pencil_case) :: p
!
      if (lfirst) then
!
! get the previous average to use in this timestep
!
         call get_old_average(p) 
!
! calculate the new average to use in next timestep
!
         call set_new_average(p)
!
      endif
!
    endsubroutine runtime_phiavg
!*******************************************************************
    subroutine get_old_average(p)
!
      use Sub, only: calc_phiavg_general,calc_phiavg_unitvects
      use Global, only: set_global,get_global
      use Mpicomm, only: stop_it
      use General, only: spline
!
      real, dimension(nx,3) :: bavg,uavg
      real, dimension(nr,3) :: bavg_coarse,uavg_coarse
      real, dimension(nx) :: rhoavg,ukepler
      real, dimension(nr) :: rhoavg_coarse,rcyl_coarse
      real :: rloop_int,rloop_ext,rmid,step
      integer :: i,ir,j
      type (pencil_case) :: p
      logical :: err
!
! Get the unit vectors
!
      call calc_phiavg_general()
      call calc_phiavg_unitvects()
!
! in the first time step, there is no average yet 
! use the pencil value then. The same will be used
! for r_int and r_ext 
!
      if (lmagnetic) then
         bavg(:,1)=p%bb(:,1)*pomx+p%bb(:,2)*pomy 
         bavg(:,2)=p%bb(:,1)*phix+p%bb(:,2)*phiy 
         bavg(:,3)=p%bb(:,3)                    
      endif
!
      uavg(:,1)=p%uu(:,1)*pomx+p%uu(:,2)*pomy
      uavg(:,2)=p%uu(:,1)*phix+p%uu(:,2)*phiy
      uavg(:,3)=p%uu(:,3)
!         
      if (it /= 1) then
!
! get the arrays of nr points calculated in set_new_average
! in the previous time-step
!
         call get_global(rhoavg_coarse,'rhoavg',nr)
         call get_global(uavg_coarse,'uavg',nr)
         if (lmagnetic) call get_global(bavg_coarse,'bavg',nr)
!
! expand it onto the pencil with spline interpolation
!
         step = (r_ext - r_int)/nr
         do ir=1,nr
            rloop_int = r_int + (ir-1)*step
            rloop_ext = r_int + ir*step
            rmid = 0.5*(rloop_int + rloop_ext)
            rcyl_coarse(ir)=rmid
         enddo   
!
         call spline(rcyl_coarse,rhoavg_coarse,rcyl_mn,rhoavg,nr,nx,err)
         do j=1,3 
            call spline(rcyl_coarse,uavg_coarse(:,j),rcyl_mn,uavg(:,j),nr,nx,err)
         enddo
         if (lmagnetic) then
            do j=1,3
               call spline(rcyl_coarse,bavg_coarse(:,j),rcyl_mn,bavg(:,j),nr,nx,err)
            enddo
         endif
!
! fill in with pencil values the parts of the array that are away from the interpolation 
!
         do i=1,nx
            if ((rcyl_mn(i).lt.rcyl_coarse(1)).or.(rcyl_mn(i).gt.rcyl_coarse(nr))) then
               rhoavg(i) = p%rho(i)
               uavg(i,1)=p%uu(i,1)*pomx(i)+p%uu(i,2)*pomy(i)
               uavg(i,2)=p%uu(i,1)*phix(i)+p%uu(i,2)*phiy(i)
               uavg(i,3)=p%uu(i,3)
               if (lmagnetic) then
                  bavg(i,1)=p%bb(i,1)*pomx(i)+p%bb(i,2)*pomy(i)
                  bavg(i,2)=p%bb(i,1)*phix(i)+p%bb(i,2)*phiy(i)
                  bavg(i,3)=p%bb(i,3)
               endif
            endif
         enddo
      endif
!
! Store the average onto a global variable, to be read on THIS
! timestep
!
      if (lmagnetic) call set_global(bavg,m,n,'bbs',nx)
      call set_global(uavg,m,n,'uus',nx)
      call set_global(rhoavg,m,n,'rhos',nx)
!
    endsubroutine get_old_average
!*******************************************************************
    subroutine set_new_average(p)
!
      use Sub, only: calc_phiavg_general,calc_phiavg_unitvects
      use Global, only: set_global
      use Mpicomm 
!
      real, dimension(nr,3) :: bavg_coarse,uavg_coarse
      real, dimension(nr) :: rhoavg_coarse,s_rho,rho_sum
      real, dimension(nr) :: s_uphi,s_urad,s_uzed
      real, dimension(nr) :: s_bphi,s_brad,s_bzed
      real, dimension(nr) :: up_sum,ur_sum,uz_sum
      real, dimension(nr) :: bp_sum,br_sum,bz_sum
      real, dimension(nx) :: uphi,urad,uzed
      real, dimension(nx) :: bphi,brad,bzed
      integer, dimension(nr) :: k,ktot
      real :: step,rloop_int,rloop_ext
      integer :: ir,i
      type (pencil_case) :: p
!
      call calc_phiavg_general()
      call calc_phiavg_unitvects()
!
      urad=p%uu(:,1)*pomx+p%uu(:,2)*pomy 
      uphi=p%uu(:,1)*phix+p%uu(:,2)*phiy 
      uzed=p%uu(:,3) 
      if (lmagnetic) then
         brad=p%bb(:,1)*pomx+p%bb(:,2)*pomy
         bphi=p%bb(:,1)*phix+p%bb(:,2)*phiy
         bzed=p%bb(:,3)
      endif
!
! number of radial zones
!
      step=(r_ext - r_int)/nr
!
! each zone has its limits rloop_int and rloop_ext
!
      s_uphi=0. ; s_urad=0. ; s_uzed=0.
      s_rho=0.
      if (lmagnetic) then 
         s_bphi=0. ; s_brad=0. ; s_bzed=0.
      endif
      k=0
!
      do ir=1,nr
         rloop_int = r_int + (ir-1)*step
         rloop_ext = r_int + ir*step
         do i=1,nx
            if ((rcyl_mn(i).le.rloop_ext).and.(rcyl_mn(i).ge.rloop_int)) then
!
               s_rho(ir)  = s_rho(ir) + p%rho(i)
!
               s_uphi(ir) = s_uphi(ir) + uphi(i)
               s_urad(ir) = s_urad(ir) + urad(i)
               s_uzed(ir) = s_uzed(ir) + uzed(i)
!               
               if (lmagnetic) then
                  s_bphi(ir) = s_bphi(ir) + bphi(i)
                  s_brad(ir) = s_brad(ir) + brad(i)
                  s_bzed(ir) = s_bzed(ir) + bzed(i)
               endif
!             
               k(ir)=k(ir)+1
!
            endif
         enddo
      enddo
!
! go filling the sums and the counter
!
      call mpireduce_sum(s_rho,rho_sum,nr)
!
      call mpireduce_sum(s_urad,ur_sum,nr)
      call mpireduce_sum(s_uphi,up_sum,nr)
      call mpireduce_sum(s_uzed,uz_sum,nr)
! 
      if (lmagnetic) then
         call mpireduce_sum(s_brad,br_sum,nr)
         call mpireduce_sum(s_bphi,bp_sum,nr)
         call mpireduce_sum(s_bzed,bz_sum,nr)
      endif
!
      call mpireduce_sum_int(k,ktot,nr)
!
! Broadcast the values
!
      call mpibcast_real(rho_sum,nr)
!
      call mpibcast_real(ur_sum,nr)
      call mpibcast_real(up_sum,nr)
      call mpibcast_real(uz_sum,nr)
!         
      if (lmagnetic) then
         call mpibcast_real(br_sum,nr)
         call mpibcast_real(bp_sum,nr)
         call mpibcast_real(bz_sum,nr)
      endif
!
      call mpibcast_int(ktot,nr)
!
! In the last point the sums are finalized. 
!
      if (llastpoint) then
!
! stop if any ktot is zero
!
         if (any(ktot == 0)) &
              call error("set_new_average","ktot=0") 

         rhoavg_coarse=rho_sum/ktot
!
         uavg_coarse(:,1)=ur_sum/ktot
         uavg_coarse(:,2)=up_sum/ktot
         uavg_coarse(:,3)=uz_sum/ktot
!
         if (lmagnetic) then
            bavg_coarse(:,1)=br_sum/ktot
            bavg_coarse(:,2)=bp_sum/ktot
            bavg_coarse(:,3)=bz_sum/ktot
         endif
!
! set the averages as global variables to use in the next timestep
!
         call set_global(rhoavg_coarse,'rhoavg',nr)
         call set_global(uavg_coarse,'uavg',nr)
         if (lmagnetic) call set_global(bavg_coarse,'bavg',nr)
!
       endif  
!       
     endsubroutine set_new_average
!***************************************************************
  endmodule Planet
  
