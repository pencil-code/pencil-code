! $Id: planet.f90,v 1.75 2006-10-12 17:53:00 wlyra Exp $
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
  real, dimension(nrcylrun)  :: rho_tmp
  real, dimension(nrcylrun,3) :: u_tmp,b_tmp
  integer, dimension(nrcylrun) :: k_tmp
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
      lpenc_requested(i_rho)=.true.
      lpenc_requested(i_uu)=.true.
      lpenc_requested(i_bb)=.true.
!
    endsubroutine pencil_criteria_planet
!***********************************************************************
    subroutine runtime_phiavg(p)
!
      type (pencil_case) :: p
!
! get the previous average to use in this timestep
!
         call get_old_average(p) 
!
! calculate the new average to use in next timestep
!
         call set_new_average(p)
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
      real, dimension(nrcylrun,3) :: bavg_coarse,uavg_coarse
      real, dimension(nx) :: rhoavg,ukepler
      real, dimension(nrcylrun) :: rhoavg_coarse,rcyl_coarse
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
      if (lhydro) then
         uavg(:,1)=p%uu(:,1)*pomx+p%uu(:,2)*pomy
         uavg(:,2)=p%uu(:,1)*phix+p%uu(:,2)*phiy
         uavg(:,3)=p%uu(:,3)
      endif
!         
      if (it /= 1) then
!
! get the arrays of nrcylrun points calculated in set_new_average
! in the previous time-step
!
         if (ldensity)  call get_global(rhoavg_coarse,'rhoavg',nrcylrun)
         if (lhydro)    call get_global(uavg_coarse,'uavg',nrcylrun)
         if (lmagnetic) call get_global(bavg_coarse,'bavg',nrcylrun)
!
! expand it onto the pencil with spline interpolation
!
         step = (r_ext - r_int)/nrcylrun
         do ir=1,nrcylrun
            rloop_int = r_int + (ir-1)*step
            rloop_ext = r_int + ir*step
            rmid = 0.5*(rloop_int + rloop_ext)
            rcyl_coarse(ir)=rmid
         enddo   
!
         if (ldensity) &
              call spline(rcyl_coarse,rhoavg_coarse,rcyl_mn,rhoavg,nrcylrun,nx,err)
         do j=1,3 
            if (lhydro) &
                 call spline(rcyl_coarse,uavg_coarse(:,j),rcyl_mn,uavg(:,j),nrcylrun,nx,err)
            if (lmagnetic) &
                 call spline(rcyl_coarse,bavg_coarse(:,j),rcyl_mn,bavg(:,j),nrcylrun,nx,err)
         enddo      
!
! fill in with pencil values the parts of the array that are away from the interpolation 
!
         do i=1,nx
            if ((rcyl_mn(i).lt.rcyl_coarse(1)).or.(rcyl_mn(i).gt.rcyl_coarse(nrcylrun))) then
               if (ldensity) rhoavg(i) = p%rho(i)
               if (lhydro) then
                  uavg(i,1)=p%uu(i,1)*pomx(i)+p%uu(i,2)*pomy(i)
                  uavg(i,2)=p%uu(i,1)*phix(i)+p%uu(i,2)*phiy(i)
                  uavg(i,3)=p%uu(i,3)
               endif
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
      if (lhydro)    call set_global(uavg,m,n,'uus',nx)
      if (ldensity)  call set_global(rhoavg,m,n,'rhos',nx)
!
    endsubroutine get_old_average
!*******************************************************************
    subroutine set_new_average(p)
!
      use Sub, only: calc_phiavg_general,calc_phiavg_unitvects
      use Global, only: set_global
      use Mpicomm 
!
      real, dimension(nrcylrun,3) :: bavg_coarse,uavg_coarse
      real, dimension(nrcylrun,3) :: s_u,s_b,u_sum,b_sum
      real, dimension(nx,3) :: uuf,bbf
      real, dimension(nrcylrun) :: rhoavg_coarse,s_rho,rho_sum,ktot1
      integer, dimension(nrcylrun) :: k,ktot
      real :: step,rloop_int,rloop_ext
      real, dimension(nrcylrun) :: rmid
      integer :: ir,i,j
      type (pencil_case) :: p
!
      call calc_phiavg_general()
      call calc_phiavg_unitvects()
!
! rad, phi and zed
!
      if (lhydro) then
         uuf(:,1)=p%uu(:,1)*pomx+p%uu(:,2)*pomy 
         uuf(:,2)=p%uu(:,1)*phix+p%uu(:,2)*phiy 
         uuf(:,3)=p%uu(:,3) 
      endif
      if (lmagnetic) then
         bbf(:,1)=p%bb(:,1)*pomx+p%bb(:,2)*pomy
         bbf(:,2)=p%bb(:,1)*phix+p%bb(:,2)*phiy
         bbf(:,3)=p%bb(:,3)
      endif
!
! number of radial zones
!
      step=(r_ext - r_int)/nrcylrun
!
! each zone has its limits rloop_int and rloop_ext
!
      k=0
      if (ldensity)  s_rho=0.
      if (lhydro)    s_u=0.
      if (lmagnetic) s_b=0.
!
      do ir=1,nrcylrun
         rloop_int = r_int + (ir-1)*step
         rloop_ext = r_int + ir*step
         rmid(ir) = 0.5*(rloop_int + rloop_ext)
         do i=1,nx
            if ((rcyl_mn(i).le.rloop_ext).and.(rcyl_mn(i).ge.rloop_int)) then
!
               k(ir)=k(ir)+1
               if (ldensity)     s_rho(ir) = s_rho(ir) + p%rho(i)
               do j=1,3
                  if (lhydro)    s_u(ir,j) = s_u(ir,j) + uuf(i,j)
                  if (lmagnetic) s_b(ir,j) = s_b(ir,j) + bbf(i,j)
               enddo
!
            endif
         enddo
      enddo
!
! go filling the sums and the counter in this processor
!
      if (lfirstpoint) then
         k_tmp = k
         if (ldensity)  rho_tmp = s_rho
         if (lhydro)    u_tmp=s_u
         if (lmagnetic) b_tmp=s_b
      else
         k_tmp   = k_tmp + k
         if (ldensity)  rho_tmp = rho_tmp+s_rho
         if (lhydro)    u_tmp=u_tmp+s_u
         if (lmagnetic) b_tmp=b_tmp+s_b
      endif
!
! In the last point the sums are finalized
!
      if (llastpoint) then
!
! Sum across processors, send to root
!
         call mpireduce_sum_int(k_tmp,ktot,nrcylrun)
         if (ldensity) call mpireduce_sum(rho_tmp,rho_sum,nrcylrun)
         do j=1,3
            if (lhydro)    call mpireduce_sum(u_tmp(:,j),u_sum(:,j),nrcylrun)
            if (lmagnetic) call mpireduce_sum(b_tmp(:,j),b_sum(:,j),nrcylrun)
         enddo
!
! Broadcast the values
!
         call mpibcast_int(ktot,nrcylrun)
         if (ldensity) call mpibcast_real(rho_sum,nrcylrun)
         do j=1,3
            if (lhydro)    call mpibcast_real(u_sum(:,j),nrcylrun)
            if (lmagnetic) call mpibcast_real(b_sum(:,j),nrcylrun)
         enddo
!
! stop if any ktot is zero
!
         if (any(ktot == 0)) &
              call error("set_new_average","ktot=0") 
!
         ktot1=1./ktot
         if (ldensity)  rhoavg_coarse=rho_sum*ktot1
         do j=1,3 
            if (lhydro)    uavg_coarse(:,j)=u_sum(:,j)*ktot1
            if (lmagnetic) bavg_coarse(:,j)=b_sum(:,j)*ktot1
         enddo
!
! set the averages as global variables to use in the next timestep
!
         if (ldensity)  call set_global(rhoavg_coarse,'rhoavg',nrcylrun)
         if (lhydro)    call set_global(uavg_coarse,'uavg',nrcylrun)
         if (lmagnetic) call set_global(bavg_coarse,'bavg',nrcylrun)
!
      endif
!       
     endsubroutine set_new_average
!***************************************************************
  endmodule Planet
  
