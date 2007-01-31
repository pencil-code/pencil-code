! $Id: planet.f90,v 1.80 2007-01-31 12:21:16 wlyra Exp $
!
!  This modules contains the routines for accretion disc and planet
!  building simulations.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lplanet = .true.
!
! PENCILS PROVIDED uavg,bavg,rhoavg
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
  include 'planet.h'
!
  real, dimension(nrcyl)  :: rho_tmp
  real, dimension(nrcyl,3) :: u_tmp,b_tmp
  integer, dimension(nrcyl) :: k_tmp
!
  real, dimension (nrcyl,3) :: bavg_coarse,uavg_coarse
  real, dimension (nrcyl) :: rhoavg_coarse
!
  contains
!
!***********************************************************************
    subroutine pencil_criteria_planet()
!
!  All pencils that the Planet module depends on are specified here.
!
!  06-nov-05/wlad: coded
!  16-nov-06/tony: pencilised coordinates and averages
!
      lpenc_requested(i_rho)=.true.
      lpenc_requested(i_uu)=.true.
      lpenc_requested(i_bb)=.true.
      lpenc_requested(i_rcyl_mn)=.true.
      if (lmagnetic.or.lhydro) then
        lpenc_requested(i_phix)=.true.
        lpenc_requested(i_phiy)=.true.
        lpenc_requested(i_pomx)=.true.
        lpenc_requested(i_pomy)=.true.
      endif
      if (lmagnetic) lpenc_requested(i_bavg)=.true.
      if (lhydro)    lpenc_requested(i_uavg)=.true.
      if (ldensity)  lpenc_requested(i_rhoavg)=.true.
!
    endsubroutine pencil_criteria_planet
!***********************************************************************
    subroutine pencil_interdep_planet(lpencil_in)
!
!  Interdependency among pencils from the Planet module is specified here.
!
!  16-nov-06/tony: coded
!
      logical, dimension(npencils) :: lpencil_in
!
!      if (lpencil_in(i_u2)) lpencil_in(i_uu)=.true.
!
    endsubroutine pencil_interdep_planet
!*******************************************************************
    subroutine calc_pencils_planet(f,p)
!
!  DOCUMENT ME
!
!  ??-???-??/wlad: coded
!
      use Mpicomm
      use General, only: spline
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension(nx) :: ukepler
      real, dimension(nrcyl) :: rcyl_coarse
      real, dimension(nrcyl) :: rmiddle
      real, dimension(nrcyl) :: s_rho,rho_sum,ktot1
      real, dimension(nrcyl,3) :: s_u,s_b,u_sum,b_sum
      real, dimension(nx,3) :: uuf,bbf
      real :: rloop_int,rloop_ext,rmid,step
      integer, dimension(nrcyl) :: k,ktot
      integer :: i,j,ir
      logical :: err
!
! in the first time step, there is no average yet
! use the pencil value then. The same will be used
! for r_int and r_ext
!
      if (lmagnetic) then
        p%bavg(:,1)=p%bb(:,1)*p%pomx+p%bb(:,2)*p%pomy
        p%bavg(:,2)=p%bb(:,1)*p%phix+p%bb(:,2)*p%phiy
        p%bavg(:,3)=p%bb(:,3)
      endif
!
      if (lhydro) then
        p%uavg(:,1)=p%uu(:,1)*p%pomx+p%uu(:,2)*p%pomy
        p%uavg(:,2)=p%uu(:,1)*p%phix+p%uu(:,2)*p%phiy
        p%uavg(:,3)=p%uu(:,3)
      endif
!
      if (.not.headtt) then
!
! expand it onto the pencil with spline interpolation
!
        step = (r_ext - r_int)/nrcyl
        do ir=1,nrcyl
          rloop_int = r_int + (ir-1)*step
          rloop_ext = r_int + ir*step
          rmid = 0.5*(rloop_int + rloop_ext)
          rcyl_coarse(ir)=rmid
        enddo
!
        if (ldensity) &
          call spline(rcyl_coarse,rhoavg_coarse,p%rcyl_mn,p%rhoavg,nrcyl,nx,err)
        do j=1,3
          if (lhydro) &
            call spline(rcyl_coarse,uavg_coarse(:,j),p%rcyl_mn,p%uavg(:,j),nrcyl,nx,err)
          if (lmagnetic) &
            call spline(rcyl_coarse,bavg_coarse(:,j),p%rcyl_mn,p%bavg(:,j),nrcyl,nx,err)
        enddo
!
! fill in with pencil values the parts of the array that are away from the interpolation
!
        do i=1,nx
          if ((p%rcyl_mn(i).lt.rcyl_coarse(1)).or.(p%rcyl_mn(i).gt.rcyl_coarse(nrcyl))) then
            if (ldensity) p%rhoavg(i) = p%rho(i)
            if (lhydro) then
              p%uavg(i,1)=p%uu(i,1)*p%pomx(i)+p%uu(i,2)*p%pomy(i)
              p%uavg(i,2)=p%uu(i,1)*p%phix(i)+p%uu(i,2)*p%phiy(i)
              p%uavg(i,3)=p%uu(i,3)
            endif
            if (lmagnetic) then
              p%bavg(i,1)=p%bb(i,1)*p%pomx(i)+p%bb(i,2)*p%pomy(i)
              p%bavg(i,2)=p%bb(i,1)*p%phix(i)+p%bb(i,2)*p%phiy(i)
              p%bavg(i,3)=p%bb(i,3)
            endif
          endif
        enddo
      endif
!
! rad, phi and zed
!
      if (lhydro) then
        uuf(:,1)=p%uu(:,1)*p%pomx+p%uu(:,2)*p%pomy
        uuf(:,2)=p%uu(:,1)*p%phix+p%uu(:,2)*p%phiy
        uuf(:,3)=p%uu(:,3)
      endif
      if (lmagnetic) then
        bbf(:,1)=p%bb(:,1)*p%pomx+p%bb(:,2)*p%pomy
        bbf(:,2)=p%bb(:,1)*p%phix+p%bb(:,2)*p%phiy
        bbf(:,3)=p%bb(:,3)
      endif
!
! number of radial zones
!
      step=(r_ext - r_int)/nrcyl
!
! each zone has its limits rloop_int and rloop_ext
!
      k=0
      if (ldensity)  s_rho=0.
      if (lhydro)    s_u=0.
      if (lmagnetic) s_b=0.
!
      do ir=1,nrcyl
        rloop_int = r_int + (ir-1)*step
        rloop_ext = r_int + ir*step
        rmiddle(ir) = 0.5*(rloop_int + rloop_ext)
        do i=1,nx
          if ((p%rcyl_mn(i).le.rloop_ext).and.(p%rcyl_mn(i).ge.rloop_int)) then
!
            k(ir)=k(ir)+1
            if (ldensity)    s_rho(ir) = s_rho(ir) + p%rho(i)
!
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
        call mpireduce_sum_int(k_tmp,ktot,nrcyl)
        if (ldensity) call mpireduce_sum(rho_tmp,rho_sum,nrcyl)
        do j=1,3
          if (lhydro)    call mpireduce_sum(u_tmp(:,j),u_sum(:,j),nrcyl)
          if (lmagnetic) call mpireduce_sum(b_tmp(:,j),b_sum(:,j),nrcyl)
        enddo
!
! Broadcast the values
!
        call mpibcast_int(ktot,nrcyl)
        if (ldensity) call mpibcast_real(rho_sum,nrcyl)
        do j=1,3
          if (lhydro)    call mpibcast_real(u_sum(:,j),nrcyl)
          if (lmagnetic) call mpibcast_real(b_sum(:,j),nrcyl)
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
      endif
!
    endsubroutine calc_pencils_planet
!***************************************************************
  endmodule Planet

