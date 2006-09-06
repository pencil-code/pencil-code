! $Id: planet.f90,v 1.67 2006-09-06 18:06:01 wlyra Exp $
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
  include 'planet.h' 
!  
! initialize variables needed for the planet module
!
!
! start and run parameters
!
!
! things needed for companion
!
  real :: gc=0. !planet's mass
  real :: b_pot=0.           !peak radius for potential
  integer :: nc=2        !exponent of smoothed potential 
  integer :: n_periods=5 !periods for ramping
  logical :: lramp=.false.,llocal_iso=.false.
  logical :: lmigrate=.false.!,lnorm=.false.
  real :: Gvalue=1. !gravity constant in same unit as density
  logical :: lcalc_turb=.false.
  integer :: nr=10
!
  namelist /planet_init_pars/ gc,nc,b_pot,&
       llocal_iso,lramp
!
  namelist /planet_run_pars/ gc,nc,b_pot,lramp, &
       llocal_iso,lmigrate,Gvalue,n_periods,nr,&
       lcalc_turb
! 
  integer :: idiag_torqint=0,idiag_torqext=0
  integer :: idiag_totenergy=0,idiag_totangmom=0
!
  contains
!
!***********************************************************************
    subroutine register_planet()
!
!  06-nov-05/wlad: coded
!
      use Mpicomm, only: stop_it
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_planet called twice')
      first = .false.
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_planet: ENTER'
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: planet.f90,v 1.67 2006-09-06 18:06:01 wlyra Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_planet: nvar > mvar')
      endif
!
    endsubroutine register_planet
!***********************************************************************
    subroutine initialize_planet(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  08-nov-05/wlad: coded
!
      use Mpicomm, only : stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      if (lroot) print*, 'initialize_planet'
!
    endsubroutine initialize_planet
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
    subroutine calc_pencils_planet(f,p)
!
! calculate keplerian velocity as a pencil, so I don't
! have to recalculate it on other routines from f,g0 and r0_pot
!
! 28-mar-06/wlad : coded 
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!   
    endsubroutine calc_pencils_planet
!***********************************************************************
    subroutine read_planet_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!            
      if (present(iostat)) then
        read(unit,NML=planet_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=planet_init_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_planet_init_pars
!***********************************************************************
    subroutine write_planet_init_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=planet_init_pars)
!                                                                         
    endsubroutine write_planet_init_pars
!***********************************************************************
    subroutine read_planet_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=planet_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=planet_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_planet_run_pars
!***********************************************************************
    subroutine write_planet_run_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=planet_run_pars)
!
    endsubroutine write_planet_run_pars
!***********************************************************************
    subroutine gravity_companion(rp_mn,rpcyl_mn,ax,ay,g0,r0_pot,n_pot,p)
!
!  calculate the gravity of a companion offcentered by (Rx,Ry,Rz)
!
!  12-aug-05/wlad       : coded
!  08-nov-05/wlad       : moved here (previously in Gravity Module)
!  21-nov-05/wlad+anders: changed orbital elements approach to
!                         particles' positions and velocities   
!
      use Sub
      use Global
      use Mpicomm
!     
      real, dimension (nx,nspar) :: rp_mn,rpcyl_mn
      real, dimension (nx,3) :: ggc,ggs
      real, dimension (nx) :: g_companion,g_star
      real, dimension (nx) :: rrs,rrs1,rrc,rrc1,rs,rp
      real, dimension (nspar) :: ax,ay
      real :: Omega_inertial,azp,gtc
      real :: axp,ayp,axs,ays,azs,g0,gp,gs,r0_pot
      integer :: n_pot
      logical :: lheader,lfirstcall=.true.
      type (pencil_case) :: p
!
!  Position of the particles
!  ipar(1) is planet, ipar(2) is star
!
      axp=ax(1); ayp=ay(1)
      axs=ax(2); ays=ax(2)
      azp=0.;azs=0
!
      lheader = lfirstcall .and. headtt .and. lroot
!
      if (headtt) print*,&
           'gravity_companion: Adding gravity of companion located at x,y,z=',&
           axp,ayp,azp
      if (lheader) print*,&
           'gravity_companion: Mass ratio of secondary-to-primary = ',gc/g0
!
      if (lheader.and.(Omega.eq.0)) print*,'gravity_companion: inertial frame'
!
      if (Omega/=0) &
           call stop_it('gravity_companion: Corotation not implemented')
!
!  Ramp up the mass of the planet for 5 periods
!
      ggc=0.
      if (gc.ne.0) then 
         call get_ramped_mass(gp,gs,g0)
!
!  Planet's gravity field - always spherical
!
         rrc =rp_mn(:,1) 
         rrc1=1./rp_mn(:,1)
!           
         g_companion=-gp*rrc*(rrc*rrc+b_pot*b_pot)**(-1.5)
!
         ggc(:,1) = (x(l1:l2)-axp)*rrc1*g_companion 
         ggc(:,2) = (y(  m  )-ayp)*rrc1*g_companion
         ggc(:,3) = (z(  n  )-azp)*rrc1*g_companion
!
      endif
!
!  Star's gravity field
!
      if (lcylindrical) then
         rrs=rpcyl_mn(:,2) 
      else
         rrs=rp_mn(:,2) 
      endif
!
      rrs1=1./rrs
!      
!  gravity_star will call the ramped mass from inside again
!
      ggs = 0.                  ! or ggs(:,3) might be uninitialized in
                                !  cylindrical case 
      call gravity_star(g0,r0_pot,n_pot,g_star,rrs)
!
      ggs(:,1) = (x(l1:l2)-axs)*rrs1*g_star
      ggs(:,2) = (y(  m  )-ays)*rrs1*g_star
      if (.not.lcylindrical) ggs(:,3) = (z(n)-azs)*rrs1*g_star
!
!  Reset gravity as global variable 
!      
      call set_global(ggs+ggc,m,n,'gg',nx)
!      
!  Stuff for calc_torque. Should maybe change it to particles_nbody
!
      rs = rpcyl_mn(:,2)
      rp = rpcyl_mn(:,1)
!
      if (ldiagnos) then
         if ((idiag_torqint/=0) .or. (idiag_torqext/=0)) &
              call calc_torque(p%rho,gp,axp,ayp,rpcyl_mn(:,1))
         
         if ((idiag_totenergy/=0).or.(idiag_totangmom/=0)) &
              call calc_monitored(rs,rp,axs,ays,axp,ayp,gs,gp,r0_pot,p)
      endif
      lfirstcall=.false.
!
    endsubroutine gravity_companion
!***********************************************************************
    subroutine calc_torque(dens,gp,ax,ay,rp)
!
! 05-nov-05/wlad : coded
!
      use Sub
!
      real, dimension(nx) :: torque,torqint,torqext
      real, dimension(nx) :: rp,rpre,dens
      real :: ax,ay,Rc,roche,gp
      integer :: i
!
! Planet's hills radius
!
      Rc = sqrt(ax*ax + ay*ay)
      roche = Rc*(gp/3.)**(1./3.)
!
      rpre = ax*y(m) - ay*x(l1:l2)
!
      torque = gp*dens*rpre*(rp*rp+b_pot*b_pot)**(-1.5)
!
      torqext=0.
      torqint=0.
!
      do i=1,nx
!
! Exclude material from inside the Roche Lobe
!
         if (rp(i).ge.roche) then
!
! External torque
!
            if (rcyl_mn(i).ge.Rc) torqext(i) = torque(i)
!
! Internal torque
!
            if (rcyl_mn(i).le.Rc) torqint(i) = torque(i)
!
         endif
      enddo
!
      call sum_lim_mn_name(torqext,idiag_torqext)
      call sum_lim_mn_name(torqint,idiag_torqint)
!
    endsubroutine calc_torque
!***********************************************************************
    subroutine get_ramped_mass(gp,gs,g0) 
!      
! Ramps up the mass of the planet from 0 to gc over
! n_periods orbits. If lramp=.false., will return gc.
! Currently just used for the comparison
! project. Called by both gravity_companion and dvvp_dt 
!
! 03-mar-06/wlad : coded
!
      real :: fgp,gp,gs,g0,tcut
      intent(out) :: gp,gs
!
      fgp = 0.5*(1 - sqrt(1-4*gc))
      gp = fgp
!
      if (lramp) then
         tcut = n_periods * 2*pi
         if (t .le. tcut) then
            gp = fgp* (sin(pi/2. * t/tcut))**2
         endif
      endif      
!
      gs = g0 - gp
!
    endsubroutine get_ramped_mass
!***********************************************************
   subroutine gravity_star(g0,r0_pot,n_pot,g_r,rr_mn)
!
! 08-nov-05/wlad : coded
!
     use Mpicomm, only : stop_it
!
     real, dimension (nx), intent(out) :: g_r
     real, dimension (nx) :: rr_mn
     integer :: i,n_pot
     real :: g0,r0_pot,gp,gs
!
     if (n_pot.ne.2) then
        print*,'planet: smoothed gravity used for star but smoothing lenght'
        print*,'is not equal to 2. Better stop and change, since it will lead'
        print*,'to boundary troubles as omega does not flatten properly.'
        call stop_it('')
     endif
!
     call get_ramped_mass(gp,gs,g0)
!
     g_r=-gs*rr_mn*(rr_mn*rr_mn+r0_pot*r0_pot)**(-1.5)
!   
   endsubroutine gravity_star
!***************************************************************
   subroutine calc_monitored(rs,rp,xs,ys,xp,yp,gs,gp,r0,p)

! calculate total energy and angular momentum
! and output their evolution as monitored variables
!
! 10-nov-05/wlad   : coded
!
     use Sub
! 
     real, dimension(nx) :: rs,rp,r,uphi
     real, dimension(nx) :: angular_momentum
     real, dimension(nx) :: kin_energy,pot_energy,total_energy
     real :: xs,ys,xp,yp  ! position of star and planet
     real :: gs,gp,r0     ! star's and planet's mass
     integer :: i
     type (pencil_case) :: p
!
     if (headtt) print*,'planet : calculate total energy'
     if (headtt) print*,'planet : calculate total angular momentum'
!
! Calculate the total energy and integrate it
!
! Kinetic energy
!
     kin_energy = p%rho * p%u2/2.
!     
! Potential energy - uses smoothed potential
!
     pot_energy = -1.*(gs*(rs*rs+r0*r0)**(-0.5) &
          + gp*(rp*rp+b_pot*b_pot)**(-0.5))*p%rho
!     
     total_energy = kin_energy + pot_energy  
!
     call sum_lim_mn_name(total_energy,idiag_totenergy)
!
! Angular momentum
!
     uphi = (-p%uu(:,1)*y(m) + p%uu(:,2)*x(l1:l2))/rcyl_mn
     angular_momentum = p%rho*rcyl_mn*uphi
!
     call sum_lim_mn_name(angular_momentum,idiag_totangmom)
!
   endsubroutine calc_monitored
!***************************************************************
    subroutine planet_phiavg(p)
!
      type (pencil_case) :: p
!
      if (lcalc_turb) then
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
    endsubroutine planet_phiavg
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
    subroutine rprint_planet(lreset,lwrite)
!
!  reads and registers monitoring quantities for planet-disk
!
!  06-nov-05/wlad: coded
!
      use Sub
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
         idiag_torqint=0
         idiag_torqext=0
         idiag_totenergy=0
         idiag_totangmom=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'rprint_gravity: run through parse list'
      do iname=1,nname
         call parse_name(iname,cname(iname),cform(iname),&
              'torqint',idiag_torqint)
         call parse_name(iname,cname(iname),cform(iname),&
              'torqext',idiag_torqext)
         call parse_name(iname,cname(iname),cform(iname),&
              'totenergy',idiag_totenergy)
         call parse_name(iname,cname(iname),cform(iname),&
              'totangmom',idiag_totangmom)
      enddo
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!  idl needs this even if everything is zero
!
      if (lwr) then
        write(3,*) 'i_torqint=',idiag_torqint
        write(3,*) 'i_torqext=',idiag_torqext
        write(3,*) 'i_totenergy=',idiag_totenergy
        write(3,*) 'i_totangmom=',idiag_totangmom
      endif
!
      if (NO_WARN) print*,lreset  !(to keep compiler quiet)
!
    endsubroutine rprint_planet
!***************************************************************
  endmodule Planet
  
