! $Id: planet.f90,v 1.57 2006-08-20 22:58:05 wlyra Exp $
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
  logical :: lcs2_global=.false.
  logical :: lmigrate=.false.!,lnorm=.false.
  real :: Gvalue=1. !gravity constant in same unit as density
  logical :: lcs2_thick=.false.
  logical :: lcalc_turb=.false.
  integer :: nr=10
!
  namelist /planet_init_pars/ gc,nc,b_pot,&
       lcs2_global,llocal_iso,lcs2_thick,lramp
!
  namelist /planet_run_pars/ gc,nc,b_pot,lramp, &
       llocal_iso,lcs2_global, &
       lmigrate,Gvalue,n_periods,lcs2_thick,nr,&
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
           "$Id: planet.f90,v 1.57 2006-08-20 22:58:05 wlyra Exp $")
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
!  02-feb-06/wlad: sound speed set as global in start time
!
!will need to add itorque as f-variable
!
      use Mpicomm, only : stop_it
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      logical :: lstarting
!
      if (lroot) print*, 'initialize_planet'
!
      if (lcs2_global) then
         if (llocal_iso) then
            do n=n1,n2
               do m=m1,m2
                  call set_soundspeed
               enddo
            enddo
         else
            print*,'planet : not isothermal, but sound speed is set'
            print*,'as global variable. Better stop and check'
            call stop_it('')
         endif
      endif
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
      real, dimension(mx,my,mz,mvar+maux) :: f
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
    subroutine gravity_companion(fp,dfp,rp_mn,rpcyl_mn,g0,r0_pot,n_pot,p)
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
      use Particles_Cdata
      use Mpicomm
!     
      real, dimension (mpar_loc,mpvar) :: fp,dfp 
      real, dimension (nx,mpar_loc) :: rp_mn,rpcyl_mn
      real, dimension (nx,3) :: ggc,ggs
      real, dimension (nx) :: g_companion,g_star
      real, dimension (nx) :: rrs,rrs1,rrc,rrc1
      real :: Omega_inertial,ax,ay,az,gtc
      real :: axs,ays,azs,g0,gp,gs,r0_pot,mdot,m2dot
      integer :: n_pot
      logical :: lheader,lfirstcall=.true.
      type (pencil_case) :: p
!
!  Position of the particles
!  1 is planet, 2 is star
!
      ax = fp(1,ixp) ; axs = fp(2,ixp) 
      ay = fp(1,iyp) ; ays = fp(2,iyp) 
      az = fp(1,izp) ; azs = fp(2,izp) 
!
      lheader = lfirstcall .and. headtt .and. lroot
!
      if (headtt) print*,&
           'gravity_companion: Adding gravity of companion located at x,y,z=',&
           ax,ay,az
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
         call get_ramped_mass(gp,gs,g0,mdot,m2dot)
!
!  Planet's gravity field - always spherical
!
         rrc =rp_mn(:,1) 
         rrc1=1./rp_mn(:,1)
!           
         g_companion=-gp*rrc*(rrc*rrc+b_pot*b_pot)**(-1.5)
!
         ggc(:,1) = (x(l1:l2)-ax)*rrc1*g_companion 
         ggc(:,2) = (y(  m  )-ay)*rrc1*g_companion
         ggc(:,3) = (z(  n  )-az)*rrc1*g_companion
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
      if (ldiagnos) then
         if ((idiag_torqint/=0) .or. (idiag_torqext/=0)) &
              call calc_torque(p%rho,gp,ax,ay,rpcyl_mn(:,1))
         
         if ((idiag_totenergy/=0).or.(idiag_totangmom/=0)) &
              call calc_monitored(rpcyl_mn,axs,ays,ax,ay,gs,gp,r0_pot,p)
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
    subroutine get_ramped_mass(gp,gs,g0,mdot,m2dot) 
!      
! Ramps up the mass of the planet from 0 to gc over
! n_periods orbits. If lramp=.false., will return gc.
! Currently just used for the comparison
! project. Called by both gravity_companion and dvvp_dt 
!
! 03-mar-06/wlad : coded
!
      real :: fgp,gp,gs,g0,tcut
      real :: mdot,m2dot
      intent(out) :: gp,gs,mdot,m2dot
!
      fgp = 0.5*(1 - sqrt(1-4*gc))
      gp = fgp
!
      mdot=0. ; m2dot=0.
!
      if (lramp) then
         tcut = n_periods * 2*pi
         if (t .le. tcut) then
            gp = fgp* (sin(pi/2. * t/tcut))**2
!            
            mdot  = 0.5*fgp* pi/tcut      * sin(pi*t/tcut)
            m2dot = 0.5*fgp*(pi/tcut)**2  * cos(pi*t/tcut)
!
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
     real :: g0,r0_pot,gp,gs,mdot,m2dot
!
     if (n_pot.ne.2) then
        print*,'planet: smoothed gravity used for star but smoothing lenght'
        print*,'is not equal to 2. Better stop and change, since it will lead'
        print*,'to boundary troubles as omega does not flatten properly.'
        call stop_it('')
     endif
!
     call get_ramped_mass(gp,gs,g0,mdot,m2dot)
!
     g_r=-gs*rr_mn*(rr_mn*rr_mn+r0_pot*r0_pot)**(-1.5)
!   
   endsubroutine gravity_star
!***************************************************************
   subroutine calc_monitored(rpcyl_mn,xs,ys,xp,yp,gs,gp,r0,p)

! calculate total energy and angular momentum
! and output their evolution as monitored variables
!
! 10-nov-05/wlad   : coded
!
     use Sub
! 
     real, dimension(nx,mpar_loc) :: rpcyl_mn
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
     rs = rpcyl_mn(:,2)
     rp = rpcyl_mn(:,1)
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
   subroutine set_soundspeed
!
! Locally isothermal structure for accretion disks. 
! The energy equation is not solved,but the variation
! of temperature with radius (T ~ r-1) is crucial for the
! treatment of ad hoc alpha viscosity.
!
! The obscure coefficients were calculated with an IDL
! routine to match a 1/r fall of cs2 with a constant value
! inside r = 0.4. Have to make it general some other day. 
!
! cs = H * Omega, being H the scale height and (H/r) = cte.
!
! This routine is called in start time only
!
! 25-nov-05/wlad : coded
! 
     use Global, only: set_global
!
      real, dimension (nx) :: rrp,cs2
      real :: g0=1.,r0_pot=0.2
      logical :: lheader,lfirstcall=.true.
!
!     It is better to set cs2 as a global variable than to recalculate it
!     at every timestep...
!
      lheader = headtt.and.lfirstcall.and.lroot
!
      if (lheader) print*,&
           'planet: setting sound speed as global variable'
!
      if (lcylindrical) then
         rrp = sqrt(x(l1:l2)**2 + y(m)**2) + tini
      else
         rrp = sqrt(x(l1:l2)**2 + y(m)**2 + z(n)**2) + tini
      endif   
!      
      if (lcs2_thick) then 
!
! This is for H=cs/Omega = Lz/2
!
         if (lheader) &
              print*,'set soundspeed: constant H=',Lxyz(3)/2.
!
         cs2 = g0**2 * (rrp**2 + r0_pot**2)**(-1.5) &
              *(Lxyz(3)/2.)**2
      else
!
! This one is for Mach Number = 20 everywhere 
!
         if (lheader) &
              print*,'set soundspeed: Mach=20 all over'
!
         where ((rrp.le.0.4).and.(rrp.ge.0.2)) 
         
            cs2 = 0.00749375        & 
                 +0.00679450*rrp    &
                 -0.0149310 *rrp**2 &
                 -0.0580586 *rrp**3 &
                 +0.0816006 *rrp**4 
         endwhere
!
         where (rrp.gt.0.4)
            cs2 = 0.05**2 * g0 / rrp
         endwhere
!
         where (rrp.lt.0.2) 
            cs2 = 0.008
         endwhere

      endif
!    
      call set_global(cs2,m,n,'cs2',nx)
!       
      lfirstcall=.false.
!
    endsubroutine set_soundspeed
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
!
      real, dimension(nx,3) :: bavg,uavg
      real, dimension(nr,3) :: bavg_coarse,uavg_coarse
      real, dimension(nx) :: rcyl,rhoavg
      real, dimension(nr) :: rhoavg_coarse
      real :: rloop_1,rloop_int,rloop_ext,rmid,rmid_1,rmid_2
      real :: dudr,dbdr,dr,step,upu,upb,dwr,u0,b0
      real :: upr,r0
      integer :: i,ir,aux
      type (pencil_case) :: p
!
! Get the unit vectors
!
      rcyl = sqrt(x(l1:l2)**2 + y(m)**2)
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
! expand it onto the pencil with linear interpolation
!
         step = (r_ext - r_int)/nr
!
         do i=1,nx
!
! this retrives ir, from 1 to nr
!
            ir = floor((rcyl(i)-r_int)/step) + 1 
            rloop_int = r_int + (ir-1)*step
            rloop_ext = r_int + ir*step
            rmid = 0.5*(rloop_int + rloop_ext)
            rmid_1 = rmid + step

            if ((rcyl(i) .gt. r_int).and.(rcyl(i) .le. r_ext)) then 
!
! gives problem for ir=nr because ir+1 will be out of bounds
!
               if (ir /= nr) then 
                  dwr = rmid_1 - rmid
                  dr = rcyl(i) - rmid
!
                  upr = rhoavg_coarse(ir+1) - rhoavg_coarse(ir)
                  r0 = rhoavg_coarse(ir)
                  rhoavg(i) = r0 + upr/dwr*dr
!
                  do aux=1,3
                     upu = uavg_coarse(ir+1,aux) - uavg_coarse(ir,aux)
                     dudr = upu/dwr
                     u0 = uavg_coarse(ir,aux)
                     uavg(i,aux) = u0 + dudr*dr
!                     
                     if (lmagnetic) then
                        upb = bavg_coarse(ir+1,aux) - bavg_coarse(ir,aux)
                        dbdr = upb/dwr
                        b0 = bavg_coarse(ir,aux)
                        bavg(i,aux) = b0 + dbdr*dr
                     endif   
                  enddo
!
               endif
!
! so do it backward for ir=nr
!
               if (ir == nr) then 
                  rmid_2 = rmid-step
                  dwr = rmid - rmid_2
                  dr = rcyl(i) - rmid
!
                  upr = rhoavg_coarse(ir) - rhoavg_coarse(ir-1)
                  r0 = rhoavg_coarse(ir)
                  rhoavg(i) = r0 + upr/dwr*dr
!
                  do aux=1,3
                     upu = uavg_coarse(ir,aux) - uavg_coarse(ir-1,aux)
                     dudr = upu/dwr
                     u0 = uavg_coarse(ir,aux)
                     uavg(i,aux) = u0 + dudr*dr
!
                     if (lmagnetic) then
                        upb = bavg_coarse(ir,aux) - bavg_coarse(ir-1,aux)
                        dbdr = upb/dwr
                        b0 = bavg_coarse(ir,aux)
                        bavg(i,aux) = b0 + dbdr*dr
                     endif
                  enddo
!
               endif
!            
            endif  !end if r_int r_ext
!
! close pencil loop
!
         enddo
!
! end (it.ne.1)
!
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
  
