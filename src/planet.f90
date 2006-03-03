! $Id: planet.f90,v 1.22 2006-03-03 15:06:28 wlyra Exp $
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
  real :: gc=0.          !location and mass
  real :: b=0.           !peak radius for potential
  integer :: nc=2        !exponent of smoothed potential 
  integer :: n_periods=5 !periods for ramping
  logical :: lramp=.false.
  logical :: lwavedamp=.false.,llocal_iso=.false.
  logical :: lsmoothlocal=.false.,lcs2_global=.false.
  logical :: lmigrate=.false.,lnorm=.false.
  real :: Gvalue=1. !gravity constant in same unit as density
!
  namelist /planet_init_pars/ gc,nc,b,lsmoothlocal,&
       lcs2_global,llocal_iso
!
  namelist /planet_run_pars/ gc,nc,b,lramp, &
       lwavedamp,llocal_iso,lsmoothlocal,lcs2_global, &
       lmigrate,lnorm,Gvalue,n_periods
! 
  integer :: idiag_torqint=0,idiag_torqext=0
  integer :: idiag_torqrocheint=0,idiag_torqrocheext=0
  integer :: idiag_totenergy=0,idiag_totangmom=0
  integer :: idiag_totmass=0
!
  contains
!
!***********************************************************************
    subroutine register_planet()
!
!  06-nov-05/wlad: coded
!
      use Cdata, only: ip,nvar,lroot !,mvar,nvar
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
           "$Id: planet.f90,v 1.22 2006-03-03 15:06:28 wlyra Exp $")
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
      use Cdata
!      
      lpenc_requested(i_lnrho)=.true.     
!
    endsubroutine pencil_criteria_planet
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
    subroutine gravity_companion(f,df,fp,dfp,g0,r0_pot,n_pot,p)
!
!  calculate the gravity of a companion offcentered by (Rx,Ry,Rz)
!
!  12-aug-05/wlad       : coded
!  08-nov-05/wlad       : moved here (previously in Gravity Module)
!  21-nov-05/wlad+anders: changed orbital elements approach to
!                         particles' positions and velocities   
!
      use Cdata
      use Sub
      use Global
      use Particles_Cdata
      use Mpicomm
!     
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp,dfp 
      real, dimension (nx,3) :: ggc,ggs
      real, dimension (nx) :: g_companion,rrc,rrs,g_star
      real :: Omega_inertial,ax,ay,az,gtc
      real :: axs,ays,azs,g0,r0_pot
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
      call get_ramped_mass(gtc)
!
!  Planet's gravity field
!
      rrc=sqrt((x(l1:l2)-ax)**2+(y(m)-ay)**2+(z(n)-az)**2) + tini
!           
      g_companion=-gtc*rrc**(nc-1) &
          *(rrc**nc+b**nc)**(-1./nc-1.)
!
      ggc(:,1) = (x(l1:l2)-ax)/rrc*g_companion 
      ggc(:,2) = (y(  m  )-ay)/rrc*g_companion
      ggc(:,3) = (z(  n  )-az)/rrc*g_companion
!
!  Star's gravity field
!
      if (lcylindrical) then
         rrs=sqrt((x(l1:l2)-axs)**2+(y(m)-ays)**2)+tini
      else
         rrs=sqrt((x(l1:l2)-axs)**2+(y(m)-ays)**2+(z(n)-azs)**2) + tini
      endif   
!      
      if (gc.ne.0) then
         call gravity_star(g0,r0_pot,n_pot,g_star,axs,ays,azs)
      else
         call gravity_star(g0,r0_pot,n_pot,g_star)
      endif
!
      ggs(:,1) = (x(l1:l2)-axs)/rrs*g_star
      ggs(:,2) = (y(  m  )-ays)/rrs*g_star
      ggs(:,3) = (z(  n  )-azs)/rrs*g_star     
      if (lcylindrical) ggs(:,3) = 0.
!
!  Reset gravity as global variable 
!  In the future, compute only potential and call grav=grad(pot)
!      
      call set_global(ggs+ggc,m,n,'gg',nx)
!
!  Add force on planet and star due to disk gravity
!
      if (lmigrate) call gravity_disk(fp,dfp,p,r0_pot)
!      
!  Stuff for calc_torque. Should maybe change it to particles_planet
!
      if (ldiagnos) then
         if ((idiag_torqint/=0) .or. (idiag_torqext/=0) &
              .or.(idiag_torqrocheint/=0) .or.(idiag_torqrocheext/=0)) &  
              call calc_torque(p%rho,gtc,ax,ay,b)
         
         if ((idiag_totenergy/=0).or.(idiag_totangmom/=0) &
              .or.(idiag_totmass/=0)) &
              call calc_monitored(f,axs,ays,ax,ay,g0,gtc,r0_pot,p)
      endif
      lfirstcall=.false.
!
    endsubroutine gravity_companion
!***********************************************************************
    subroutine calc_torque(dens,gtc,ax,ay,b)
!
! 05-nov-05/wlad : coded
!
      use Sub
      use Cdata
!
      real, dimension(nx) :: torqint,torqext,torque
      real, dimension(nx) :: torqint_rc,torqext_rc
      real, dimension(nx) :: r,re,rpre,dens
      real :: b,ax,ay,Rc,roche,gtc
      real :: sumtorqext,sumtorqint
      real :: sumtorqext_rc,sumtorqint_rc
      integer :: i
!
! Planet's hills radius
!
      Rc = sqrt(ax**2 + ay**2)
      roche = Rc*(gtc/3.)**(1./3.) 
!
      r  = sqrt(x(l1:l2)**2 + y(m)**2)
      re = sqrt((x(l1:l2)-ax)**2 + (y(m)-ay)**2)
      rpre = ax*y(m) - ay*x(l1:l2)
!
      torque = dens*rpre*(re**2+b**2)**(-1.5)
      torque = torque*gtc*dx*dy
!
      torqint = 0. ; torqint_rc = 0. 
      torqext = 0. ; torqext_rc = 0.
!
      print*,'dx,dy',dx,dy
      print*,'Rc,roche,gtc,b,ax,ay,sma',Rc,roche,gtc,b,ax,ay,sqrt(ax**2+ay**2)
!
      do i=1,nx
!         
! External torque
!
         if ((r(i).ge.Rc).and.(r(i).le.r_ext)) then
            if (re(i).ge.roche) then
               torqext(i)    = torque(i)
            else
               torqext_rc(i) = torque(i)
            endif
         endif
!  
! Internal torque
!
         if ((r(i).le.Rc).and.(r(i).ge.r_int)) then
            if (re(i).ge.roche) then
               torqint(i)    = torque(i)
            else
               torqint_rc(i) = torque(i)
            endif
         endif
!
      enddo
!
      sumtorqext    = sum(torqext)    ; sumtorqint    = sum(torqint)
      sumtorqext_rc = sum(torqext_rc) ; sumtorqint_rc = sum(torqint_rc)
!
      call surf_mn_name(sumtorqext,idiag_torqext)
      call surf_mn_name(sumtorqint,idiag_torqint)
      call surf_mn_name(sumtorqint_rc,idiag_torqrocheint)
      call surf_mn_name(sumtorqext_rc,idiag_torqrocheext)
!      
    endsubroutine calc_torque
!***********************************************************************
    subroutine local_isothermal(cs20,cs2)
!
!22-aug-05/wlad: coded
!08-nov-05/wlad: moved here (previously in the EoS module)
!25-nov-05/wlad: changed to be just a call of a global variable
!
      use Cdata
      use Global, only: get_global
!
      real, intent(in)  :: cs20
      real, dimension (nx), intent(out) :: cs2
!
      if (headtt) print*,&
           'planet: local isothermal equation of state for accretion disk'
!  
      call get_global(cs2,m,n,'cs2')
!
!
    endsubroutine local_isothermal
!***********************************************************
    subroutine get_ramped_mass(gtc)
!      
! Ramps up the mass of the planet from 0 to gc over
! n_period orbits. If lramp=.false., will return gc.
! Currently just used for the comparison
! project. Called by both gravity_companion and dvvp_dt 
!
! 03-mar-06/wlad : coded
!
      use Cdata
! 
      real :: gtc,tcut
!
      intent(out) :: gtc
!
      gtc = gc
      if (lramp) then
         tcut = n_period * 2*pi
         if (t .le. tcut) then
            gtc = gc* (sin(pi/2. * t/tcut)))**2   
         endif
      endif      
!
    endsubroutine get_ramped_mass
!***********************************************************
    subroutine wave_damping(f,df)
!
! 05-nov-05/wlad : coded
!
! Wave killing zone. Its only purpose is to have the same 
! specifications as Miguel's setup for the comparison project.
!
! As far as I understand, this thing is just like udamping 
! but with a varying pdamp, in the form of
! a parabola y=a*x**2 + c. Performs like freezing in the sense
! that it does not change the derivative in the boundary of the
! physical zone. But it does not stop the material completely
! in the outer boundary of the freezing ring.
!
!
      use Cdata
      use Global
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(mx,my,mz,mvar) :: df
      integer ider,j,k,i,ii
      real, dimension(nx,3) :: gg_mn
      real, dimension(nx) :: r,pdamp,aux0,velx0,vely0
      real :: tau
!
      if (headtt) print*,&
           'wave_damping: damping motions for inner and outer boundary'
!
      tau = 2*pi/(0.4)**(-1.5)
      r = sqrt(x(l1:l2)**2 + y(m)**2)
!     
! for 0.4 : 1 ; 0.5 : 0
!
      pdamp = -11.111111*r**2 + 2.77777778  !parabolic function R
!
! for 0.4 : 0 ; 0.5 : 1
! pdamp = 11.1111111*r**2 - 1.777777778  
!      
      where (r .le. 0.4) 
         pdamp = 1.
      endwhere
      where (r .ge. 0.5)
         pdamp = 0.
      endwhere
!      
      call get_global(gg_mn,m,n,'gg')
      aux0 = 1./r*sqrt(gg_mn(:,1)**2+gg_mn(:,2)**2+gg_mn(:,3)**2)
!
! aux0 is omega2
!           
! aux0 = g0*r**(n_pot-2)*(r**n_pot+r0_pot**n_pot)**(-1./n_pot-1.) 
! aux is gravity, aux0 is omega**2
!      
      velx0 =  -y(  m  ) * sqrt(aux0)   !initial conditions
      vely0 =  x(l1:l2) * sqrt(aux0)
!      
      do i=l1,l2
         ii = i-l1+1
         if ((r(ii).le.0.5).and.(r(ii).gt.0.4)) then
            df(i,m,n,ilnrho) = df(i,m,n,ilnrho) - (f(i,m,n,ilnrho) - 1.       )/tau * pdamp(ii) 
            df(i,m,n,iux)    = df(i,m,n,iux)    - (f(i,m,n,iux)    - velx0(ii))/tau * pdamp(ii)
            df(i,m,n,iuy)    = df(i,m,n,iuy)    - (f(i,m,n,iuy)    - vely0(ii))/tau * pdamp(ii)
         endif
      enddo
!     
! Outer boundary
!     
     tau = 2*pi/(2.5)**(-1.5)
!     
! for 2.1 : 0 , 2.5 : 1
!
     pdamp = 0.543478*r**2 - 2.3967391  !parabolic function R
!
! for 2.1 : 1, 2.5 : 0
! pdamp = -0.543478*r**2 + 3.3967375
!     
     where (r .ge. 2.5) 
        pdamp = 1.
     endwhere
     where (r .le. 2.1)
        pdamp = 0.
     endwhere
!     
     do i=l1,l2
        ii = i-l1+1
        if ((r(ii) .ge. 2.1).and.(r(ii).le.2.5)) then
           df(i,m,n,ilnrho) = df(i,m,n,ilnrho) - (f(i,m,n,ilnrho) - 1.       )/tau * pdamp(ii) 
           df(i,m,n,iux)    = df(i,m,n,iux)    - (f(i,m,n,iux)    - velx0(ii))/tau * pdamp(ii)
           df(i,m,n,iuy)    = df(i,m,n,iuy)    - (f(i,m,n,iuy)    - vely0(ii))/tau * pdamp(ii)
        endif
     enddo
!
   endsubroutine wave_damping
!***************************************************************
   subroutine gravity_star(g0,r0_pot,n_pot,g_r,xstar,ystar,zstar)
!
! 08-nov-05/wlad : coded
!
     use Cdata
     use Mpicomm, only : stop_it
!
     real, dimension (nx), intent(out) :: g_r
     real, dimension (nx) :: rr_mn
     real, optional :: xstar,ystar,zstar !initial position of star
     integer :: i,n_pot
     real :: g0,axs,ays,azs,r0_pot
!
     if (present(xstar)) then;axs=xstar;else;axs=0.;endif
     if (present(ystar)) then;ays=ystar;else;ays=0.;endif
     if (present(zstar)) then;azs=zstar;else;azs=0.;endif
!  
        if (lcylindrical) then
           rr_mn = sqrt((x(l1:l2)-axs)**2 + (y(m)-ays)**2) + tini
        else
           rr_mn = sqrt((x(l1:l2)-axs)**2 + (y(m)-ays)**2 + (z(n)-azs)**2) + tini
        endif
!
     if (n_pot.ne.2) then
        print*,'planet: smoothed gravity used for star but smoothing lenght'
        print*,'is not equal to 2. Better stop and change, since it will lead'
        print*,'to boundary troubles as omega does not flatten properly.'
        call stop_it('')
     endif
!
     g_r=-g0*rr_mn**(n_pot-1) &
          *(rr_mn**n_pot+r0_pot**n_pot)**(-1./n_pot-1.)
!   
   endsubroutine gravity_star
!***************************************************************
   subroutine gravity_disk(fp,dfp,p,r0_pot)
!
! This routine takes care of the migration process, calculating
! the backreaction of the disk onto planet and star.
!
! 01-feb-06/wlad+anders: coded
!
! WL: Apparently, there is a difference between a real number and 
! an array of dimension=1. As mpireduce deals with arrays, I
! have to define the sums as arrays as well. 
!
     use Mpicomm
     use Particles_cdata
!
     real, dimension (mpar_loc,mpvar) :: fp,dfp 
     real, dimension(nx) :: re,grav_gas
     real, dimension(1) :: sumx_loc,sumy_loc,sumz_loc
     real, dimension(1) :: sumx,sumy,sumz
     integer :: k
     real :: r_smooth,r0_pot,dv
     type (pencil_case) :: p
!
     do k=1,npar
        if (k==2) then
           r_smooth = r0_pot !star
        else if (k==1) then
           r_smooth = b      !planet
        else 
           call stop_it('gravity_disk - more than 2 planet particles?')
        endif   
!
        re = sqrt((x(l1:l2) - fp(k,ixp))**2 &
             +    (y(  m  ) - fp(k,iyp))**2 &
             +    (z(  n  ) - fp(k,izp))**2)
!       
        dv = dx*dy
        if (nzgrid/=1) dv = dv*dz
!
        grav_gas = Gvalue*p%rho*dv*re/ &
             (re**2 + r_smooth**2)**(-1.5)
!                  
        sumx_loc(1) = sum(grav_gas * (x(l1:l2) - fp(k,ixp))/re)
        sumy_loc(1) = sum(grav_gas * (y(  m  ) - fp(k,iyp))/re)
        sumz_loc(1) = sum(grav_gas * (z(  n  ) - fp(k,izp))/re)
!                  
        call mpireduce_sum(sumx_loc,sumx,1) 
        call mpireduce_sum(sumy_loc,sumy,1) 
        call mpireduce_sum(sumz_loc,sumz,1)
!        
        if (lroot) then
           dfp(k,ivpx) = dfp(k,ivpx) + sumx(1)  
           dfp(k,ivpy) = dfp(k,ivpy) + sumy(1)
           dfp(k,ivpz) = dfp(k,ivpz) + sumz(1)
        endif
!       
     enddo               
!
   endsubroutine gravity_disk
!***************************************************************
   subroutine calc_monitored(f,xs,ys,xp,yp,g0,gp,r0_pot,p)

! calculate total energy and angular momentum
! and output their evolution as monitored variables
!
! 10-nov-05/wlad   : coded
!
     use Sub
     use Cdata
! 
     real, dimension(mx,my,mz,mvar+maux) :: f
     real, dimension(nx) :: rstar,rplanet,vel2,r,uphi
     real, dimension(nx) :: angular_momentum
     real, dimension(nx) :: kin_energy,pot_energy,total_energy
     real :: xs,ys,xp,yp  !position of star and planet
     real :: g0,gp,r0_pot !star's and planet's mass
     real :: ang_tot,mass_tot,energy_tot
     integer :: i
     type (pencil_case) :: p
!
     if (headtt) print*,'planet : calculate total energy'
     if (headtt) print*,'planet : calculate total angular momentum'
!
! Calculate the total energy and integrate it
!
     rstar   = sqrt((x(l1:l2)-xs)**2 + (y(m)-ys)**2)+tini
     rplanet = sqrt((x(l1:l2)-xp)**2 + (y(m)-yp)**2)+tini
!
! Kinetic energy
!
     vel2 = f(l1:l2,m,n,iux)**2 + f(l1:l2,m,n,iuy)**2
     kin_energy = p%rho * vel2/2.
!     
! Potential energy - uses smoothed potential
!
     pot_energy = -1.*(g0*(rstar**2+r0_pot**2)**(-1./2) &
          + gp*(rplanet**2+b**2)**(-1./2))*p%rho
!     
     total_energy = kin_energy + pot_energy  
!     
! integrate it
!
     energy_tot = -pi*(r_ext - r_int)
     call sum_lim_mn_name(total_energy,idiag_totenergy, &
          energy_tot,lnorm)
!
! Angular momentum
!
     r = sqrt(x(l1:l2)**2 + y(m)**2)+tini !this is correct: baricenter
     uphi = (-f(l1:l2,m,n,iux)*y(m) + f(l1:l2,m,n,iuy)*x(l1:l2))/r
     angular_momentum = p%rho * r * uphi
!     
! integrate it
!
     ang_tot = 0.8*pi*(r_ext**2.5 - r_int**2.5)
     call sum_lim_mn_name(angular_momentum,idiag_totangmom, &
          ang_tot,lnorm)
!
! Total mass
!
     mass_tot= pi*(r_ext**2 - r_int**2)
     call sum_lim_mn_name(p%rho,idiag_totmass, &
          mass_tot,lnorm)
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
! 25-nov-05/wlad : coded
! 
     use Cdata
      use Global, only: set_global
!
      real, dimension (nx) :: rrp,cs2
      real :: g0=1.
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
!    
      call set_global(cs2,m,n,'cs2',nx)
!       
      lfirstcall=.false.
!
    endsubroutine set_soundspeed
!***************************************************************
    subroutine rprint_planet(lreset,lwrite)
!
!  reads and registers monitoring quantities for planet-disk
!
!  06-nov-05/wlad: coded
!
      use Cdata
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
         idiag_torqrocheint=0
         idiag_torqrocheext=0
         idiag_totenergy=0
         idiag_totangmom=0
         idiag_totmass=0
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
              'torqrocheint',idiag_torqrocheint)
         call parse_name(iname,cname(iname),cform(iname),&
              'torqrocheext',idiag_torqrocheext)
         call parse_name(iname,cname(iname),cform(iname),&
              'totenergy',idiag_totenergy)
         call parse_name(iname,cname(iname),cform(iname),&
              'totangmom',idiag_totangmom)
         call parse_name(iname,cname(iname),cform(iname),&
              'totmass',idiag_totmass)
      enddo
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!  idl needs this even if everything is zero
!
      if (lwr) then
        write(3,*) 'i_torqint=',idiag_torqint
        write(3,*) 'i_torqext=',idiag_torqext
        write(3,*) 'i_torqrocheint=',idiag_torqrocheint
        write(3,*) 'i_torqrocheext=',idiag_torqrocheext
        write(3,*) 'i_totenergy=',idiag_totenergy
        write(3,*) 'i_totangmom=',idiag_totangmom
        write(3,*) 'i_totmass=',idiag_totmass
      endif
!
      if (NO_WARN) print*,lreset  !(to keep compiler quiet)
!
    endsubroutine rprint_planet
!***************************************************************
  endmodule Planet
  
