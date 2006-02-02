! $Id: planet.f90,v 1.13 2006-02-02 12:26:53 ajohan Exp $
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
  
  use Cdata
  use Cparam
  use Messages
  
  implicit none
  
  include 'planet.h' 
  
  !initialize variables needed for the planet module

  integer :: dummy1
  namelist /planet_init_pars/ dummy1
    
  !
  ! run parameters
  !

  !things needed for companion
  !real :: Rx=0.,Ry=0.,Rz=0.
  real :: gc=0.    !location and mass
  real :: b=0.      !peak radius for potential
  integer :: nc=2   !exponent of smoothed potential 
  logical :: lramp=.false.
  logical :: lwavedamp=.false.,llocal_iso=.false.
  logical :: lsmoothlocal=.false.,lcs2_global=.false.
  logical :: lmigrate=.false.,lnorm=.false.
  real :: Gvalue=1. !gravity constant in same unit as density
!
  namelist /planet_run_pars/ gc,nc,b,lramp, &
       lwavedamp,llocal_iso,lsmoothlocal,lcs2_global, &
       lmigrate,lnorm,Gvalue
! 
  integer :: idiag_torqint=0,idiag_torqext=0
  integer :: idiag_torqrocheint=0,idiag_torqrocheext=0
  integer :: idiag_totalenergy=0,idiag_angularmomentum=0
  integer :: idiag_totmass2=0

  contains

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
           "$Id: planet.f90,v 1.13 2006-02-02 12:26:53 ajohan Exp $")
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
!will need to add itorque as f-variable
!
      use Mpicomm, only : stop_it
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      logical :: lstarting
!
      if (lroot) print*, 'initialize_planet'
      if ((lcs2_global).and.(.not.llocal_iso)) then
         print*,'planet : not isothermal, but sound speed is set'
         print*,'as global variable. Better stop and check'
         call stop_it('')
      endif
!
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
                                                                                                   
      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit
                                                                                                   
    endsubroutine read_planet_init_pars
!***********************************************************************
    subroutine write_planet_init_pars(unit)
      integer, intent(in) :: unit
                                                                                                   
      if (NO_WARN) print*,unit
                                                                                                   
    endsubroutine write_planet_init_pars
!***********************************************************************
    subroutine read_planet_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=planet_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=planet_run_pars,ERR=99)
      endif
                                                                                                                                                                                                
99    return
    endsubroutine read_planet_run_pars
!***********************************************************************
    subroutine write_planet_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=planet_run_pars)

    endsubroutine write_planet_run_pars
!***********************************************************************
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
         idiag_totalenergy=0
         idiag_angularmomentum=0
         idiag_totmass2=0
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
              'totalenergy',idiag_totalenergy)
         call parse_name(iname,cname(iname),cform(iname),&
              'angularmomentum',idiag_angularmomentum)
         call parse_name(iname,cname(iname),cform(iname),&
              'totmass2',idiag_totmass2)
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
        write(3,*) 'i_totalenergy=',idiag_totalenergy
        write(3,*) 'i_angularmomentum=',idiag_angularmomentum
        write(3,*) 'i_totmass2=',idiag_totmass2
      endif
!
      if (NO_WARN) print*,lreset  !(to keep compiler quiet)
!
    endsubroutine rprint_planet
!***********************************************************************
    subroutine gravity_companion(f,df,fp,dfp,g0,r0_pot,n_pot,p)
!
!  calculate the gravity of a companion offcentered by (Rx,Ry,Rz)
!
!  12-aug-05/wlad       : coded
!  08-nov-05/wlad       : moved here (previously in Gravity Module)
!  21-nov-05/wlad+anders: changed orbital elements approach to
!                         particles' positions and velocities   

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
      real, dimension (nx) :: g_companion,rrc,rrs,g_star,dens
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
      az = 0.        ; azs = 0. 
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
      gtc = gc
      if (lramp) then
          !ramping is for the comparison project
         if (t .le. 2*pi) then
            if (lheader) print*,&
                 'gravity_companion: Ramping up the mass of the companion'
            gtc = gc* (sin(t/4.))**2   !20 = pi/10*Period=2pi
         endif
      endif
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
      rrs=sqrt((x(l1:l2)-axs)**2+(y(m)-ays)**2+(z(n)-azs)**2) + tini
!      
      call gravity_star(g0,r0_pot,n_pot,g_star,axs,ays)
!
      ggs(:,1) = (x(l1:l2)-axs)/rrs*g_star
      ggs(:,2) = (y(  m  )-ays)/rrs*g_star
      ggs(:,3) = (z(  n  )-azs)/rrs*g_star     
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
        if ((idiag_torqint/=0) .or. (idiag_torqext/=0) .or. &
             (idiag_torqrocheint/=0) .or.(idiag_torqrocheext/=0)) then  
           dens = p%rho
           call calc_torque(dens,gtc,ax,ay,b)
        endif

        if ((idiag_totalenergy/=0).or.(idiag_angularmomentum/=0)) & 
             call calc_monitored(f,axs,ays,ax,ay,g0,gtc,r0_pot,p)

        lfirstcall = .false.
      endif
!
    endsubroutine gravity_companion
!***********************************************************************
    subroutine calc_torque(dens,gtc,ax,ay,b)

! 05-nov-05/wlad : coded

      use Sub
      use Cdata

      real, dimension(nx) :: r,torqint,torqext
      real, dimension(nx) :: torqrocheint,torqrocheext
      real, dimension(nx) :: gridmass
      real, dimension(nx) :: dist,rpre,dens
      real :: b,ax,ay,Rc,roche,gtc,r_lim
      integer :: i

      r_lim = 2.5

      r = sqrt(x(l1:l2)**2 + y(m)**2)

      dist = sqrt((x(l1:l2)-ax)**2 + (y(m)-ay)**2)
      rpre = ax*y(m) - ay*x(l1:l2)

      Rc = sqrt(ax**2 + ay**2)

      torqint = 0.   ; torqext=0.
      torqrocheint=0.; torqrocheext=0.; gridmass=0.
      
      roche = Rc*(gtc/3.)**(1./3.) !Jupiter roche

      !gridmass = dens*dx*dy/(pi*r_ext**2)

      do i=1,nx
         
         !external torque, excluding roche lobe
         if ((r(i).ge.Rc).and.(r(i).le.r_lim).and.(dist(i).ge.roche)) then 
            torqext(i) = gtc*dens(i)*rpre(i)*(dist(i)**2+b**2)**(-1.5)*dx*dy
         endif 
      
         !internal torque, excluding roche lobe
         if ((r(i).le.Rc).and.(r(i).ge.r_int).and.(dist(i).ge.roche)) then
            torqint(i) = gtc*dens(i)*rpre(i)*(dist(i)**2+b**2)**(-1.5)*dx*dy
         endif

         !external torque, roche lobe
         if ((r(i).ge.Rc).and.(dist(i).ge.0.5*roche).and.(dist(i).le.roche)) then 
            torqrocheext(i) = gtc*dens(i)*rpre(i)*(dist(i)**2+b**2)**(-1.5)*dx*dy
         endif 
      
         !internal torque, roche lobe
         if ((r(i).le.Rc).and.(dist(i).ge.0.5*roche).and.(dist(i).le.roche)) then 
            torqrocheint(i) = gtc*dens(i)*rpre(i)*(dist(i)**2+b**2)**(-1.5)*dx*dy
         endif


      enddo

      !call sum_mn_name(torque,idiag_torqint)

      call sum_mn_name(torqint,idiag_torqint)
      call sum_mn_name(torqext,idiag_torqext)
      call sum_mn_name(torqrocheint,idiag_torqrocheint)
      call sum_mn_name(torqrocheext,idiag_torqrocheext)
      
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

    endsubroutine local_isothermal
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

      use Cdata
      use Global

      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(mx,my,mz,mvar) :: df
      integer ider,j,k,i,ii
      real, dimension(nx,3) :: gg_mn
      real, dimension(nx) :: r,pdamp,aux0,velx0,vely0
      real :: tau,r_lim

      if (headtt) print*,&
           'wave_damping: damping motions for inner and outer boundary'

      r_lim=2.5
      
      tau = 2*pi/(0.4)**(-1.5)
      r = sqrt(x(l1:l2)**2 + y(m)**2)
     
      !for 0.4 : 1 ; 0.5 : 0
      pdamp = -11.111111*r**2 + 2.77777778  !parabolic function R
      !for 0.4 : 0 ; 0.5 : 1
      !pdamp = 11.1111111*r**2 - 1.777777778  
      
      where (r .le. 0.4) 
         pdamp = 1.
      endwhere
      where (r .ge. 0.5)
         pdamp = 0.
      endwhere
      
      call get_global(gg_mn,m,n,'gg')
      aux0 = 1./r*sqrt(gg_mn(:,1)**2+gg_mn(:,2)**2+gg_mn(:,3)**2)
      !aux0 is omega2
           
      !aux0 = g0*r**(n_pot-2)*(r**n_pot+r0_pot**n_pot)**(-1./n_pot-1.) 
      !aux is gravity, aux0 is omega**2
      
      velx0 =  -y(  m  ) * sqrt(aux0)   !initial conditions
      vely0 =  x(l1:l2) * sqrt(aux0)
      
      do i=l1,l2
         ii = i-l1+1
         if ((r(ii).le.0.5).and.(r(ii).gt.0.4)) then
            df(i,m,n,ilnrho) = df(i,m,n,ilnrho) - (f(i,m,n,ilnrho) - 1.       )/tau * pdamp(ii) 
            df(i,m,n,iux)    = df(i,m,n,iux)    - (f(i,m,n,iux)    - velx0(ii))/tau * pdamp(ii)
            df(i,m,n,iuy)    = df(i,m,n,iuy)    - (f(i,m,n,iuy)    - vely0(ii))/tau * pdamp(ii)
         endif
      enddo
      
     
     !outer boundary
     
     tau = 2*pi/(r_lim)**(-1.5)
     
     !!for 2.1 : 0 , 2.5 : 1
     pdamp = 0.543478*r**2 - 2.3967391  !parabolic function R
     !!for 2.1 : 1, 2.5 : 0
     !pdamp = -0.543478*r**2 + 3.3967375
     
     
     where (r .ge. 2.5) 
        pdamp = 1.
     endwhere
     where (r .le. 2.1)
        pdamp = 0.
     endwhere
     
     do i=l1,l2
        ii = i-l1+1
        if ((r(ii) .ge. 2.1).and.(r(ii).le.2.5)) then
           df(i,m,n,ilnrho) = df(i,m,n,ilnrho) - (f(i,m,n,ilnrho) - 1.       )/tau * pdamp(ii) 
           df(i,m,n,iux)    = df(i,m,n,iux)    - (f(i,m,n,iux)    - velx0(ii))/tau * pdamp(ii)
           df(i,m,n,iuy)    = df(i,m,n,iuy)    - (f(i,m,n,iuy)    - vely0(ii))/tau * pdamp(ii)
        endif
     enddo

   endsubroutine wave_damping
!***************************************************************
   subroutine gravity_star(g0,r0_pot,n_pot,g_r,xstar,ystar)
!
! 08-nov-05/wlad : coded
! 25-nov-05/wlad : call global cs2 here, if local isothermal
!
     use Cdata
     use Mpicomm, only : stop_it
!
     real, dimension (nx), intent(out) :: g_r
     real, dimension (nx) :: rr_mn
     real, optional :: xstar,ystar !initial position of star
     integer :: i,n_pot
     real :: g0,axs,ays,r0_pot
!
     if (present(xstar)) then;axs=xstar;else;axs=0.;endif
     if (present(ystar)) then;ays=ystar;else;ays=0.;endif
!  
     rr_mn = sqrt((x(l1:l2)-axs)**2 + (y(m)-ays)**2) + tini
!
     if (.not.lsmoothlocal) then 
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
     else

     !smooth the potential locally (only for r < r0)
     !
     !needed when the gravity profile 
     !should be like rr_mn in the inner boundary
     !but one does not want the ridiculous
     !smoothing that smooth_newton gives to n=2 

     !my gravity profile for r0 = 0.4
        where ((rr_mn.le.0.4).and.(rr_mn.ge.0.2))
           g_r = 2.01669                &
                -1.69235   *rr_mn**(-1) &
                -0.821300  *rr_mn**(-2) &
                +0.0812857 *rr_mn**(-3) &
                -0.00195754*rr_mn**(-4) 
        endwhere
        where (rr_mn.gt.0.4) 
           g_r = -g0/rr_mn**2
        endwhere
        where (rr_mn.lt.0.2) 
           g_r = -g0*rr_mn**(n_pot-1) &
                *(rr_mn**n_pot+r0_pot**n_pot)**(-1./n_pot-1.)
        endwhere
        
     endif

     if (lcs2_global.and.llocal_iso) then
        call set_soundspeed(g0,r0_pot,n_pot)
        if (llastpoint) lcs2_global=.false.
     endif

   endsubroutine gravity_star
!***************************************************************
   subroutine gravity_disk(fp,dfp,p,r0_pot)
!
! This routine takes care of the migration process, calculating
! the backreaction of the disk onto planet and star.
!
! 01-feb-06/wlad : coded
!
! ps: Apparently, there is a difference between a real number and 
! an array of dimension=1. As mpireduce deals with arrays, I
! have to define the sums as arrays as well. 
!
     use Mpicomm
     use Particles_cdata
!
     real, dimension (mpar_loc,mpvar) :: fp,dfp 
     real, dimension(nx) :: re,grav_gas
     real, dimension(1) :: sumx_loc,sumy_loc,sumx,sumy
     integer :: k
     real :: r_smooth,r0_pot
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

        re = sqrt((x(l1:l2) - fp(k,ixp))**2 +  (y(m) - fp(k,iyp))**2)
!
        grav_gas = Gvalue*p%rho*dx*dy*re/ &
             (re**2 + r_smooth**2)**(-1.5)
!                  
        sumx_loc(1) = sum(grav_gas * (x(l1:l2) - fp(k,ixp))/re)
        sumy_loc(1) = sum(grav_gas * (y(  m  ) - fp(k,iyp))/re)
!                  
        call mpireduce_sum(sumx_loc,sumx,1) 
        call mpireduce_sum(sumy_loc,sumy,1) 
        !call mpireduce_sum(sumz_loc,sumz,1)
!        
        if (lroot) then
           dfp(k,ivpx) = dfp(k,ivpx) + sumx(1)  
           dfp(k,ivpy) = dfp(k,ivpy) + sumy(1)
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
! 10-nov-05/wlad : coded


     use Sub
     use Cdata
 
     real, dimension(mx,my,mz,mvar+maux) :: f
     real, dimension(nx) :: rstar,rplanet,vel2,r,uphi
     real, dimension(nx) :: angular_momentum
     real, dimension(nx) :: kin_energy,pot_energy,total_energy
     real :: xs,ys,xp,yp !position of star and planet
     real :: g0,gp,r0_pot !star's and planet's mass
     real :: ang_tot,mass_tot,energy_tot
     integer :: i
     type (pencil_case) :: p

     if (headtt) print*,'planet : calculate total energy'
     if (headtt) print*,'planet : calculate total angular momentum'

     !calculate the total energy and integrate it

     rstar   = sqrt((x(l1:l2)-xs)**2 + (y(m)-ys)**2)+tini
     rplanet = sqrt((x(l1:l2)-xp)**2 + (y(m)-yp)**2)+tini

     !kinetic energy
     vel2 = f(l1:l2,m,n,iux)**2 + f(l1:l2,m,n,iuy)**2
     kin_energy = f(l1:l2,m,n,ilnrho) * vel2/2.
     
     !potential energy - uses smoothed potential
     pot_energy = -1.*(g0*(rstar**2+r0_pot**2)**(-1./2) &
          + gp*(rplanet**2+b**2)**(-1./2))*f(l1:l2,m,n,ilnrho)

     total_energy = kin_energy + pot_energy

     !integrate it
     energy_tot = -pi*(r_ext - r_int)
     call sum_lim_mn_name(total_energy,idiag_totalenergy, &
          energy_tot,lnorm)

     !angular momentum
     r = sqrt(x(l1:l2)**2 + y(m)**2)+tini !this is correct: baricenter
     uphi = (-f(l1:l2,m,n,iux)*y(m) + f(l1:l2,m,n,iuy)*x(l1:l2))/r
     angular_momentum = f(l1:l2,m,n,ilnrho) * r * uphi

     !integrate it
     ang_tot = 0.8*pi*(r_ext**2.5 - r_int**2.5)
     call sum_lim_mn_name(angular_momentum,idiag_angularmomentum, &
          ang_tot,lnorm)

     mass_tot= pi*(r_ext**2 - r_int**2)
     call sum_lim_mn_name(p%rho,idiag_totmass2, &
          mass_tot,lnorm)
          

   endsubroutine calc_monitored
!***************************************************************
   subroutine set_soundspeed(g0,r0_pot,n_pot)
!
! Locally isothermal structure for accretion disks. 
! The energy equation is not solved,but the variation
! of temperature with radius (T ~ r-1) is crucial for the
! treatment of ad hoc alpha viscosity.
!
! cs = H * Omega, being H the scale height and (H/r) = cte.
!
      use Cdata
      use Global, only: set_global
!
      real, dimension (nx) :: rrp,cs2
      real :: g0,r0_pot
      integer :: n_pot
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
      rrp = sqrt(x(l1:l2)**2 + y(m)**2)
!      
      where ((rrp.le.0.4).and.(rrp.ge.0.2)) 
         
         cs2 = 0.00749375        & 
              +0.00679450*rrp    &
              -0.0149310 *rrp**2 &
              -0.0580586 *rrp**3 &
              +0.0816006 *rrp**4 
      endwhere
      where (rrp.gt.0.4)
         !cs = 0.05 * Omega * r
         cs2 = 0.05**2 * g0 / rrp
      endwhere
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
  endmodule Planet
  
