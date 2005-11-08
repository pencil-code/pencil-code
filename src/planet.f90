! $Id: planet.f90,v 1.2 2005-11-08 23:42:05 wlyra Exp $
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
  real :: Rx=0.,Ry=0.,Rz=0.,gc=0.  !location and mass
  real :: b=0.      !peak radius for potential
  integer :: nc=2   !exponent of smoothed potential 
  logical :: lcompanion=.false.,lramp=.false.,llocal_iso=.false.
  logical :: lwavedamp=.false.

  namelist /planet_run_pars/ Rx,Ry,Rz,gc,lcompanion,nc,b,lramp, &
       llocal_iso,lwavedamp
! 

  integer :: idiag_torqint=0,idiag_torqext=0
  integer :: idiag_torqrocheint=0,idiag_torqrocheext=0
  integer :: idiag_totalmass=0

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
           "$Id: planet.f90,v 1.2 2005-11-08 23:42:05 wlyra Exp $")
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

      real, dimension (mx,my,mz,mvar+maux) :: f
      logical :: lstarting

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
         idiag_totalmass=0
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
              'totalmass',idiag_totalmass)

      enddo
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!  idl needs this even if everything is zero
!
      if (lwr) then
        write(3,*) 'i_torqint=',idiag_torqint
        write(3,*) 'i_torqext=',idiag_torqext
        write(3,*) 'i_torqrocheint',idiag_torqrocheint
        write(3,*) 'i_torqrocheext',idiag_torqrocheext
        write(3,*) 'i_totalmass',idiag_totalmass
        write(3,*) 'nname=',nname
      endif
!
      if (NO_WARN) print*,lreset  !(to keep compiler quiet)
!
    endsubroutine rprint_planet
!***********************************************************************
    subroutine gravity_companion(f,df,p)
!
!  add duu/dt according to the gravity of a companion offcentered by (Rx,Ry,Rz)
!
!  12-aug-05/wlad: coded
!  08-nov-05/wlad: moved here (previously in Gravity Module)
!

      use Cdata
      use Sub
      use Global
      use Gravity
!     
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: ggc,ggs
      real, dimension (nx) :: g_companion,rrc,rrs,g_star,dens
      real :: Omega_inertial,Rc,phase,phi,ax,ay,az,gtc
      real :: Rs, axs,ays,azs,mc,ms,phis

      gtc = gc
      if (lramp) then
         !ramp up the mass of the planet for 5 periods
         !ramping is for the comparison project
         if (t .le. 10*pi) then
            if (headtt) print*,&
                 'gravity_companion: Ramping up the mass of the companion'
            gtc = gc* (sin(t/20))**2  !20 = pi/10*Period=2pi
         endif
      endif

      if (headtt) print*,&
           'gravity_companion: Adding gravity of companion located at x,y,z=',Rx,Ry,Rz
      if (headtt) print*,&
           'gravity_companion: Mass ratio of secondary-to-primary = ',gc/g0

      if (Omega /= 0) then  
         if (headtt) print*,'gravity_companion: corotational frame'
         ax = Rx  
         ay = Ry            
      else
         if (headtt) print*,'gravity_companion: inertial frame'
         !add these three to grav_run_pars later
         Rc = sqrt(Rx**2+Ry**2+Rz**2)
         Omega_inertial = sqrt(g0/Rc**3)
         phase = acos(Rx/Rc)

         phi = Omega_inertial*t + phase  
         ax = Rc * cos(phi) 
         ay = Rc * sin(phi) 
      endif
            
      az = Rz !just for notational consistency

      rrc=sqrt((x(l1:l2)-ax)**2+(y(m)-ay)**2+(z(n)-az)**2)


      !where (rrc.ge.b)
      !   g_companion = -gtc/rrc**2
      !elsewhere
      !   g_companion = -gtc*(    &
      !        -11.25*rrc**5/b**7 &
      !        +35.00*rrc**4/b**6 &
      !        -31.50*rrc**3/b**5 &
      !        + 8.75*rrc   /b**3 )
      !endwhere
      
      
      g_companion=-gtc*rrc**(nc-1) &
          *(rrc**nc+b**nc)**(-1./nc-1.)

         ggc(:,1) = (x(l1:l2)-ax)/(rrc)*g_companion 
         ggc(:,2) = (y(  m  )-ay)/(rrc)*g_companion
         ggc(:,3) = (z(  n  )-az)/(rrc)*g_companion
      ! 

      df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + ggc
      
      if (lwavedamp) call wave_damping(f,df)

      !gravity star

      phis = phi + pi
      mc = gtc ; ms = 1 - gc ; Rs = Rc*(gc/ms) 
      
      axs = Rs * cos(phis) 
      ays = Rs * sin(phis) 
      azs = 0. !notational consistency

      rrs=sqrt((x(l1:l2)-axs)**2+(y(m)-ays)**2+(z(n)-azs)**2)
      
      g_star=-g0*rrs**(n_pot-1) &
          *(rrs**n_pot+r0_pot**n_pot)**(-1./n_pot-1.)

      ggs(:,1) = (x(l1:l2)-axs)/(rrs)*g_star
      ggs(:,2) = (y(  m  )-ays)/(rrs)*g_star
      ggs(:,3) = (z(  n  )-azs)/(rrs)*g_star
      
      !reset the gravity from the star as global variable 
      
      call set_global(ggs,m,n,'gg',nx)
      
      !stuff for calc_torque
      dens = f(l1:l2,m,n,ilnrho)
      call calc_torque(dens,gtc,ax,ay,b,Rc,phi)

    endsubroutine gravity_companion
!***********************************************************************
    subroutine calc_torque(dens,gtc,ax,ay,b,Rc,phi)

! 05-nov-05/wlad : coded

      use Sub
      use Cdata

      real, dimension(nx) :: r,torqint,torqext
      real, dimension(nx) :: torqrocheint,torqrocheext
      real, dimension(nx) :: dist,rpre,dens
      real :: b,ax,ay,Rc,phi,roche,gtc
      integer :: i

      r = sqrt(x(l1:l2)**2 + y(m)**2)

      dist = sqrt((x(l1:l2)-ax)**2 + (y(m)-ay)**2)
      rpre = ax*y(m) - ay*x(l1:l2)

      torqint = 0.   ; torqext=0.
      torqrocheint=0.; torqrocheext=0.


      roche = Rc*(gtc/3.)**(1./3.) !Jupiter roche

      do i=1,nx
         !external torque, excluding roche lobe
         if ((r(i).ge.Rc).and.(r(i).le.r_ext).and.(dist(i).ge.roche)) then 
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

      if (idiag_torqint/=0) call sum_mn_name(torqint,idiag_torqint)
      if (idiag_torqext/=0) call sum_mn_name(torqext,idiag_torqext)
      if (idiag_torqrocheint/=0) call sum_mn_name(torqrocheint,idiag_torqrocheint)
      if (idiag_torqrocheext/=0) call sum_mn_name(torqrocheext,idiag_torqrocheext)

    endsubroutine calc_torque
!***********************************************************************
!    subroutine local_isothermal(cs20,cs2)
!!
!!22-aug-05/wlad: coded
!!08-nov-05/wlad: moved here (previously in the EoS module)
!!
!! Locally isothermal structure for accretion disks. 
!! The energy equation is not solved,but the variation
!! of temperature with radius (T ~ r-1) is crucial for the
!! treatment of ad hoc alpha viscosity.
!!
!! cs = H * Omega, being H the scale height and (H/r) = cte.
!!
!      use Cdata
!      use Global, only: get_global
!
!      real, dimension(nx,3) :: gg_mn
!      real, dimension(nx) :: rr_mn,gr,Om
!      real :: Mach,plaw
!
!      real, intent(in)  :: cs20
!      real, dimension (nx), intent(out) :: cs2
!!
!
!      if (headtt) print*,&
!           'planet: local isothermal equation of state for accretion disk'
!      
!      plaw = 0.
!
!      Mach = sqrt(1./cs20)
!
!      rr_mn = sqrt(x(l1:l2)**2 + y(m)**2 + z(n)**2) + epsi
!
!      call get_global(gg_mn,m,n,'gg')
!      gr = sqrt(gg_mn(:,1)**2+gg_mn(:,2)**2+gg_mn(:,3)**2)
!      
!      Om = sqrt(gr/rr_mn * (1 + plaw/Mach**2)**(-1))
!      
!      !0.5 is plaw
!
!      cs2 = (Om * rr_mn / Mach)**2
!
!    endsubroutine local_isothermal
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
! physical zone, but it does not stop the material completely
! in the outer boundary of the freezing ring.
!

      use Cdata
      use Gravity

      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(mx,my,mz,mvar) :: df
      integer ider,j,k,i,ii
      real, dimension(nx) :: r,pdamp,aux0,velx0,vely0
      real :: tau
      
      tau = 2*pi/(r_int)**(-1.5)
      r = sqrt(x(l1:l2)**2 + y(m)**2)
     
      !for 0.4 : 1 ; 0.5 : 0
      pdamp = -11.111111*r**2 + 2.77777778  !parabolic function R
      !for 0.4 : 0 ; 0.5 : 1
      !pdamp = 11.1111111*r**2 - 1.777777778  
      
      where (r .le. r_int) 
         pdamp = 1.
      endwhere
      where (r .ge. 0.5)
         pdamp = 0.
      endwhere
      
      aux0 = g0*r**(n_pot-2)*(r**n_pot+r0_pot**n_pot)**(-1./n_pot-1.) !omega**2
      
      velx0 =  -y(  m  ) * sqrt(aux0)   !initial conditions
      vely0 =  x(l1:l2) * sqrt(aux0)
      
      do i=l1,l2
         ii = i-l1+1
         if ((r(ii).le.0.5).and.(r(ii).ge.0.4)) then
            df(i,m,n,ilnrho) = df(i,m,n,ilnrho) - (f(i,m,n,ilnrho) - 1.       )/tau * pdamp(ii) 
            df(i,m,n,iux)    = df(i,m,n,iux)    - (f(i,m,n,iux)    - velx0(ii))/tau * pdamp(ii)
            df(i,m,n,iuy)    = df(i,m,n,iuy)    - (f(i,m,n,iuy)    - vely0(ii))/tau * pdamp(ii)
         endif
      enddo
      
     
     !outer boundary
     
     tau = 2*pi/(r_ext)**(-1.5)
     
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
  endmodule Planet
  
