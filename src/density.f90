! $Id: density.f90,v 1.91 2003-06-16 04:41:10 brandenb Exp $

!  This module is used both for the initial condition and during run time.
!  It contains dlnrho_dt and init_lnrho, among other auxiliary routines.

module Density

  use Cparam
  use Cdata

  implicit none

  real :: cs0=1., rho0=1.
  real :: cs20, lnrho0
  real :: ampllnrho=0., gamma=5./3., widthlnrho=.1
  real :: rho_left=1., rho_right=1., cdiffrho=0., lnrho_const=0.
  real :: cs2bot=1., cs2top=1., gamma1,amplrho=0
  real :: radius_lnrho=.5,kx_lnrho=1.,ky_lnrho=1.,kz_lnrho=1.
  real :: eps_planet=.5
  real :: b_ell=1., q_ell=5., hh0=0., rbound=1.
  real :: mpoly=1.5
  real :: mpoly0=1.5,mpoly1=1.5,mpoly2=1.5
  real :: frec_lnrho=1,ampl_osc_lnrho=1e-3
  integer:: isothtop=0
  character (len=labellen) :: initlnrho='zero', initlnrho2='zero'

  namelist /density_init_pars/ &
       cs0,rho0,ampllnrho,gamma,initlnrho,widthlnrho, &
       rho_left,rho_right,lnrho_const,cs2bot,cs2top, &
       initlnrho2,radius_lnrho,eps_planet, &
       b_ell,q_ell,hh0,rbound, &
       mpoly, &
       kx_lnrho,ky_lnrho,kz_lnrho,amplrho

  namelist /density_run_pars/ &
       cs0,rho0,gamma,cdiffrho,cs2bot,cs2top,frec_lnrho,ampl_osc_lnrho

  ! diagnostic variables (needs to be consistent with reset list below)
  integer :: i_ekin=0,i_rhom=0,i_ekintot=0

  contains

!***********************************************************************
    subroutine register_density()
!
!  Initialise variables which should know that we solve the
!  compressible hydro equations: ilnrho; increase nvar accordingly.
!
!   4-jun-02/axel: adapted from hydro
!
      use Mpicomm, only: stop_it
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_density called twice')
      first = .false.
!
      ldensity = .true.
!
      ilnrho = nvar+1           ! indix to access lam
      nvar = nvar+1             ! added 1 variable
!
      if ((ip<=8) .and. lroot) then
        print*, 'Register_density:  nvar = ', nvar
        print*, 'ilnrho = ', ilnrho
      endif
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: density.f90,v 1.91 2003-06-16 04:41:10 brandenb Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_hydro: nvar > mvar')
      endif
!
    endsubroutine register_density
!***********************************************************************
    subroutine initialize_density()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded 
!
!  do nothing
!
    endsubroutine initialize_density
!***********************************************************************
    subroutine init_lnrho(f,xx,yy,zz)
!
!  initialise lnrho; called from start.f90
!
!  7-nov-01/wolf: coded
! 28-jun-02/axel: added isothermal
!
      use Mpicomm
      use Sub
      use IO
      use Global
      use Gravity
      use Initcond
      use Initcond_spec
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz,tmp,pot,prof
      real :: lnrhoint,cs2int,pot0,lnrho_left,lnrho_right
      real :: zbot,ztop
!
!  define bottom and top height
!
      zbot=xyz0(3)
      ztop=xyz0(3)+Lxyz(3)
!
!  Set default values for sound speed at top and bottom.
!  These may be updated in one of the following initialization routines.
!
      cs2top=cs20; cs2bot=cs20
!
!  different initializations of lnrho (called from start).
!  If initrho does't match, f=0 is assumed (default).
!
      lnrho0      = alog(rho0)
      lnrho_left  = alog(rho_left)
      lnrho_right = alog(rho_right)
      select case(initlnrho)

      case('zero', '0'); if(lroot) print*,'zero lnrho'
      case('const_lnrho'); f(:,:,:,ilnrho)=lnrho_const
      case('constant'); f(:,:,:,ilnrho)=alog(rho_left)
      case('isothermal'); call isothermal_density(f)
      case('stratification'); call stratification(amplrho,f,xx,yy,zz)
      case('polytropic_simple'); call polytropic_simple(f)
      case('hydrostatic-z', '1'); print*,'use polytropic_simple instead!'
      case('xjump'); call jump(f,ilnrho,lnrho_left,lnrho_right,widthlnrho,'x')
      case('yjump'); call jump(f,ilnrho,lnrho_left,lnrho_right,widthlnrho,'y')
      case('zjump'); call jump(f,ilnrho,lnrho_left,lnrho_right,widthlnrho,'z')
      case('gaussian3d'); call gaussian3d(ampllnrho,f,ilnrho,xx,yy,zz,radius_lnrho)
      case('gaussian-noise'); call gaunoise(ampllnrho,f,ilnrho,ilnrho)
      case('gaussian-noise-x')
        !
        !  noise, but just x-dependent
        !
        call gaunoise(ampllnrho,f,ilnrho,ilnrho)
        f(:,:,:,ilnrho)=spread(spread(f(:,4,4,ilnrho),2,my),3,mz) !(watch 1-D)

      case('rho-jump-z', '2')
        !
        !  density jump (for shocks)
        !
        if (lroot) print*,'density jump; rho_left,right=',rho_left,rho_right
        if (lroot) print*,'density jump; widthlnrho=',widthlnrho
        prof=.5*(1.+tanh(zz/widthlnrho))
        f(:,:,:,ilnrho)=alog(rho_left)+alog(rho_left/rho_right)*prof

      case ('hydrostatic-z-2', '3')
        !
        !  hydrostatic density stratification for isentropic atmosphere
        !
        if (lgravz) then
          if (lroot) print*,'vertical density stratification'
          !        f(:,:,:,ilnrho)=-zz
          ! isentropic case:
          !        zmax = -cs20/gamma1/gravz
          !        print*, 'zmax = ', zmax
          !        f(:,:,:,ilnrho) = 1./gamma1 * alog(abs(1-zz/zmax))
          ! linear entropy gradient;
          f(:,:,:,ilnrho) = -grads0*zz &
                            + 1./gamma1*alog( 1 + gamma1*gravz/grads0/cs20 &
                                                  *(1-exp(-grads0*zz)) )
          if (notanumber(f(:,:,:,ilnrho))) then
            STOP "INIT_LNRHO: Imaginary density values"
          endif
        endif

      case ('hydrostatic-r')
        !
        !  hydrostatic radial density stratification for isentropic (or
        !  isothermal) sphere
        !
        if (lgravr) then
          if (lroot) print*, &
               'radial density stratification (assumes s=const)'
!ajwm - here's the init call that needs sorting!
          call initialize_gravity()     ! get coefficients cpot(1:5)

          call potential(xx,yy,zz,pot,POT0=pot0) ! gravity potential
          call output(trim(directory)//'/pot.dat',pot,1)
          !
          ! rho0, cs0 are the values in the centre, pot0
          !
          if (gamma /= 1) then  ! isentropic
            f(:,:,:,ilnrho) = lnrho0 &
                              + alog(1 - gamma1*(pot-pot0)/cs20) / gamma1
          else                  ! isothermal
            f(:,:,:,ilnrho) = lnrho0 - (pot-pot0)/cs20
          endif
        endif

      case ('piecew-poly', '4')
        !
        !  piecewise polytropic for stellar convection models
        !
        if (lroot) print*,'piecewise polytropic vertical stratification (lnrho)'
        ! top region

        cs2int = cs0**2
        lnrhoint = lnrho0
        f(:,:,:,ilnrho) = lnrho0 ! just in case
        call polytropic_lnrho_z(f,mpoly2,zz,tmp,zref,z2,ztop+Lz, &
                            isothtop,cs2int,lnrhoint)
          ! unstable layer
        call polytropic_lnrho_z(f,mpoly0,zz,tmp,z2,z1,z2,0,cs2int,lnrhoint)
          ! stable layer
        call polytropic_lnrho_z(f,mpoly1,zz,tmp,z1,z0,z1,0,cs2int,lnrhoint)
        !
        !  calculate cs2bot and cs2top for run.x (boundary conditions)
        !
        cs2bot = cs2int + gamma/gamma1*gravz/(mpoly2+1)*(zbot-z0  )
        if (isothtop /= 0) then
          cs2top = cs20
        else
          cs2top = cs20 + gamma/gamma1*gravz/(mpoly0+1)*(ztop-zref)
        endif

      case ('polytropic', '5')
        !
        !  polytropic stratification
        !  cs0, rho0 and ss0=0 refer to height z=zref
        !
        if (lroot) print*,'polytropic vertical stratification (lnrho)'
        !
        cs20 = cs0**2
        cs2int = cs20
        lnrhoint = lnrho0
        f(:,:,:,ilnrho) = lnrho0 ! just in case
        ! only one layer
        call polytropic_lnrho_z(f,mpoly0,zz,tmp,zref,z0,z0+2*Lz, &
             0,cs2int,lnrhoint)
        !
        !  calculate cs2bot and cs2top for run.x (boundary conditions)
        !
        cs2bot = cs20 + gamma*gravz/(mpoly0+1)*(zbot-zref)
        cs2top = cs20 + gamma*gravz/(mpoly0+1)*(ztop-zref)

      case('sound-wave', '11')
        !
        !  sound wave (should be consistent with hydro module)
        !
        if (lroot) print*,'x-wave in lnrho; ampllnrho=',ampllnrho
        f(:,:,:,ilnrho)=lnrho_const+ampllnrho*sin(kx_lnrho*xx)

      case('sound-wave2')
        !
        !  sound wave (should be consistent with hydro module)
        !
        if (lroot) print*,'x-wave in lnrho; ampllnrho=',ampllnrho
        f(:,:,:,ilnrho)=lnrho_const+ampllnrho*cos(kx_lnrho*xx)

      case('shock-tube', '13')
        !
        !  shock tube test (should be consistent with hydro module)
        !  
        if (lroot) print*,'polytopic standing shock'
        prof=.5*(1.+tanh(xx/widthlnrho))
        f(:,:,:,ilnrho)=alog(rho_left)+(alog(rho_right)-alog(rho_left))*prof

      case('sin-xy')
        !
        !  sin profile in x and y
        !  
        if (lroot) print*,'lnrho: sin(x)*sin(y)'
        f(:,:,:,ilnrho) = &
             alog(rho0) + ampllnrho*sin(kx_lnrho*xx)*sin(ky_lnrho*yy)

      case('sin-xy-rho')
        !
        !  sin profile in x and y, but in rho, not ln(rho)
        !  
        if (lroot) print*,'rho: sin(x)*sin(y)'
        f(:,:,:,ilnrho) = &
             alog(rho0*(1+ampllnrho*sin(kx_lnrho*xx)*sin(ky_lnrho*yy)))

      case('linear')
        !
        !  linear profile in kk.xxx
        !  
        if (lroot) print*,'lnrho: linear profile'
        f(:,:,:,ilnrho) = alog(rho0) &
             + ampllnrho*(kx_lnrho*xx+ky_lnrho*yy+kz_lnrho*zz) &
               / sqrt(kx_lnrho**2+ky_lnrho**2+kz_lnrho**2)

      case('planet')
        !
        !  planet solution of Goodman, Narayan & Goldreich (1987)
        !
        call planet(rbound,f,xx,yy,zz,eps_planet,radius_lnrho,gamma,cs20,widthlnrho,hh0)

      case('kepvor')
        !
        !  planet solution of Chavanis (2000)
        !
        call kepvor(f,xx,yy,zz,b_ell,q_ell,gamma,cs20,hh0)


      case('enthblob')
        !
        !  attempt to produce vortex from enthalpy blob
        !
        call enthblob(f,xx,yy,zz,b_ell,q_ell,gamma,cs20,hh0)


      case('Ferriere'); call ferriere(f)

      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*,'No such value for initlnrho: ', trim(initlnrho)
        call stop_it("")

      endselect
!
!  sanity check
!
      if (notanumber(f(:,:,:,ilnrho))) then
        STOP "INIT_LNRHO: Imaginary density values"
      endif
!
!  check that cs2bot,cs2top are ok
!
      if(lroot) print*,'init_lnrho: cs2bot,cs2top=',cs2bot,cs2top
      if(lroot) print*,'e.g. for ionization runs: cs2bot,cs2top not yet set' 
!
!  different initializations of lnrho (called from start).
!  If initrho does't match, f=0 is assumed (default).
!
      select case(initlnrho2)

      case('addblob')
        !
        if (lroot) print*,'init_lnrho: add blob'
        f(:,:,:,ilnrho)=f(:,:,:,ilnrho) &
          +ampllnrho*exp(-(xx**2+yy**2+zz**2)/radius_lnrho**2)
      endselect
!
      if(ip==0) print*, yy(1,1,1) ! keep compiler quiet
!
    endsubroutine init_lnrho
!***********************************************************************
    subroutine polytropic_lnrho_z( &
         f,mpoly,zz,tmp,zint,zbot,zblend,isoth,cs2int,lnrhoint)
!
!  Implement a polytropic profile in ss above zbot. If this routine is
!  called several times (for a piecewise polytropic atmosphere), on needs
!  to call it from top to bottom.
!
!  zint    -- height of (previous) interface, where cs2int and lnrhoint
!             are set
!  zbot    -- z at bottom of layer
!  zblend  -- smoothly blend (with width whcond) previous ss (for z>zblend)
!             with new profile (for z<zblend)
!  isoth   -- flag for isothermal stratification;
!  lnrhoint -- value of lnrho at the interface, i.e. at the zint on entry,
!             at the zbot on exit
!  cs2int  -- same for cs2
!
      use Sub, only: step
      use Gravity, only: gravz
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: tmp,p,zz
      real, dimension (mz) :: stp
      real :: mpoly,zint,zbot,zblend,beta1,cs2int,lnrhoint
      integer :: isoth
!
      ! NB: beta1 is not dT/dz, but dcs2/dz = (gamma-1)c_p dT/dz
      if (isoth /= 0) then ! isothermal layer
        beta1 = 0.
        tmp = gamma*gravz/cs2int*(zz-zint)
        lnrhoint = lnrhoint + gamma*gravz/cs2int*(zbot-zint)
      else
        beta1 = gamma*gravz/(mpoly+1)
        tmp = 1 + beta1*(zz-zint)/cs2int
        tmp = max(tmp,epsi)  ! ensure arg to log is positive
        tmp = lnrhoint + mpoly*alog(tmp)
        lnrhoint = lnrhoint + mpoly*alog(1 + beta1*(zbot-zint)/cs2int)
      endif
      cs2int = cs2int + beta1*(zbot-zint) ! cs2 at layer interface (bottom)
      !
      ! smoothly blend the old value (above zblend) and the new one (below
      ! zblend) for the two regions:
      !
      stp = step(z,zblend,widthlnrho)
      p = spread(spread(stp,1,mx),2,my)
      f(:,:,:,ilnrho) = p*f(:,:,:,ilnrho) + (1-p)*tmp
!
    endsubroutine polytropic_lnrho_z
!***********************************************************************
    subroutine dlnrho_dt(f,df,uu,glnrho,divu,lnrho)
!
!  continuity equation
!  calculate dlnrho/dt = - u.gradlnrho - divu
!
!   7-jun-02/axel: incoporated from subroutine pde
!
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f,df
      real, dimension (nx,3) :: uu,glnrho
      real, dimension (nx) :: lnrho,divu,uglnrho,glnrho2
      real, dimension (nx) :: del2lnrho
      real :: diffrho
!
      intent(in)  :: f,uu,divu
      intent(out) :: df,glnrho,lnrho
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'SOLVE dlnrho_dt'
      if (headtt) call identify_bcs('lnrho',ilnrho)
!
!  define lnrho; calculate density gradient and avection term
!
      lnrho=f(l1:l2,m,n,ilnrho)
      call grad(f,ilnrho,glnrho)
      call dot_mn(uu,glnrho,uglnrho)
!
!  continuity equation
!
      df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)-uglnrho-divu
!
!  mass diffusion, in units of dxmin*cs0
!
      if (cdiffrho /= 0.) then
        diffrho=cdiffrho*dxmin*cs0
        call del2(f,ilnrho,del2lnrho)
        call dot2_mn(glnrho,glnrho2)
        df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)+diffrho*(del2lnrho+glnrho2)
        maxdiffus=amax1(maxdiffus,diffrho)
      endif
!
    endsubroutine dlnrho_dt
!***********************************************************************
    subroutine rprint_density(lreset)
!
!  reads and registers print parameters relevant for compressible part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Sub
!
      integer :: iname
      logical :: lreset
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_ekin=0; i_rhom=0; i_ekintot=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ekintot',i_ekintot)
        call parse_name(iname,cname(iname),cform(iname),'ekin',i_ekin)
        call parse_name(iname,cname(iname),cform(iname),'rhom',i_rhom)
      enddo
!
!  write column where which magnetic variable is stored
!
      write(3,*) 'i_ekintot=',i_ekintot
      write(3,*) 'i_ekin=',i_ekin
      write(3,*) 'i_rhom=',i_rhom
      write(3,*) 'nname=',nname
      write(3,*) 'ilnrho=',ilnrho
!
    endsubroutine rprint_density
!***********************************************************************
!  Here comes a collection of different density stratification routines
!***********************************************************************
    subroutine isothermal_density(f)
!
!  Isothermal stratification (for lnrho and ss)
!  This routine should be independent of the gravity module used.
!  When entropy is present, this module also initializes entropy.
!
!  Sound speed (and hence Temperature), and density (at infinity) are 
!  initialised to their respective reference values:
!           sound speed: cs^2_0            from start.in  
!           density: rho0 = exp(lnrho0)
!
!   8-jul-02/axel: incorporated/adapted from init_lnrho
!  11-jul-02/axel: fixed sign; should be tmp=gamma*pot/cs20
!  02-apr-03/tony: made entropy explicit rather than using tmp/-gamma  
!  11-jun-03/tony: moved entropy initialisation to separate routine 
!                  to allow isothermal condition for arbitrary density
!
      use Gravity
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) :: pot,tmp
!
!  Stratification depends on the gravity potential
!
      if (lroot) print*,'isothermal stratification'
      do n=n1,n2
      do m=m1,m2
        call potential(x(l1:l2),y(m),z(n),pot)
        tmp=-gamma*pot/cs20
        f(l1:l2,m,n,ilnrho)=lnrho0+tmp
!        if(lentropy) f(l1:l2,m,n,ient)= -gamma1/gamma*tmp
!                                      = gamma1*pot/cs20
!   MOVED to isothermal_entropy routine

        if(lentropy) f(l1:l2,m,n,ient)= &
             -gamma1*(f(l1:l2,m,n,ilnrho)-lnrho0)/gamma
      enddo
      enddo

!
!  cs2 values at top and bottom may be needed to boundary conditions.
!  The values calculated here may be revised in the entropy module.
!
      cs2bot=cs20
      cs2top=cs20
!
    endsubroutine isothermal_density
!***********************************************************************
    subroutine polytropic_simple(f)
!
!  Polytropic stratification (for lnrho and ss)
!  This routine should be independent of the gravity module used.
!
!  To maintain continuity with respect to the isothermal case,
!  one may want to specify cs20 (=1), and so zinfty is calculated from that.
!  On the other hand, for polytropic atmospheres it may be more
!  physical to specify zinfty (=1), ie the top of the polytropic atmosphere.
!  This is done if zinfty is different from 0.
!
!   8-jul-02/axel: incorporated/adapted from init_lnrho
!
      use Gravity
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) :: pot,dlncs2,ptop,pbot,zero=0.
      real :: ggamma,ztop,zbot,zinfty2
!
!  identifier
!
      if (lroot) print*,'polytropic_simple: mpoly=',mpoly
!
!  zinfty is calculated such that rho=rho0 and cs2=cs20 at z=zref.
!  Note: gravz is normally negative!
!
      if (grav_profile=='const') then
        if(lroot.and.gravz==0.) print*,'polytropic_simple: divide by gravz=0'
        zinfty=zref+(mpoly+1.)*cs20/(-gamma*gravz)
      elseif (grav_profile=='const_zero') then
        if(lroot.and.gravz==0.) print*,'polytropic_simple: divide by gravz=0'
        zinfty=zref+(mpoly+1.)*cs20/(-gamma*gravz)
      elseif (grav_profile=='linear') then
        if(lroot.and.gravz==0.) print*,'polytropic_simple: divide by gravz=0'
        zinfty2=zref**2+(mpoly+1.)*cs20/(-.5*gamma*gravz)
        if(zinfty2<0) then
          if(lroot) print*,'polytropic_simple: zinfty**2<0 is not ok'
          zinfty2=0. !(and see what happens)
        endif
        zinfty=sqrt(zinfty2)
      else
        if(lroot) print*,'polytropic_simple: zinfty not prepared!'
      endif
!
!  check whether zinfty lies outside the domain (otherwise density
!  would vanish within the domain). At the moment we are not properly
!  testing the lower boundary on the case of a disk (commented out).
!
      ztop=xyz0(3)+Lxyz(3)
      zbot=xyz0(3)
      !-- if(zinfty<amin1(ztop,zgrav) .or. (-zinfty)>amin1(zbot,zgrav)) then
      if(zinfty<amin1(ztop,zgrav)) then
        if(lroot) print*,'polytropic_simple: domain too big; zinfty=',zinfty
        call stop_it('rho and cs2 will vanish within domain')
      endif
!
      ggamma=1.+1./mpoly
!
      do n=n1,n2
      do m=m1,m2
        call potential(x(l1:l2),y(m),z(n),pot)
        dlncs2=alog(-gamma*pot/((mpoly+1.)*cs20))
        f(l1:l2,m,n,ilnrho)=lnrho0+mpoly*dlncs2
        if(lentropy) f(l1:l2,m,n,ient)=mpoly*(ggamma/gamma-1.)*dlncs2
      enddo
      enddo
!
!  cs2 values at top and bottom may be needed to boundary conditions.
!  In spherical geometry, ztop is z at the outer edge of the box,
!  so this calculation still makes sense.
!
      call potential(zero,0.,ztop,ptop)
      cs2top=-gamma/(mpoly+1.)*ptop(1)
!
!  In spherical geometry ztop should never be used.
!  Even in slab geometry ztop is not normally used.
!
      call potential(zero,0.,zbot,pbot)
      cs2bot=-gamma/(mpoly+1.)*pbot(1)
!
    endsubroutine polytropic_simple
!***********************************************************************
    subroutine ferriere(f)
!
!  density profile from K. Ferriere, ApJ 497, 759, 1998, 
!   eqns (6), (7), (9), (13), (14) [with molecular H, n_m, neglected]
!   at solar radius.  (for interstellar runs)
!  entropy is set via pressure, assuming a constant T for each gas component
!   (cf. eqn 15) and crudely compensating for non-thermal sources.
!  [an alternative treatment of entropy, based on hydrostatic equilibrium,
!   might be preferable. This was implemented in serial (in e.g. r1.59)
!   but abandoned as overcomplicated to adapt for nprocz /= 0.]
!
      use Mpicomm
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: absz,n_c,n_w,n_i,n_h
!  T in K, k_B s.t. pp is in code units ( = 9.59e-15 erg/cm/s^2)
!  (i.e. k_B = 1.381e-16 (erg/K) / 9.59e-15 (erg/cm/s^2) )
      real :: T_c=500.0,T_w=8.0e3,T_i=8.0e3,T_h=1.0e6,k_B=0.0144
      real :: rho,lnrho,pp,pp0 
      real, dimension(2) :: fmpi2
!      real, dimension(1) :: fmpi1
!      integer :: iproctop
!
      if (lroot) print*,'Ferriere density profile'
!
!  first define reference values of pp, cs2, at midplane.  
!  pressure is set to 6 times thermal pressure, this factor roughly
!  allowing for other sources, as modelled by Ferriere.
!
      pp0=6.0*k_B*(rho0/1.38) *                                               &
       (1.09*0.340*T_c + 1.09*0.226*T_w + 2.09*0.025*T_i + 2.27*0.00048*T_h)
      cs20=gamma*pp0/rho0
      cs0=sqrt(cs20)
!      ss0=alog(gamma*pp0/cs20/rho0)/gamma   !ss0=zero  (not needed)
!
      do n=n1,n2            ! nb: don't need to set ghost-zones here
      absz=abs(z(n))
      do m=m1,m2
!  cold gas profile n_c (eq 6)
        n_c=0.340*(0.859*exp(-min((z(n)/0.127)**2,70.)) +         &
                   0.047*exp(-min((z(n)/0.318)**2,70.)) +         &
                   0.094*exp(-min(absz/0.403,70.)))     
!  warm gas profile n_w (eq 7)
        n_w=0.226*(0.456*exp(-min((z(n)/0.127)**2,70.)) +         &
                   0.403*exp(-min((z(n)/0.318)**2,70.)) +         &
                   0.141*exp(-min(absz/0.403,70.)))
!  ionized gas profile n_i (eq 9)
        n_i=0.0237*exp(-absz) + 0.0013* exp(-min(absz/0.150,70.))
!  hot gas profile n_h (eq 13)
        n_h=0.00048*exp(-absz/1.5)         
!  normalised s.t. rho0 gives mid-plane density directly (in 10^-24 g/cm^3)
        rho=rho0/(0.340+0.226+0.025+0.00048)*(n_c+n_w+n_i+n_h)
        lnrho=alog(rho)
        f(l1:l2,m,n,ilnrho)=lnrho
!  define entropy via pressure, assuming fixed T for each component
        if(lentropy) then
!  thermal pressure (eq 13)
          pp=6.0*k_B*(rho0/1.38) *                                        &
           (1.09*n_c*T_c + 1.09*n_w*T_w + 2.09*n_i*T_i + 2.27*n_h*T_h)
          f(l1:l2,m,n,ient)=alog(gamma*pp/cs20)/gamma +                   &
                                     gamma1/gamma*lnrho0 - lnrho
!  calculate cs2bot,top: needed for a2/c2 b.c.s (fixed T)
          if (n == n1 .and. m == m1) cs2bot=gamma*pp/rho
          if (n == n2 .and. m == m1) cs2top=gamma*pp/rho
        endif
       enddo
      enddo
!
!  broadcast cs2bot, top
!
      if (lentropy) then
!  just use symmetry to set cs2top=cs2bot, and broadcast from root
        cs2top=cs2bot
        fmpi2=(/ cs2bot, cs2top /)
        call mpibcast_real(fmpi2,2)
        cs2bot=fmpi2(1); cs2top=fmpi2(2)
!!  or do directly from the right processor
!        fmpi1=(/ cs2bot /)
!        call mpibcast_real(fmpi1,1)    ! this case can use mpibcast_real
!        cs2bot=fmpi1(1)
!        iproctop=(nprocz-1)*nprocy     ! ipz=nprocz-1,ipy=0
!        fmpi1=(/ cs2top /)
!        call mpibcast_real_nonroot(fmpi1,1,iproctop)
!        cs2top=fmpi1(1)
      endif
      if (lroot) print*, 'ferriere: cs2bot=',cs2bot, ' cs2top=',cs2top
!
    endsubroutine ferriere
!***********************************************************************
    subroutine bc_lnrho_temp_z(f,topbot)
!
!  boundary condition for lnrho *and* ss: constant temperature
!
!  27-sep-2002/axel: coded
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: tmp
      integer :: i
!
      if(ldebug) print*,'ENTER: bc_lnrho_temp_z, cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (lroot) print*,'bc_lnrho_temp_z: bot not yet implemented'
        call stop_it("")
!
!  top boundary
!
      case('top')
        if (ldebug) print*,'set z top temperature: cs2top=',cs2top
        if (cs2top<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
!
!  set first boundary value for entropy, and then ghost points antisymmetrically
!
        f(:,:,n2,ient) = 0.5*tmp - gamma1/gamma*(f(:,:,n2,ilnrho)-lnrho0)
        do i=1,nghost; f(:,:,n2+i,ient) = 2*f(:,:,n2,ient)-f(:,:,n2-i,ient); enddo
!
!  set density in the ghost zones so that dlnrho/dz + ds/dz = gz/cs2top
!  for the time being, we don't worry about lnrho0 (assuming that it is 0)
!
        tmp=gravz/cs2top
        do i=1,nghost
          f(:,:,n2+i,ilnrho) = f(:,:,n2-i,ilnrho) +f(:,:,n2-i,ient) &
                                                  -f(:,:,n2+i,ient) +2*i*dz*tmp
        enddo

      case default
        if(lroot) print*,"invalid argument for 'bc_ss_temp_z'"
        call stop_it("")
      endselect
!
    endsubroutine bc_lnrho_temp_z
!***********************************************************************
    subroutine bc_lnrho_pressure_z(f,topbot)
!
!  boundary condition for lnrho: constant pressure
!
!   4-apr-2003/axel: coded
!   1-may-2003/axel: added the same for top boundary
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: i,l
!
      if(ldebug) print*,'ENTER: bc_lnrho_temp_z, cs20,cs0=',cs20,cs0
!
!  Constant pressure, i.e. antisymmetric
!  This assumes that the entropy is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('top')
        if (ldebug) print*,'bc_lnrho_pressure_z: lnrho_top,ss_top=',lnrho_top,ss_top
!
!  fix entropy if inflow (uz>0); otherwise leave s unchanged
!  afterwards set s antisymmetrically about boundary value
!
        if(lentropy) then
!         do m=m1,m2
!         do l=l1,l2
!           if (f(l,m,n1,iuz)>=0) then
!             f(l,m,n1,ient)=ss_bot
!           else
!             f(l,m,n1,ient)=f(l,m,n1+1,ient)
!           endif
!         enddo
!         enddo
          f(:,:,n2,ient)=ss_top
          do i=1,nghost; f(:,:,n2+i,ient)=2*f(:,:,n2,ient)-f(:,:,n2-i,ient); enddo
!
!  set density value such that pressure is constant at the bottom
!
          f(:,:,n2,ilnrho)=lnrho_top+ss_top-f(:,:,n2,ient)
        else
          f(:,:,n2,ilnrho)=lnrho_top
        endif
!
!  make density antisymmetric about boundary
!  another possibility might be to enforce hydrostatics
!  ie to set dlnrho/dz=-g/cs^2, assuming zero entropy gradient
!
        do i=1,nghost
          f(:,:,n2+i,ilnrho)=2*f(:,:,n2,ilnrho)-f(:,:,n2-i,ilnrho)
        enddo
!
!  top boundary
!
      case('bot')
        if (ldebug) print*,'bc_lnrho_pressure_z: lnrho_bot,ss_bot=',lnrho_bot,ss_bot
!
!  fix entropy if inflow (uz>0); otherwise leave s unchanged
!  afterwards set s antisymmetrically about boundary value
!
        if(lentropy) then
!         do m=m1,m2
!         do l=l1,l2
!           if (f(l,m,n1,iuz)>=0) then
!             f(l,m,n1,ient)=ss_bot
!           else
!             f(l,m,n1,ient)=f(l,m,n1+1,ient)
!           endif
!         enddo
!         enddo
          f(:,:,n1,ient)=ss_bot
          do i=1,nghost; f(:,:,n1-i,ient)=2*f(:,:,n1,ient)-f(:,:,n1+i,ient); enddo
!
!  set density value such that pressure is constant at the bottom
!
          f(:,:,n1,ilnrho)=lnrho_bot+ss_bot-f(:,:,n1,ient)
        else
          f(:,:,n1,ilnrho)=lnrho_bot
        endif
!
!  make density antisymmetric about boundary
!  another possibility might be to enforce hydrostatics
!  ie to set dlnrho/dz=-g/cs^2, assuming zero entropy gradient
!
        do i=1,nghost
          f(:,:,n1-i,ilnrho)=2*f(:,:,n1,ilnrho)-f(:,:,n1+i,ilnrho)
        enddo
!
      case default
        if(lroot) print*,"invalid argument for 'bc_lnrho_pressure_z'"
        call stop_it("")
      endselect
!
    endsubroutine bc_lnrho_pressure_z
!***********************************************************************

endmodule Density
