! $Id: density.f90,v 1.38 2002-07-20 17:43:53 dobler Exp $

module Density

!  This module is used both for the initial condition and during run time.
!  It contains dlnrho_dt and init_lnrho, among other auxiliary routines.

  use Cparam
  use Cdata

  implicit none

  real :: cs0=1., rho0=1.
  real :: cs20, lnrho0
  real :: ampllnrho=0., gamma=5./3., widthlnrho=.1
  real :: rho_left=1., rho_right=1., cdiffrho=0.
  real :: cs2bot=1., cs2top=1., gamma1
  real :: radius_lnrho=.5,kx_lnrho=0.,ky_lnrho=0.,kz_lnrho=0.
  real :: mpoly=1.5
  real :: eps_planet=.5
  character (len=labellen) :: initlnrho='zero', initlnrho2='zero'

  namelist /density_init_pars/ &
       cs0,rho0,ampllnrho,gamma,initlnrho,widthlnrho, &
       rho_left,rho_right,cs2bot,cs2top, &
       initlnrho2,radius_lnrho,eps_planet, &
       mpoly, &
       kx_lnrho,ky_lnrho,kz_lnrho

  namelist /density_run_pars/ &
       cs0,rho0,gamma,cdiffrho,cs2bot,cs2top

  ! other variables (needs to be consistent with reset list below)
  integer :: i_eth=0,i_ekin=0,i_rhom=0

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
        print*, 'Register_dnesity:  nvar = ', nvar
        print*, 'ilnrho = ', ilnrho
      endif
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: density.f90,v 1.38 2002-07-20 17:43:53 dobler Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_hydro: nvar > mvar')
      endif
!
    endsubroutine register_density
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
      use Global
      use Gravity
      use Initcond
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz,tmp,pot,prof
      real :: lnrhoint,cs2int,pot0
!
!  Set default values for sound speed at top and bottom.
!  These may be updated in one of the following initialization routines.
!
      cs2top=cs20; cs2bot=cs20
!
!  different initializations of lnrho (called from start).
!  If initrho does't match, f=0 is assumed (default).
!
      lnrho0 = alog(rho0)
      select case(initlnrho)

      case('zero', '0'); f(:,:,:,ilnrho)=0.
      case('isothermal'); call isothermal(f)
      case('polytropic_simple'); call polytropic_simple(f)
      case('hydrostatic-z', '1'); print*,'use polytropic_simple instead!'
      case('gaussian-noise'); call gaunoise(ampllnrho,f,ilnrho,ilnrho)
      case('gaussian-noise-x')
        !
        !  noise, but just x-dependent
        !
        call gaunoise(ampllnrho,f,ilnrho,ilnrho)
        f(:,:,:,ilnrho)=spread(spread(f(:,4,4,ilnrho),2,my),3,mz) !(watch 1-D)

      case('rho-jump', '2')
        !
        !  density jump (for shocks?)
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
          call setup_grav()     ! get coefficients cpot(1:5)
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
        !  piecewise polytropic for solar convection stuff
        !
        if (lroot) print*,'piecewise polytropic vertical stratification (lnrho)'
        ! top region

        cs2int = cs0**2
        lnrhoint = lnrho0
        f(:,:,:,ilnrho) = lnrho0 ! just in case
        call polytropic_lnrho_z(f,mpoly2,zz,tmp,zref,z2,z0+2*Lz, &
                            isothtop,cs2int,lnrhoint)
          ! unstable layer
        call polytropic_lnrho_z(f,mpoly0,zz,tmp,z2,z1,z2,0,cs2int,lnrhoint)
          ! stable layer
        call polytropic_lnrho_z(f,mpoly1,zz,tmp,z1,z0,z1,0,cs2int,lnrhoint)
        !
        !  calculate cs2bot and cs2top for run.x (boundary conditions)
        !
        cs2bot = cs2int + gamma/gamma1*gravz/(mpoly2+1)*(z(n1)-z0  )
        if (isothtop /= 0) then
          cs2top = cs20
        else
          cs2top = cs20 + gamma/gamma1*gravz/(mpoly0+1)*(z(n2)-zref)
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
        cs2bot = cs20 + gamma*gravz/(mpoly0+1)*(z(n1)-zref)
        cs2top = cs20 + gamma*gravz/(mpoly0+1)*(z(n2)-zref)

      case('sound-wave', '11')
        !
        !  sound wave (should be consistent with hydro module)
        !
        if (lroot) print*,'x-wave in lnrho; ampllnrho=',ampllnrho
        f(:,:,:,ilnrho)=ampllnrho*sin(xx)

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
        call planet(ampllnrho,f,xx,yy,zz,eps_planet,radius_lnrho,gamma)

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
      real, dimension (mx,my,mz,mvar) :: f
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
      stp = step(z,zblend,whcond)
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
      real, dimension (mx,my,mz,mvar) :: f,df
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
        i_eth=0;i_ekin=0;i_rhom=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'eth',i_eth)
        call parse_name(iname,cname(iname),cform(iname),'ekin',i_ekin)
        call parse_name(iname,cname(iname),cform(iname),'rhom',i_rhom)
      enddo
!
!  write column where which magnetic variable is stored
!
      write(3,*) 'i_eth=',i_eth
      write(3,*) 'i_ekin=',i_ekin
      write(3,*) 'i_rhom=',i_rhom
      write(3,*) 'nname=',nname
      write(3,*) 'ilnrho=',ilnrho
!
    endsubroutine rprint_density
!***********************************************************************
    subroutine bc_lnrho_db(f,topbot)
!
!  boundary condition for density
!
!  may-2002/nils: coded
!  11-jul-2002/nils: moved into the density module
!
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my) :: fder
      integer :: i
      !
      !  Set ghost zone to reproduce one-sided boundary condition 
      !  (2nd order):
      !  Finding the derivatives on the boundary using a one 
      !  sided final difference method. This derivative is being 
      !  used to calculate the boundary points. This will probably
      !  only be used for ln(rho)
      !
    select case(topbot)
!
! Bottom boundary
!
    case('bot')
       do i=1,nghost
          fder=(-3*f(:,:,n1-i+1,ilnrho)+4*f(:,:,n1-i+2,ilnrho)&
               -f(:,:,n1-i+3,ilnrho))/(2*dz)
          f(:,:,n1-i,ilnrho)=f(:,:,n1-i+2,ilnrho)-2*dz*fder
       end do
    case('top')
       do i=1,nghost
          fder=(3*f(:,:,n2+i-1,ilnrho)-4*f(:,:,n2+i-2,ilnrho)&
               +f(:,:,n2+i-3,ilnrho))/(2*dz)
          f(:,:,n2+i,ilnrho)=f(:,:,n2+i-2,ilnrho)+2*dz*fder
       end do
    case default
       if(lroot) print*,"invalid argument for 'bc_ss_flux'"
    endselect
!
  end subroutine bc_lnrho_db
!***********************************************************************
!  Here comes a collection of different density stratification routines
!***********************************************************************
    subroutine isothermal(f)
!
!  Isothermal stratification (for lnrho and ss)
!  This routine should be independent of the gravity module used.
!  When entropy is present, this module also initializes entropy.
!
!   8-jul-02/axel: incorporated/adapted from init_lnrho
!  11-jul-02/axel: fixed sign; should be tmp=gamma*pot/cs20
!
      use Gravity
!
      real, dimension (mx,my,mz,mvar) :: f
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
        if(lentropy) f(l1:l2,m,n,ient)=-gamma1/gamma*tmp
      enddo
      enddo

!
!  cs2 values at top and bottom may be needed to boundary conditions.
!  The values calculated here may be revised in the entropy module.
!
      cs2bot=cs20
      cs2top=cs20
!
    endsubroutine isothermal
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
      real, dimension (mx,my,mz,mvar) :: f
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
        zinfty=zref+(mpoly+1.)*cs20/(-gamma*gravz)
      elseif (grav_profile=='linear') then
        zinfty2=zref**2+(mpoly+1.)*cs20/(-.5*gamma*gravz)
        if(zinfty2<0) then
          if(lroot) print*,'polytropic_simple: zinfty**2<0 is not ok'
          zinfty2=0. !(see what happens)
        endif
        zinfty=sqrt(zinfty2)
      else
        if(lroot) print*,'polytropic_simple: zinfty not prepared!'
      endif
!
!  check whether zinfty lies outside the domain (otherwise density
!  would vanish within the domain)
!
      ztop=xyz0(3)+Lxyz(3)
      zbot=xyz0(3)
      if(zinfty<ztop .or. (-zinfty)>zbot) then
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

endmodule Density
