! $Id: density.f90,v 1.161 2004-05-18 15:38:24 dobler Exp $

!  This module is used both for the initial condition and during run time.
!  It contains dlnrho_dt and init_lnrho, among other auxiliary routines.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Density

  use Cparam
  use Cdata
  use Ionization, only: cs0,cs20,lnrho0,rho0,lcalc_cp,gamma,gamma1,cs2top,cs2bot,cp

  implicit none


  real :: ampllnrho=0., widthlnrho=.1
  real :: rho_left=1., rho_right=1., cdiffrho=0., diffrho=0., diffrho_shock=0.
  real :: lnrho_const=0., rho_const=1.
  real :: amplrho=0,cs2cool=0.
  real :: radius_lnrho=.5,kx_lnrho=1.,ky_lnrho=1.,kz_lnrho=1.
  real :: eps_planet=.5
  real :: b_ell=1., q_ell=5., hh0=0., rbound=1.
  real :: dlnrhobdx=0.,co1_ss=0.,co2_ss=0.,Sigma1=150.
  real :: mpoly=1.5
  real :: mpoly0=1.5,mpoly1=1.5,mpoly2=1.5
  real, dimension(3) :: gradlnrho0=(/0.,0.,0./)
  integer :: isothtop=0
  logical :: lupw_lnrho=.false.
  character (len=labellen), dimension(ninit) :: initlnrho='nothing'
  character (len=labellen) :: strati_type='lnrho_ss',initlnrho2='nothing'
  character (len=4) :: iinit_str
  complex :: coeflnrho=0.

  namelist /density_init_pars/ &
       cs0,rho0,ampllnrho,gamma,initlnrho,initlnrho2,widthlnrho, &
       rho_left,rho_right,lnrho_const,cs2bot,cs2top, &
       radius_lnrho,eps_planet, &
       b_ell,q_ell,hh0,rbound, &
       mpoly,strati_type, &
       kx_lnrho,ky_lnrho,kz_lnrho,amplrho,coeflnrho, &
       dlnrhobdx,co1_ss,co2_ss,Sigma1,cp

  namelist /density_run_pars/ &
       cs0,rho0,gamma,cdiffrho,diffrho,diffrho_shock,gradlnrho0, &
       cs2bot,cs2top,lupw_lnrho,cp
  ! diagnostic variables (needs to be consistent with reset list below)
  integer :: i_ekin=0,i_rhom=0,i_ekintot=0
  integer :: i_lnrhomphi=0,i_rhomphi=0

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
        print*, 'register_density: nvar = ', nvar
        print*, 'register_density: ilnrho = ', ilnrho
      endif
!
!  Put variable name in array
!
      varname(ilnrho) = 'lnrho'
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: density.f90,v 1.161 2004-05-18 15:38:24 dobler Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_density: nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',lnrho $'
          if (nvar == mvar) write(4,*) ',lnrho'
        else
          write(4,*) ',lnrho $'
        endif
        write(15,*) 'lnrho = fltarr(mx,my,mz)*one'
      endif
!
    endsubroutine register_density
!***********************************************************************
    subroutine initialize_density()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  For compatibility with other applications, we keep the possibility
!  of giving diffrho units of dxmin*cs0, but cs0 is not well defined general
!
!  24-nov-02/tony: coded 
!  31-aug-03/axel: normally, diffrho should be given in absolute units
!
    if(diffrho==0.) then
      diffrho=cdiffrho*dxmin*cs0
    endif
!
!   make sure all relevant parameters are set for spherical shell problems
!
    select case(initlnrho(1))
      case('geo-kws')
        if (lroot) print*,'initialize_density: set reference sound speed for spherical shell problem'
!ajwm Shouldn't be done here they should be set as input parameters in
!ajwm start.in.  There may be a problem at present with run since these
!ajwm vales will not get set consistantly.
!ajwm cs0 and lnrho0 are parameters of the EOS not of density
!        cs20=gamma
!        cs0=sqrt(cs20)
!        cp=5./2.
!ajwm Routine should really be called initialize_eos eventually
!ajwm important thing is the parameters of the Equation of State 
!ajwm have changed so we'd better recalculate everything consistently
!ajwm but that won't work either since 
!        call initialize_ionization()
    endselect
!
    endsubroutine initialize_density
!***********************************************************************
    subroutine init_lnrho(f,xx,yy,zz)
!
!  initialise lnrho; called from start.f90
!
!  7-nov-01/wolf: coded
! 28-jun-02/axel: added isothermal
! 15-oct-03/dave: added spherical shell (kws)
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
      real :: pot_ext,lnrho_ext,cs2_ext
      real :: zbot,ztop
      logical :: lnothing
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

      lnothing=.true.

      do iinit=1,ninit

      if (initlnrho(iinit)/='nothing') then

      lnothing=.false.

      call chn(iinit,iinit_str)

      select case(initlnrho(iinit))

      case('zero', '0'); f(:,:,:,ilnrho)=0.
      case('const_lnrho'); f(:,:,:,ilnrho)=lnrho_const
      case('constant'); f(:,:,:,ilnrho)=alog(rho_left)
      case('mode'); call modes(ampllnrho,coeflnrho,f,ilnrho,kx_lnrho,ky_lnrho,kz_lnrho,xx,yy,zz)
      case('blob'); call blob(ampllnrho,f,ilnrho,radius_lnrho,0.,0.,0.)
      case('isothermal'); call isothermal_density(f)
      case('stratification'); call stratification(f,strati_type)
      case('polytropic_simple'); call polytropic_simple(f)
      case('hydrostatic-z', '1'); print*, &
                              'init_lnrho: use polytropic_simple instead!'
      case('xjump'); call jump(f,ilnrho,lnrho_left,lnrho_right,widthlnrho,'x')
      case('yjump'); call jump(f,ilnrho,lnrho_left,lnrho_right,widthlnrho,'y')
      case('zjump'); call jump(f,ilnrho,lnrho_left,lnrho_right,widthlnrho,'z')
      case('soundwave-x'); call soundwave(ampllnrho,f,ilnrho,kx=1.)
      case('soundwave-y'); call soundwave(ampllnrho,f,ilnrho,ky=1.)
      case('soundwave-z'); call soundwave(ampllnrho,f,ilnrho,kz=1.)
      case('sinwave-x'); call sinwave(ampllnrho,f,ilnrho,kx=kx_lnrho)
      case('sinwave-y'); call sinwave(ampllnrho,f,ilnrho,ky=ky_lnrho)
      case('sinwave-z'); call sinwave(ampllnrho,f,ilnrho,kz=kz_lnrho)
      case('coswave-x'); call coswave(ampllnrho,f,ilnrho,kx=kx_lnrho)
      case('coswave-y'); call coswave(ampllnrho,f,ilnrho,ky=ky_lnrho)
      case('coswave-z'); call coswave(ampllnrho,f,ilnrho,kz=kz_lnrho)
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
        if (lroot) print*, &
               'init_lnrho: density jump; rho_left,right=',rho_left,rho_right
        if (lroot) print*, &
               'init_lnrho: density jump; widthlnrho=',widthlnrho
        prof=.5*(1.+tanh(zz/widthlnrho))
        f(:,:,:,ilnrho)=alog(rho_left)+alog(rho_left/rho_right)*prof

      case ('hydrostatic-z-2', '3')
        !
        !  hydrostatic density stratification for isentropic atmosphere
        !
        if (lgravz) then
          if (lroot) print*,'init_lnrho: vertical density stratification'
          !        f(:,:,:,ilnrho)=-zz
          ! isentropic case:
          !        zmax = -cs20/gamma1/gravz
          !        print*, 'zmax = ', zmax
          !        f(:,:,:,ilnrho) = 1./gamma1 * alog(abs(1-zz/zmax))
          ! linear entropy gradient;
          f(:,:,:,ilnrho) = -grads0*zz &
                            + 1./gamma1*alog( 1 + gamma1*gravz/grads0/cs20 &
                                                  *(1-exp(-grads0*zz)) )
        endif

      case ('hydrostatic-r')
        !
        !  hydrostatic radial density stratification for isentropic (or
        !  isothermal) sphere
        !
        if (lgravr) then
          if (lroot) print*, &
               'init_lnrho: radial density stratification (assumes s=const)'
!ajwm - here's the init call that needs sorting!
!          call initialize_gravity()     ! get coefficients cpot(1:5)

          call potential(xx,yy,zz,pot,POT0=pot0) ! gravity potential
          call output(trim(directory)//'/pot.dat',pot,1)
          !
          ! rho0, cs0, pot0 are the values in the centre
          !
          if (gamma /= 1) then  ! isentropic
            f(:,:,:,ilnrho) = lnrho0 &
                              + alog(1 - gamma1*(pot-pot0)/cs20) / gamma1
          else                  ! isothermal
            f(:,:,:,ilnrho) = lnrho0 - (pot-pot0)/cs20
          endif
        endif

      case ('isentropic-star')
        !
        !  isentropic/isothermal hydrostatic sphere"
        !    ss  = 0       for r<R,
        !    cs2 = const   for r>R
        !
        !  Only makes sense if both initlnrho=initss='isentropic-star'
        !
        if (lgravr) then
          if (lroot) print*, &
               'init_lnrho: isentropic star with isothermal atmosphere'
!          call initialize_gravity()     ! get coefficients cpot(1:5)
          call potential(xx,yy,zz,pot,POT0=pot0) ! gravity potential
          call output(trim(directory)//'/pot.dat',pot,1)
          !
          ! rho0, cs0, pot0 are the values in the centre
          !
          if (gamma /= 1) then
            ! Note:
            ! (a) `where' is expensive, but this is only done at
            !     initialization.
            ! (b) Comparing pot with pot_ext instead of r with r_ext will
            !     only work if grav_r<=0 everywhere -- but that seems
            !     reasonable.
            call potential(R=r_ext,POT=pot_ext) ! get pot_ext=pot(r_ext)
            lnrho_ext = lnrho0 + alog(1 - gamma1*(pot_ext-pot0)/cs20) / gamma1
            cs2_ext   = cs20*(1 - gamma1*(pot_ext-pot0)/cs20)
            ! Adjust for given cs2cool (if given) or set cs2cool (otherwise)
            if (cs2cool/=0) then
              lnrho_ext = lnrho_ext - alog(cs2cool/cs2_ext)
            else
              cs2cool   = cs2_ext
            endif
            ! Add temperature and entropy jump (such that pressure
            ! remains continuous) if cs2cool was specified in start.in:
            ! where (sqrt(xx**2+yy**2+zz**2) <= r_ext) ! isentropic for r<r_ext
            where (pot <= pot_ext) ! isentropic for r<r_ext
              f(:,:,:,ilnrho) = lnrho0 &
                                + alog(1 - gamma1*(pot-pot0)/cs20) / gamma1
            elsewhere           ! isothermal for r>r_ext
              f(:,:,:,ilnrho) = lnrho_ext - gamma*(pot-pot_ext)/cs2cool
            endwhere
          else                  ! gamma=1 --> simply isothermal (I guess [wd])
            f(:,:,:,ilnrho) = lnrho0 - (pot-pot0)/cs20
          endif
        endif

      case ('piecew-poly', '4')
        !
        !  piecewise polytropic for stellar convection models
        !
        if (lroot) print*, &
           'init_lnrho: piecewise polytropic vertical stratification (lnrho)'
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

      case ('piecew-disc', '41')
        !
        !  piecewise polytropic for accretion discs
        !
        if (lroot) print*, &
             'init_lnrho: piecewise polytropic disc stratification (lnrho)'
        ! bottom region

        cs2int = cs0**2
        lnrhoint = lnrho0
        f(:,:,:,ilnrho) = lnrho0 ! just in case
        call polytropic_lnrho_disc(f,mpoly1,zz,tmp,zref,z1,z1, &
                            0,cs2int,lnrhoint)
          ! unstable layer
        call polytropic_lnrho_disc(f,mpoly0,zz,tmp,z1,z2,z2,0,cs2int,lnrhoint)
          ! stable layer (top)
        call polytropic_lnrho_disc(f,mpoly2,zz,tmp,z2,ztop,ztop, &
                                   isothtop,cs2int,lnrhoint)
        !
        !  calculate cs2bot and cs2top for run.x (boundary conditions)
        !
        !cs2bot = cs2int + gamma/gamma1*gravz*nu_epicycle**2/(mpoly2+1)* &
        !         (zbot**2-z0**2)/2.
        cs2bot = cs20
        if (isothtop /= 0) then
          cs2top = cs20
        else
          cs2top = cs20 + gamma/gamma1*gravz*nu_epicycle**2/(mpoly0+1)* &
                  (ztop**2-zref**2)/2.
        endif

      case ('polytropic', '5')
        !
        !  polytropic stratification
        !  cs0, rho0 and ss0=0 refer to height z=zref
        !
        if (lroot) print*, &
                  'init_lnrho: polytropic vertical stratification (lnrho)'
        !
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
        if (lroot) print*,'init_lnrho: x-wave in lnrho; ampllnrho=',ampllnrho
        f(:,:,:,ilnrho)=lnrho_const+ampllnrho*sin(kx_lnrho*xx)

      case('sound-wave-exp')
        !
        !  sound wave (should be consistent with hydro module)
        !
        if (lroot) print*,'init_lnrho: x-wave in rho; ampllnrho=',ampllnrho
        f(:,:,:,ilnrho)=alog(rho_const+amplrho*sin(kx_lnrho*xx))

      case('sound-wave2')
        !
        !  sound wave (should be consistent with hydro module)
        !
        if (lroot) print*,'init_lnrho: x-wave in lnrho; ampllnrho=',ampllnrho
        f(:,:,:,ilnrho)=lnrho_const+ampllnrho*cos(kx_lnrho*xx)

      case('shock-tube', '13')
        !
        !  shock tube test (should be consistent with hydro module)
        !  
        if (lroot) print*,'init_lnrho: polytopic standing shock'
        prof=.5*(1.+tanh(xx/widthlnrho))
        f(:,:,:,ilnrho)=alog(rho_left)+(alog(rho_right)-alog(rho_left))*prof

      case('sin-xy')
        !
        !  sin profile in x and y
        !  
        if (lroot) print*,'init_lnrho: lnrho=sin(x)*sin(y)'
        f(:,:,:,ilnrho) = &
             alog(rho0) + ampllnrho*sin(kx_lnrho*xx)*sin(ky_lnrho*yy)

      case('sin-xy-rho')
        !
        !  sin profile in x and y, but in rho, not ln(rho)
        !  
        if (lroot) print*,'init_lnrho: rho=sin(x)*sin(y)'
        f(:,:,:,ilnrho) = &
             alog(rho0*(1+ampllnrho*sin(kx_lnrho*xx)*sin(ky_lnrho*yy)))

      case('linear')
        !
        !  linear profile in kk.xxx
        !  
        if (lroot) print*,'init_lnrho: linear profile'
        f(:,:,:,ilnrho) = alog(rho0) &
             + ampllnrho*(kx_lnrho*xx+ky_lnrho*yy+kz_lnrho*zz) &
               / sqrt(kx_lnrho**2+ky_lnrho**2+kz_lnrho**2)

      case('planet')
        !
        !  planet solution of Goodman, Narayan & Goldreich (1987)
        !  (Simple 3-D)
        !
        call planet(rbound,f,xx,yy,zz,eps_planet,radius_lnrho, &
            gamma,cs20,rho0,widthlnrho,hh0)

      case('planet_hc')
        !
        !  planet solution of Goodman, Narayan & Goldreich (1987)
        !  (3-D with hot corona)
        !
        call planet_hc(amplrho,f,xx,yy,zz,eps_planet,radius_lnrho, &
            gamma,cs20,rho0,widthlnrho)

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


      case('baroclinic')
        !
        !  Baroclinic initial condition
        !
        call baroclinic(f,xx,yy,zz,gamma,rho0,dlnrhobdx,co1_ss,co2_ss,cs20)

      
      case('Ferriere'); if(lroot) print*,'init_lnrho: Ferriere set in entropy'

      case('geo-kws')
      !
      ! radial hydrostatic profile in shell region only
      !
        if (lroot) print*,'init_lnrho: kws hydrostatic in spherical shell region'
        call shell_lnrho(f)

      case('geo-kws-constant-T')
      !
      ! radial hydrostatic profile throughout box, which is consistent
      ! with constant temperature in exterior regions, and gives continuous
      ! density at shell boundaries 
      !
        if (lroot) print*,'init_lnrho: kws hydrostatic in spherical shell and exterior'
        call shell_lnrho(f)

      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*,'init_lnrho: No such value for initlnrho(' &
                          //trim(iinit_str)//'): ',trim(initlnrho(iinit))
        call stop_it('')

      endselect

      if (lroot) print*,'init_lnrho: initlnrho(' &
                        //trim(iinit_str)//') = ',trim(initlnrho(iinit))

      endif

      enddo

      if (lnothing.and.lroot) print*,'init_lnrho: zero density'
!
!  sanity check
!
      if (notanumber(f(l1:l2,m1:m2,n1:n2,ilnrho))) then
        STOP "init_lnrho: Imaginary density values"
      endif
!
!  check that cs2bot,cs2top are ok
!  for runs with ionization or fixed ionization, don't print them
!
      if (lionization .or. lionization_fixed) then
        cs2top=impossible
        cs2bot=impossible
      else
        if(lroot) print*,'init_lnrho: cs2bot,cs2top=',cs2bot,cs2top
      endif
!
!  Add some structures to the lnrho initialized above
!  Tobi: This is only kept for backwards compatibility
!
      select case(initlnrho2)
    
      case('nothing')
    
        if (lroot.and.ip<=5) print*,"init_lnrho: initlnrho2='nothing'"

      case('addblob')

        if (lroot) print*,'init_lnrho: add blob'
        if (lroot) print*,'initlnrho2 is deprecated. '// &
                          "use initlnrho='initcond1','initcond2' instead"
        f(:,:,:,ilnrho)=f(:,:,:,ilnrho) &
          +ampllnrho*exp(-(xx**2+yy**2+zz**2)/radius_lnrho**2)

      case default

        if (lroot) print*,'init_lnrho: No such value for initlnrho2: ', &
                          trim(initlnrho2)
        call stop_it("")

      endselect
!
    endsubroutine init_lnrho
!***********************************************************************
    subroutine calculate_vars_rho(f,rho1)
!
!   Calculation of rho1
!
!   08-febr-04/bing: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f       
      real, dimension (nx) :: rho1

      intent(in) :: f
      intent(out) :: rho1
 
      rho1=exp(-f(l1:l2,m,n,ilnrho))  
      
    endsubroutine calculate_vars_rho
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
!  lnrhoin -- value of lnrho at the interface, i.e. at the zint on entry,
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
    subroutine polytropic_lnrho_disc( &
         f,mpoly,zz,tmp,zint,zbot,zblend,isoth,cs2int,lnrhoint)
!
!  Implement a polytropic profile in a disc. If this routine is
!  called several times (for a piecewise polytropic atmosphere), on needs
!  to call it from bottom (middle of disc) upwards.
!
!  zint    -- height of (previous) interface, where cs2int and lnrhoint
!             are set
!  zbot    -- z at top of layer (name analogous with polytropic_lnrho_z)
!  zblend  -- smoothly blend (with width whcond) previous ss (for z>zblend)
!             with new profile (for z<zblend)
!  isoth   -- flag for isothermal stratification;
!  lnrhoint -- value of lnrho at the interface, i.e. at the zint on entry,
!             at the zbot on exit
!  cs2int  -- same for cs2
!
!  24-jun-03/ulf:  coded
!
      use Sub, only: step
      use Gravity, only: gravz, nu_epicycle
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: tmp,p,zz
      real, dimension (mz) :: stp
      real :: mpoly,zint,zbot,zblend,beta1,cs2int,lnrhoint,nu_epicycle2
      integer :: isoth
!
      ! NB: beta1 is not dT/dz, but dcs2/dz = (gamma-1)c_p dT/dz
      nu_epicycle2 = nu_epicycle**2
      if (isoth /= 0) then ! isothermal layer
        beta1 = 0.
        tmp = gamma*gravz*nu_epicycle2/cs2int*(zz**2-zint**2)/2.
        lnrhoint = lnrhoint + gamma*gravz*nu_epicycle2/cs2int* &
                   (zbot**2-zint**2)/2.
      else
        beta1 = gamma*gravz*nu_epicycle2/(mpoly+1)
        tmp = 1 + beta1*(zz**2-zint**2)/cs2int/2.
        tmp = max(tmp,epsi)  ! ensure arg to log is positive
        tmp = lnrhoint + mpoly*alog(tmp)
        lnrhoint = lnrhoint + mpoly*alog(1 + beta1*(zbot**2-zint**2)/cs2int/2.)
      endif
      cs2int = cs2int + beta1*(zbot**2-zint**2)/2. ! cs2 at layer interface (bottom)
      !
      ! smoothly blend the old value (above zblend) and the new one (below
      ! zblend) for the two regions:
      !
      stp = step(z,zblend,widthlnrho)
      p = spread(spread(stp,1,mx),2,my)
      f(:,:,:,ilnrho) = p*f(:,:,:,ilnrho) + (1-p)*tmp
!
    endsubroutine polytropic_lnrho_disc
!***********************************************************************
    subroutine shell_lnrho(f)
!
!  Initialize density based on specified radial profile in
!  a spherical shell
!
!  22-oct-03/dave -- coded
!
      use Gravity, only: g0,r0_pot,n_pot,smoothpotential
      use Sub, only: calc_unitvects_sphere
!
      real, dimension (mx,my,mz,mvar+maux), intent(inout) :: f
      real, dimension (nx) :: pot
      real :: beta1,lnrho_int,lnrho_ext,c_int,c_ext
!
      beta1 = g0/(mpoly+1)
!
!     densities at shell boundaries
      lnrho_int = mpoly*log(1+beta1*(1/r_int-1))
      lnrho_ext = lnrho0
!
!     constants for continuity of density across boundaries
!     when using smoothed potential
      c_int = lnrho_int-(g0*exp(-mpoly*lnrho_int))*(r_int**n_pot+r0_pot**n_pot)**(-1./float(n_pot))
      c_ext = lnrho_ext-(g0*exp(-mpoly*lnrho_ext))*(r_ext**n_pot+r0_pot**n_pot)**(-1./float(n_pot))
      
        do imn=1,ny*nz
          n=nn(imn)
          m=mm(imn)
!
          call calc_unitvects_sphere()
!
          ! in the fluid shell
          where (r_mn < r_ext .AND. r_mn > r_int) f(l1:l2,m,n,ilnrho) = mpoly*log(1+beta1*(1/r_mn-1))
          ! outside the fluid shell
          if (initlnrho(1)=='geo-kws') then
            where (r_mn >= r_ext) f(l1:l2,m,n,ilnrho) = lnrho_ext
            where (r_mn <= r_int) f(l1:l2,m,n,ilnrho) = lnrho_int
          elseif (initlnrho(1)=='geo-kws-constant-T') then
            call smoothpotential(RMN=r_mn,POT=pot)
            where (r_mn >= r_ext) f(l1:l2,m,n,ilnrho) = pot+c_ext
            where (r_mn <= r_int) f(l1:l2,m,n,ilnrho) = pot+c_int
          endif
        enddo 
!      
    endsubroutine shell_lnrho
!***********************************************************************
    subroutine dlnrho_dt(f,df,uu,glnrho,divu,lnrho,shock,gshock)
!
!  continuity equation
!  calculate dlnrho/dt = - u.gradlnrho - divu
!
!   7-jun-02/axel: incoporated from subroutine pde
!
      use Sub
      use Special, only: special_calc_density
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: uu,glnrho,gshock
      real, dimension (nx) :: lnrho,divu,uglnrho,gshockglnrho,glnrho2,shock
      real, dimension (nx) :: del2lnrho
      integer :: j
!
      intent(in)  :: f,uu,divu,shock,gshock
      intent(out) :: df,glnrho,lnrho
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dlnrho_dt: SOLVE dlnrho_dt'
      if (headtt) call identify_bcs('lnrho',ilnrho)
!
!  define lnrho; calculate density gradient and avection term
!
      lnrho=f(l1:l2,m,n,ilnrho)
      call grad(f,ilnrho,glnrho)
      call u_dot_gradf(f,ilnrho,glnrho,uu,uglnrho,UPWIND=lupw_lnrho)
!
!  continuity equation
!
      df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)-uglnrho-divu
!
!  mass diffusion and shock diffusion, in absolute units (similar to nu, chi, and eta)
!
      if ((diffrho/=0.) .or. (diffrho_shock/=0.)) then
      
        call del2(f,ilnrho,del2lnrho)
        call dot2_mn(glnrho,glnrho2)
     
        if (diffrho/=0.) then
          if(headtt) print*,'dlnrho_dt: diffrho=',diffrho
          df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)+diffrho*(del2lnrho+glnrho2)
          call max_for_dt(diffrho,maxdiffus)
        endif
!
        if (diffrho_shock/=0.) then
          call dot_mn(gshock,glnrho,gshockglnrho)
          if(headtt) print*,'dlnrho_dt: diffrho_shock=',diffrho_shock
          df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)+diffrho_shock*shock*(del2lnrho+glnrho2)+diffrho_shock*gshockglnrho
          call max_for_dt(diffrho_shock*maxval(shock),maxdiffus)
        endif
      endif
!
!  add advection of imposed constant gradient of lnrho (called gradlnrho0)
!  makes sense really only for periodic boundary conditions
!  This gradient can have arbitary direction.
!
      do j=1,3
        if (gradlnrho0(j)/=0.) then
          df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)-gradlnrho0(j)*uu(:,j)
        endif
      enddo

      if (lspecial) call special_calc_density(f,df,uu,glnrho,divu,lnrho)
!
!  phi-averages
!  Note that this does not necessarily happen with ldiagnos=.true.
!
      if (l2davgfirst) then
        if (i_lnrhomphi/=0) call phisum_mn_name_rz(lnrho,i_lnrhomphi)
        if (i_rhomphi/=0) call phisum_mn_name_rz(exp(lnrho),i_rhomphi)
      endif
!
    endsubroutine dlnrho_dt
!***********************************************************************
    subroutine rprint_density(lreset,lwrite)
!
!  reads and registers print parameters relevant for compressible part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Sub
!
      integer :: iname,irz
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
        i_ekin=0; i_rhom=0; i_ekintot=0
        i_lnrhomphi=0; i_rhomphi=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'rprint_density: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ekintot',i_ekintot)
        call parse_name(iname,cname(iname),cform(iname),'ekin',i_ekin)
        call parse_name(iname,cname(iname),cform(iname),'rhom',i_rhom)
      enddo
!
!  check for those quantities for which we want phi-averages
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'lnrhomphi',i_lnrhomphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'rhomphi',i_rhomphi)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_ekintot=',i_ekintot
        write(3,*) 'i_ekin=',i_ekin
        write(3,*) 'i_rhom=',i_rhom
        write(3,*) 'nname=',nname
        write(3,*) 'ilnrho=',ilnrho
        write(3,*) 'i_lnrhomphi=',i_lnrhomphi
        write(3,*) 'i_rhomphi=',i_rhomphi
      endif
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
!      use Ionization, only: lionization_fixed, isothermal_density_ion 
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) :: pot,tmp
!
!  Stratification depends on the gravity potential
!
      if (lroot) print*,'isothermal_density: isothermal stratification'
      do n=n1,n2
        do m=m1,m2
          call potential(x(l1:l2),y(m),z(n),pot)
!        if (lionization_fixed) then
!           call isothermal_density_ion(pot,tmp)
!            tmp=0.
!        else
          tmp=-gamma*pot/cs20
!        endif
          f(l1:l2,m,n,ilnrho)=lnrho0+tmp
!        if(lentropy) f(l1:l2,m,n,iss)= -gamma1/gamma*tmp
!                                      = gamma1*pot/cs20
!   MOVED to isothermal_entropy routine

          if(lentropy) f(l1:l2,m,n,iss)= &
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
        call stop_it( &
                 'polytropic_simply: rho and cs2 will vanish within domain')
      endif
!
      ggamma=1.+1./mpoly
!
      do n=n1,n2
      do m=m1,m2
        call potential(x(l1:l2),y(m),z(n),pot)
        dlncs2=alog(-gamma*pot/((mpoly+1.)*cs20))
        f(l1:l2,m,n,ilnrho)=lnrho0+mpoly*dlncs2
        if(lentropy) f(l1:l2,m,n,iss)=mpoly*(ggamma/gamma-1.)*dlncs2
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
      if(ldebug) print*,'bc_lnrho_temp_z: cs20,cs0=',cs20,cs0
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
        if (ldebug) print*, &
                 'bc_lnrho_temp_z: set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0. .and. lroot) print*, &
                 'bc_lnrho_temp_z: cannot have cs2bot<=0'
        tmp = 2/gamma*alog(cs2bot/cs20)
!
!  set boundary value for entropy, then extrapolate ghost pts by antisymmetry
!
        f(:,:,n1,iss) = 0.5*tmp - gamma1/gamma*(f(:,:,n1,ilnrho)-lnrho0)
        do i=1,nghost; f(:,:,n1-i,iss) = 2*f(:,:,n1,iss)-f(:,:,n1+i,iss); enddo
!
!  set density in the ghost zones so that dlnrho/dz + ds/dz = gz/cs2bot
!  for the time being, we don't worry about lnrho0 (assuming that it is 0)
!
        tmp=-gravz/cs2bot
        do i=1,nghost
          f(:,:,n1-i,ilnrho) = f(:,:,n1+i,ilnrho) +f(:,:,n1+i,iss) &
                                                  -f(:,:,n1-i,iss) +2*i*dz*tmp
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                    'bc_lnrho_temp_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0. .and. lroot) print*, &
                    'bc_lnrho_temp_z: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
!
!  set boundary value for entropy, then extrapolate ghost pts by antisymmetry
!
        f(:,:,n2,iss) = 0.5*tmp - gamma1/gamma*(f(:,:,n2,ilnrho)-lnrho0)
        do i=1,nghost; f(:,:,n2+i,iss) = 2*f(:,:,n2,iss)-f(:,:,n2-i,iss); enddo
!
!  set density in the ghost zones so that dlnrho/dz + ds/dz = gz/cs2top
!  for the time being, we don't worry about lnrho0 (assuming that it is 0)
!
        tmp=gravz/cs2top
        do i=1,nghost
          f(:,:,n2+i,ilnrho) = f(:,:,n2-i,ilnrho) +f(:,:,n2-i,iss) &
                                                  -f(:,:,n2+i,iss) +2*i*dz*tmp
        enddo

      case default
        if(lroot) print*,"bc_lnrho_temp_z: invalid argument"
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
      integer :: i
!
      if(ldebug) print*,'bc_lnrho_pressure_z: cs20,cs0=',cs20,cs0
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
!             f(l,m,n1,iss)=ss_bot
!           else
!             f(l,m,n1,iss)=f(l,m,n1+1,iss)
!           endif
!         enddo
!         enddo
          f(:,:,n2,iss)=ss_top
          do i=1,nghost; f(:,:,n2+i,iss)=2*f(:,:,n2,iss)-f(:,:,n2-i,iss); enddo
!
!  set density value such that pressure is constant at the bottom
!
          f(:,:,n2,ilnrho)=lnrho_top+ss_top-f(:,:,n2,iss)
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
!             f(l,m,n1,iss)=ss_bot
!           else
!             f(l,m,n1,iss)=f(l,m,n1+1,iss)
!           endif
!         enddo
!         enddo
          f(:,:,n1,iss)=ss_bot
          do i=1,nghost; f(:,:,n1-i,iss)=2*f(:,:,n1,iss)-f(:,:,n1+i,iss); enddo
!
!  set density value such that pressure is constant at the bottom
!
          f(:,:,n1,ilnrho)=lnrho_bot+ss_bot-f(:,:,n1,iss)
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
        if(lroot) print*,"bc_lnrho_pressure_z: invalid argument"
        call stop_it("")
      endselect
!
    endsubroutine bc_lnrho_pressure_z
!***********************************************************************

endmodule Density
