! $Id: entropy_fixed.f90,v 1.3 2004-06-11 08:07:35 ajohan Exp $

!  This module takes care of entropy (initial condition
!  and time advance)

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 1
!
!***************************************************************

module Entropy

  use Cparam
  use Cdata
  use Hydro
  use Interstellar
  use Viscosity
  use Density, only: lcalc_cp,cs2cool

  implicit none

!AB: need to avoid Kappa0, because it collides with kappa0 in ionization

  !real, dimension (nx) :: cs2,TT1
  real :: radius_ss=0.1,ampl_ss=0.,widthss=2*epsi,epsilon_ss
  real :: luminosity=0.,wheat=0.1,cool=0.,rcool=1.,wcool=0.1
  real :: TT_int,TT_ext,cs2_int,cs2_ext,cool_int=0.,cool_ext=0.
  real :: chi=0.,chi_t=0.,chi_shock=0.,Kappa0=0.,Kisotr=0.
  real :: Kgperp=0.,Kgpara=0.
  real :: ss_left=1.,ss_right=1.
  real :: ss0=0.,khor_ss=1.,ss_const=0.
  real :: tau_ss_exterior=0.,T0=1.
  !parameters for Sedov type initial condition
  real :: center1_x=0., center1_y=0., center1_z=0.
  real :: center2_x=0., center2_y=0., center2_z=0.
  real :: kx_ss
  real :: thermal_background=0., thermal_peak=0., thermal_scaling=1.
  real :: hcond0=impossible
  real :: hcond1=impossible,hcond2=impossible
  real :: Fbot=impossible,FbotKbot=impossible,Kbot=impossible
  real :: Ftop=impossible,FtopKtop=impossible,Ktop=impossible
  real :: tau_cor=0.,TT_cor=0.,z_cor=0.
  real :: tauheat_buffer=0.,TTheat_buffer=0.,zheat_buffer=0.,dheat_buffer1=0.
  real :: heat_uniform=0.
  real :: deltaT_poleq=0.
  logical :: lturbulent_heat=.false.
  logical :: lcalc_heatcond_simple=.false.,lmultilayer=.true.
  logical :: lcalc_heatcond=.false.,lcalc_heatcond_constchi=.false.
  logical :: lupw_ss=.false.
  logical :: lgaspressuregradient=.true.
  character (len=labellen), dimension(ninit) :: initss='nothing'
  character (len=labellen) :: pertss='zero'
  character (len=labellen) :: cooltype='Temp',iheatcond='K-const'

  ! input parameters
  namelist /entropy_init_pars/ &
       initss,pertss,grads0,radius_ss,ampl_ss,widthss,epsilon_ss, &
       ss_left,ss_right,ss_const,mpoly0,mpoly1,mpoly2,isothtop, &
       khor_ss,thermal_background,thermal_peak,thermal_scaling,cs2cool, &
       center1_x, center1_y, center1_z, center2_x, center2_y, center2_z, &
       T0, kx_ss

  ! run parameters
  namelist /entropy_run_pars/ &
       hcond0,hcond1,hcond2,widthss, &
       luminosity,wheat,cooltype,cool,cs2cool,rcool,wcool,Fbot, &
       chi_t,chi_shock,chi,Kappa0,Kisotr,iheatcond, &
       Kgperp,Kgpara, &
       lcalc_heatcond_simple,lcalc_heatcond,lcalc_heatcond_constchi,&
       tau_ss_exterior,lmultilayer,Kbot,tau_cor,TT_cor,z_cor, &
       tauheat_buffer,TTheat_buffer,zheat_buffer,dheat_buffer1, &
       heat_uniform,lupw_ss,lcalc_cp,cool_int,cool_ext, &
       lturbulent_heat, deltaT_poleq, lgaspressuregradient

  ! other variables (needs to be consistent with reset list below)
  integer :: i_dtc=0,i_eth=0,i_ethdivum=0,i_ssm=0,i_ugradpm=0, i_ethtot=0
  integer :: i_dtchi=0
  integer :: i_ssmphi=0

  contains

!***********************************************************************
    subroutine register_entropy()
!
!  initialise variables which should know that we solve an entropy
!  equation: iss, etc; increase nvar accordingly
!
!  6-nov-01/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_ent called twice')
      first = .false.
!
      lentropy_fixed = .true.
!
      iss = mvar+naux+1             ! index to access entropy
      naux = naux+1
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_entropy: naux = ', naux
        print*, 'register_entropy: iss = ', iss
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: entropy_fixed.f90,v 1.3 2004-06-11 08:07:35 ajohan Exp $")
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
        call stop_it('register_entropy: naux > maux')
      endif
!
!  Put variable name in array
!
      varname(iss) = 'ss'
!
!  Writing files for use with IDL
!
      if (naux < maux) aux_var(aux_count)=',ss $'
      if (naux == maux) aux_var(aux_count)=',ss'
      aux_count=aux_count+1
      if (lroot) write(15,*) 'ss = fltarr(mx,my,mz)*one'
!
    endsubroutine register_entropy
!***********************************************************************
    subroutine initialize_entropy(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  21-jul-2002/wolf: coded
!  11-jun-04/anders: call init_ss again for entropy_fixed with not lstarting
!
      use Cdata
      use Gravity, only: gravz,g0
      use Density, only: mpoly
      use Ionization, only: lnTT0,get_soundspeed
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: beta1
      logical :: lstarting
!
      lneed_sij = .true.   !let Hydro module know to precalculate some things
      lneed_glnrho = .true.
!
!  radiative diffusion: initialize flux etc
!
      !
      !  Kbot and hcond0 are used interchangibly, so if one is 
      !  =impossible, set it to the other's value
      !
      if (hcond0 == impossible) then
        if (Kbot == impossible) then
          hcond0 = 0.
          Kbot = 0.
        else                    ! Kbot = possible
          hcond0 = Kbot
        endif
      else                      ! hcond0 = possible
        if (Kbot == impossible) then
          Kbot = hcond0
        else
          if (lroot) print*, 'WARNING: You should not set Kbot and hcond0 at the same time'
        endif
      endif
      !
      if (lgravz) then
        if (lmultilayer) then
          !
          !  calculate hcond1,hcond2 if they have not been set in run.in
          !
          if (hcond1==impossible) hcond1 = (mpoly1+1.)/(mpoly0+1.)
          if (hcond2==impossible) hcond2 = (mpoly2+1.)/(mpoly0+1.)
          !
          !  calculate Fbot if it has not been set in run.in
          !
          if (Fbot==impossible) then
            if (bcz1(iss)=='c1') then
              Fbot=-gamma/(gamma-1)*hcond0*gravz/(mpoly0+1)
              if (lroot) print*, &
                      'initialize_entropy: Calculated Fbot = ', Fbot
            else
              Fbot=0.
            endif
          endif
          if (hcond0*hcond1 /= 0.) then
            FbotKbot=Fbot/(hcond0*hcond1)
          else
            FbotKbot=0.
          endif
          !
          !  calculate Ftop if it has not been set in run.in
          !
          if (Ftop==impossible) then
            if (bcz2(iss)=='c1') then
              Ftop=-gamma/(gamma-1)*hcond0*gravz/(mpoly0+1)
              if (lroot) print*, &
                      'initialize_entropy: Calculated Ftop = ',Ftop
            else
              Ftop=0.
            endif
          endif
          if (hcond0*hcond2 /= 0.) then
            FtopKtop=Ftop/(hcond0*hcond2)
          else
            FtopKtop=0.
          endif
!
        else
          !
          !  Wolfgang, in future we should define chiz=chi(z) or Kz=K(z) here.
          !  calculate hcond and FbotKbot=Fbot/K
          !  (K=hcond is radiative conductivity)
          !
          !  calculate Fbot if it has not been set in run.in
          !
          if (Fbot==impossible) then
            if (bcz1(iss)=='c1') then
              Fbot=-gamma/(gamma-1)*hcond0*gravz/(mpoly+1)
              if (lroot) print*, &
                   'initialize_entropy: Calculated Fbot = ', Fbot
              Kbot=gamma1/gamma*(mpoly+1.)*Fbot
              FbotKbot=gamma/gamma1/(mpoly+1.)
              if(lroot) print*,'initialize_entropy: Calculated Fbot,Kbot=', &
                   Fbot,Kbot
            ! else
            !! Don't need Fbot in this case (?)
            !  Fbot=-gamma/(gamma-1)*hcond0*gravz/(mpoly+1)
            !  if (lroot) print*, &
            !       'initialize_entropy: Calculated Fbot = ', Fbot
            endif
          endif
          !
          !  calculate Ftop if it has not been set in run.in
          !
          if (Ftop==impossible) then
            if (bcz2(iss)=='c1') then
              Ftop=-gamma/(gamma-1)*hcond0*gravz/(mpoly+1)
              if (lroot) print*, &
                      'initialize_entropy: Calculated Ftop = ', Ftop
              Ktop=gamma1/gamma*(mpoly+1.)*Ftop
              FtopKtop=gamma/gamma1/(mpoly+1.)
              if(lroot) print*,'initialize_entropy: Ftop,Ktop=',Ftop,Ktop
            ! else
            !! Don't need Ftop in this case (?)
            !  Ftop=0.
            endif
          endif
!
        endif
      endif
!
!   make sure all relevant parameters are set for spherical shell problems
!
      select case(initss(1))
        case('geo-kws')
          if (lroot) then
            print*,'initialize_entropy: set boundary temperatures for spherical shell problem'
            if (abs(exp(lnTT0)-T0) > epsi) then
              print*,'initialize_entropy: T0 is not consistent with cs20; using cs20'
              T0=exp(lnTT0)
            endif
          endif
!
!         temperatures at shell boundaries
          beta1=g0/(mpoly+1)
          TT_ext=T0
          TT_int=1+beta1*(1/r_int-1)
!         TT_ext=gamma/gamma1*T0
!         TT_int=gamma/gamma1*(1+beta1*(1/r_int-1))
!         set up cooling parameters for spherical shell in terms of
!         sound speeds
          call get_soundspeed(log(TT_ext),cs2_ext)
          call get_soundspeed(log(TT_int),cs2_int)
!
      endselect
!  
!  Rewrite initial condition into f for fixed entropy
!        
      if (.not. lstarting) call init_ss(f,spread(spread(x,2,my),3,mz), &
                  spread(spread(y,1,mx),3,mz),spread(spread(z,1,mx),2,my))
!
    endsubroutine initialize_entropy
!***********************************************************************
    subroutine init_ss(f,xx,yy,zz)
!
!  initialise entropy; called from start.f90
!  07-nov-2001/wolf: coded
!  24-nov-2002/tony: renamed for consistancy (i.e. init_[variable name]) 
!
      use Cdata
      use Mpicomm
      use Sub
      use Gravity
      use Initcond
      use Ionization
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz,tmp,pot
      real :: cs2int,ssint,ztop,ss_ext,cd2_ext,pot0,pot_ext
      logical :: lnothing=.true.
!
      intent(in) :: xx,yy,zz
      intent(inout) :: f
!
      do iinit=1,ninit
!
!  if we pretend that ss in in reality g1lnTT, we initialize the background
!  of lnTT/gamma such that it corresponds to ss=0.
!
      if (pretend_lnTT) f(:,:,:,iss)=(f(:,:,:,ilnrho)*gamma1-alog(gamma1))/gamma
!
      if (initss(iinit)/='nothing') then
!
      lnothing=.false.
      call chn(iinit,iinit_str)
!
!  select different initial conditions
!
      select case(initss(iinit))

        case('zero', '0'); f(:,:,:,iss) = 0.
        case('const_ss'); f(:,:,:,iss)=f(:,:,:,iss)+ss_const
        case('blob'); call blob(ampl_ss,f,iss,radius_ss,0.,0.,0.)
        case('isothermal'); call isothermal_entropy(f,T0)
        case('isothermal_lnrho_ss')
          print*, 'init_ss: Isothermal density and entropy stratification'
          call isothermal_lnrho_ss(f,T0,rho0)
        case('hydrostatic-isentropic')
          call hydrostatic_isentropic(f,lnrho_bot,ss_const)
        case('wave'); f(:,:,:,iss) = ampl_ss*sin(kx_ss*xx(:,:,:) + pi)
        case('Ferriere'); call ferriere(f) 
        case('xjump'); call jump(f,iss,ss_left,ss_right,widthss,'x')
        case('yjump'); call jump(f,iss,ss_left,ss_right,widthss,'y')
        case('zjump'); call jump(f,iss,ss_left,ss_right,widthss,'z')
        case('hor-fluxtube'); call htube(ampl_ss,f,iss,iss,xx,yy,zz,radius_ss,epsilon_ss)
        case('hor-tube'); call htube2(ampl_ss,f,iss,iss,xx,yy,zz,radius_ss,epsilon_ss)

        case('sedov') 
          if (lroot) print*,'init_ss: sedov - thermal background with gaussian energy burst'
        call blob(thermal_peak,f,iss,radius_ss,center1_x,center1_y,center1_z)
      !   f(:,:,:,iss) = f(:,:,:,iss) + (alog(f(:,:,:,iss) + thermal_background)+alog(thermal_scaling))/gamma 

        case('sedov-dual') 
          if (lroot) print*,'init_ss: sedov - thermal background with gaussian energy burst'
        call blob(thermal_peak,f,iss,radius_ss,center1_x,center1_y,center1_z)
        call blob(thermal_peak,f,iss,radius_ss,center2_x,center2_y,center2_z)
      !   f(:,:,:,iss) = (alog(f(:,:,:,iss) + thermal_background)+alog(thermal_scaling))/gamma 
   
        case('shock2d') 
          call shock2d(f,xx,yy,zz)

        case('isobaric')
          !
          !  ss = - ln(rho/rho0)
          !
          if (lroot) print*,'init_ss: isobaric stratification'
          f(:,:,:,iss) = -(f(:,:,:,ilnrho)-lnrho0)

        case('isentropic', '1')
          !
          !  ss = const.
          !
          if (lroot) print*,'init_ss: isentropic stratification'
          ! ss0=alog(-gamma1*gravz*zinfty)/gamma
          ! print*,'init_ss: isentropic stratification; ss=',ss0
          f(:,:,:,iss)=0.
          if (ampl_ss/=0.) then
            print*,'init_ss: put bubble: radius_ss,ampl_ss=',radius_ss,ampl_ss
            tmp=xx**2+yy**2+zz**2
            f(:,:,:,iss)=f(:,:,:,iss)+ampl_ss*exp(-tmp/amax1(radius_ss**2-tmp,1e-20))
          !f(:,:,:,iss)=f(:,:,:,iss)+ampl_ss*exp(-tmp/radius_ss**2)
          endif

        case('linprof', '2')
          !
          !  linear profile of ss, centered around ss=0.
          !
          if (lroot) print*,'init_ss: linear entropy profile'
          f(:,:,:,iss) = grads0*zz

        case('isentropic-star')
          !
          !  isentropic/isothermal hydrostatic sphere"
          !    ss  = 0       for r<R,
          !    cs2 = const   for r>R
          !
          !  Only makes sense if both initlnrho=initss='isentropic-star'
          !
          if (.not. ldensity) &
               call stop_it('isentropic-star requires density.f90')
          if (initlnrho(1) /= initss(1)) &
               call stop_it('isentropic star requires initlnrho=initss')
          if (lgravr) then
            if (lroot) print*, &
                 'init_lnrho: isentropic star with isothermal atmosphere'
            ! call initialize_gravity()     ! already done by init_lnrho
            call potential(xx,yy,zz,POT=pot,POT0=pot0) ! gravity potential
            !
            ! rho0, cs0,pot0 are the values in the centre
            !
            if (gamma /= 1) then
              ! Note:
              ! (a) `where' is expensive, but this is only done at
              !     initialization.
              ! (b) Comparing pot with pot_ext instead of r with r_ext will
              !     only work if grav_r<=0 everywhere -- but that seems
              !     reasonable.
              call potential(R=r_ext,POT=pot_ext) ! get pot_ext=pot(r_ext)
              cs2_ext   = cs20*(1 - gamma1*(pot_ext-pot0)/cs20)
              !
              ! Make sure init_lnrho (or start.in) has already set cs2cool:
              !
              if (cs2cool == 0) &
                   call stop_it("Inconsistency: cs2cool can't be 0")
              ss_ext = 0. + alog(cs2cool/cs2_ext)
              ! where (sqrt(xx**2+yy**2+zz**2) <= r_ext) ! isentropic f. r<r_ext
              where (pot <= pot_ext) ! isentropic for r<r_ext
                f(:,:,:,iss) = 0.
              elsewhere           ! isothermal for r>r_ext
                f(:,:,:,iss) = ss_ext + gamma1*(pot-pot_ext)/cs2cool
              endwhere
            else                  ! gamma=1 --> simply isothermal (I guess [wd])
              ! [NB: Never tested this..]
              f(:,:,:,iss) = -gamma1/gamma*(f(:,:,:,ilnrho)-lnrho0)
            endif
          endif

        case('piecew-poly', '4')
          !
          !  piecewise polytropic convection setup
          !  cs0, rho0 and ss0=0 refer to height z=zref
          !
          if (lroot) print*, &
                 'init_ss: piecewise polytropic vertical stratification (ss)'
          !
!         !  override hcond1,hcond2 according to polytropic equilibrium
!         !  solution
!         !
!         hcond1 = (mpoly1+1.)/(mpoly0+1.)
!         hcond2 = (mpoly2+1.)/(mpoly0+1.)
!         if (lroot) &
!              print*, &
!              'Note: mpoly{1,2} override hcond{1,2} to ', hcond1, hcond2
        !
          cs2int = cs0**2
          ss0 = 0.              ! reference value ss0 is zero
          ssint = ss0
          f(:,:,:,iss) = 0.    ! just in case
          ! top layer
          call polytropic_ss_z(f,mpoly2,zz,tmp,zref,z2,z0+2*Lz, &
                               isothtop,cs2int,ssint)
          ! unstable layer
          call polytropic_ss_z(f,mpoly0,zz,tmp,z2,z1,z2,0,cs2int,ssint)
          ! stable layer
          call polytropic_ss_z(f,mpoly1,zz,tmp,z1,z0,z1,0,cs2int,ssint)

        case('piecew-disc', '41')
          !
          !  piecewise polytropic convective disc
          !  cs0, rho0 and ss0=0 refer to height z=zref
          !
          if (lroot) print*,'init_ss: piecewise polytropic disc'
          !
!         !  override hcond1,hcond2 according to polytropic equilibrium
!         !  solution
!         !
!         hcond1 = (mpoly1+1.)/(mpoly0+1.)
!         hcond2 = (mpoly2+1.)/(mpoly0+1.)
!         if (lroot) &
!              print*, &
!        'init_ss: Note: mpoly{1,2} override hcond{1,2} to ', hcond1, hcond2
        !
          ztop = xyz0(3)+Lxyz(3)
          cs2int = cs0**2
          ss0 = 0.              ! reference value ss0 is zero
          ssint = ss0
          f(:,:,:,iss) = 0.    ! just in case
          ! bottom (middle) layer
          call polytropic_ss_disc(f,mpoly1,zz,tmp,zref,z1,z1, &
                               0,cs2int,ssint)
          ! unstable layer
          call polytropic_ss_disc(f,mpoly0,zz,tmp,z1,z2,z2,0,cs2int,ssint)
          ! stable layer (top)
          call polytropic_ss_disc(f,mpoly2,zz,tmp,z2,ztop,ztop,&
                               isothtop,cs2int,ssint)

        case('polytropic', '5')
          !
          !  polytropic stratification
          !  cs0, rho0 and ss0=0 refer to height z=zref
          !
          if (lroot) print*,'init_ss: polytropic vertical stratification'
          !
          cs20 = cs0**2
          ss0 = 0.              ! reference value ss0 is zero
          f(:,:,:,iss) = ss0   ! just in case
          cs2int = cs20
          ssint = ss0
          ! only one layer
          call polytropic_ss_z(f,mpoly0,zz,tmp,zref,z0,z0+2*Lz,0,cs2int,ssint)
          ! reset mpoly1, mpoly2 (unused) to make IDL routine `thermo.pro' work
          mpoly1 = mpoly0
          mpoly2 = mpoly0

        case ('geo-kws')
          !
          ! radial temperature profiles for spherical shell problem
          !
          if (lroot) print*,'init_ss: kws temperature in spherical shell'
          call shell_ss(f)

        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'init_ss: No such value for initss(' &
                           //trim(iinit_str)//'): ',trim(initss(iinit))
          call stop_it("")

      endselect

      if (lroot) print*,'init_ss: initss(' &
                        //trim(iinit_str)//') = ',trim(initss(iinit))

      endif

      enddo

      if (lnothing.and.lroot) print*,'init_ss: zero entropy'
!
!  if ss_const/=0, add this constant to entropy
!  (ss_const is already taken care of)
!
!     if (ss_const/=0) f(:,:,:,iss)=f(:,:,:,iss)+ss_const
!
!  no entropy initialization when lgravr=.true.
!  why?
!
!  The following seems insane, so I comment this out.
!     if (lgravr) then
!       f(:,:,:,iss) = -0.
!     endif

!
!  Add perturbation(s)

!
!      if (lgravz)

      select case (pertss)

      case('zero', '0')
        ! Don't perturb

      case ('hexagonal', '1')
        !
        !  hexagonal perturbation
        !
        if (lroot) print*,'init_ss: adding hexagonal perturbation to ss'
        f(:,:,:,iss) = f(:,:,:,iss) &
                        + ampl_ss*(2*cos(sqrt(3.)*0.5*khor_ss*xx) &
                                    *cos(0.5*khor_ss*yy) &
                                   + cos(khor_ss*yy) &
                                  ) * cos(pi*zz)

      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*,'init_ss: No such value for pertss:', pertss
        call stop_it("")

      endselect
!
      if(ip==0) print*,xx,yy  !(to keep compiler quiet)
!
    endsubroutine init_ss
!***********************************************************************
    subroutine polytropic_ss_z( &
         f,mpoly,zz,tmp,zint,zbot,zblend,isoth,cs2int,ssint)
!
!  Implement a polytropic profile in ss above zbot. If this routine is
!  called several times (for a piecewise polytropic atmosphere), on needs
!  to call it from top to bottom.
!
!  zint    -- z at top of layer
!  zbot    -- z at bottom of layer
!  zblend  -- smoothly blend (with width widthss) previous ss (for z>zblend)
!             with new profile (for z<zblend)
!  isoth   -- flag for isothermal stratification;
!  ssint   -- value of ss at interface, i.e. at the top on entry, at the
!             bottom on exit
!  cs2int  -- same for cs2
!
      use Sub, only: step
      use Gravity, only: gravz
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: tmp,p,zz
      real, dimension (mz) :: stp
      real :: mpoly,zint,zbot,zblend,beta1,cs2int,ssint
      integer :: isoth
!
      ! NB: beta1 is not dT/dz, but dcs2/dz = (gamma-1)c_p dT/dz
      if (isoth /= 0) then ! isothermal layer
        beta1 = 0.
        tmp = ssint - gamma1*gravz*(zz-zint)/cs2int
        ssint = -gamma1*gravz*(zbot-zint)/cs2int ! ss at layer interface
      else
        beta1 = gamma*gravz/(mpoly+1)
        tmp = 1 + beta1*(zz-zint)/cs2int
        tmp = max(tmp,epsi)  ! ensure arg to log is positive
        tmp = ssint + (1-mpoly*gamma1)/gamma &
                      * alog(tmp)
        ssint = ssint + (1-mpoly*gamma1)/gamma & ! ss at layer interface
                        * alog(1 + beta1*(zbot-zint)/cs2int)
      endif
      cs2int = cs2int + beta1*(zbot-zint) ! cs2 at layer interface (bottom)

      !
      ! smoothly blend the old value (above zblend) and the new one (below
      ! zblend) for the two regions:
      !
      stp = step(z,zblend,widthss)
      p = spread(spread(stp,1,mx),2,my)
      f(:,:,:,iss) = p*f(:,:,:,iss)  + (1-p)*tmp
!
    endsubroutine polytropic_ss_z
!***********************************************************************
    subroutine polytropic_ss_disc( &
         f,mpoly,zz,tmp,zint,zbot,zblend,isoth,cs2int,ssint)
!
!  Implement a polytropic profile in ss for a disc. If this routine is
!  called several times (for a piecewise polytropic atmosphere), on needs
!  to call it from bottom (middle of disc) to top.
!
!  zint    -- z at bottom of layer
!  zbot    -- z at top of layer (naming convention follows polytropic_ss_z)
!  zblend  -- smoothly blend (with width widthss) previous ss (for z>zblend)
!             with new profile (for z<zblend)
!  isoth   -- flag for isothermal stratification;
!  ssint   -- value of ss at interface, i.e. at the bottom on entry, at the
!             top on exit
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
      real :: mpoly,zint,zbot,zblend,beta1,cs2int,ssint, nu_epicycle2
      integer :: isoth
!
      ! NB: beta1 is not dT/dz, but dcs2/dz = (gamma-1)c_p dT/dz
      nu_epicycle2 = nu_epicycle**2
      if (isoth /= 0) then ! isothermal layer
        beta1 = 0.
        tmp = ssint - gamma1*gravz*nu_epicycle2*(zz**2-zint**2)/cs2int/2.
        ssint = -gamma1*gravz*nu_epicycle2*(zbot**2-zint**2)/cs2int/2. 
              ! ss at layer interface
      else
        beta1 = gamma*gravz*nu_epicycle2/(mpoly+1)
        tmp = 1 + beta1*(zz**2-zint**2)/cs2int/2.
        tmp = max(tmp,epsi)  ! ensure arg to log is positive
        tmp = ssint + (1-mpoly*gamma1)/gamma &
                      * alog(tmp)
        ssint = ssint + (1-mpoly*gamma1)/gamma & ! ss at layer interface
                        * alog(1 + beta1*(zbot**2-zint**2)/cs2int/2.)
      endif
      cs2int = cs2int + beta1*(zbot**2-zint**2)/2. 
             ! cs2 at layer interface (bottom)

      !
      ! smoothly blend the old value (above zblend) and the new one (below
      ! zblend) for the two regions:
      !
      stp = step(z,zblend,widthss)
      p = spread(spread(stp,1,mx),2,my)
      f(:,:,:,iss) = p*f(:,:,:,iss)  + (1-p)*tmp
!
    endsubroutine polytropic_ss_disc
!***********************************************************************
    subroutine hydrostatic_isentropic(f,lnrho_bot,ss_const)
!
!  Hydrostatic initial condition at constant entropy.
!  Full ionization equation of state.
!
!  Solves dlnrho/dz=gravz/cs2 using 2nd order Runge-Kutta.
!  Currently only works for vertical gravity field.
!  Starts at bottom boundary where the density has to be set in the gravity
!  module.
!
!  This should be done in the density module but entropy has to initialize
!  first.
!
!
!  20-feb-04/tobi: coded
!
      use Ionization, only: pressure_gradient
      use Gravity, only: gravz
      use Mpicomm, only: ipz,stop_it

      real, dimension (mx,my,mz,mvar+maux), intent(inout) :: f
      real, intent(in) :: lnrho_bot,ss_const
      real :: cs2,cp1tilde,lnrho,lnrho_m

      if (.not. lgravz) then
        call stop_it("hydrostatic_isentropic: Currently only works for vertical gravity field")
      endif

      !
      ! In case this processor is not located at the very bottom
      ! perform integration through lower lying processors
      !
      lnrho=lnrho_bot
      do n=1,nz*ipz
        call pressure_gradient(lnrho,ss_const,cs2,cp1tilde)
        lnrho_m=lnrho+dz*gravz/cs2/2
        call pressure_gradient(lnrho_m,ss_const,cs2,cp1tilde)
        lnrho=lnrho+dz*gravz/cs2
      enddo

      !
      ! Do the integration on this processor
      !
      f(:,:,n1,ilnrho)=lnrho
      do n=n1+1,n2
        call pressure_gradient(lnrho,ss_const,cs2,cp1tilde)
        lnrho_m=lnrho+dz*gravz/cs2/2
        call pressure_gradient(lnrho_m,ss_const,cs2,cp1tilde)
        lnrho=lnrho+dz*gravz/cs2
        f(:,:,n,ilnrho)=lnrho
      enddo

      !
      ! Entropy is simply constant
      !
      f(:,:,:,iss)=ss_const

    endsubroutine hydrostatic_isentropic
!***********************************************************************
    subroutine shell_ss(f)
!
!  Initialize entropy based on specified radial temperature profile in
!  a spherical shell
!
!  20-oct-03/dave -- coded
!
      use Gravity, only: g0
      use Density, only: mpoly
      use Ionization, only: getentropy
      use Sub, only: calc_unitvects_sphere

      real, dimension (mx,my,mz,mvar+maux), intent(inout) :: f
      real, dimension (nx) :: lnrho,lnTT,TT,ss
      real :: beta1
!
      beta1 = g0/(mpoly+1)
      if (initss(1)=='geo-kws') then
        do m=m1,m2
        do n=n1,n2
!
          call calc_unitvects_sphere()
!
          where (r_mn >= r_ext) TT = TT_ext
          where (r_mn < r_ext .AND. r_mn > r_int) TT = 1+beta1*(1/r_mn-1)
!         where (r_mn < r_ext .AND. r_mn > r_int) TT = gamma/gamma1*(1+beta1*(1/r_mn-1))
!         goes with alternate scaling in initialize_entropy
          where (r_mn <= r_int) TT = TT_int
!
          lnrho=f(l1:l2,m,n,ilnrho)
          lnTT=log(TT)
          call getentropy(lnrho,lnTT,ss)
          f(l1:l2,m,n,iss)=ss
!
        enddo 
        enddo 
      endif
!      
   end subroutine shell_ss
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
      use Ionization
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension(nx) :: absz,n_c,n_w,n_i,n_h
!  T in K, k_B s.t. pp is in code units ( = 9.59e-15 erg/cm/s^2)
!  (i.e. k_B = 1.381e-16 (erg/K) / 9.59e-15 (erg/cm/s^2) )
      real, parameter :: T_c_cgs=500.0,T_w_cgs=8.0e3,T_i_cgs=8.0e3,T_h_cgs=1.0e6 
      real :: T_c,T_w,T_i,T_h
      real, dimension(nx) :: rho,pp,lnrho,ss,lnTT,yH
      real :: cp1tilde,mu 
!      real, dimension(nx) :: pp 
      double precision :: pp0 
!      real, dimension(2) :: fmpi2
      real, dimension(1) :: fmpi1
      real :: kpc, rhoscale
!      integer :: iproctop
!
      if (lroot) print*,'ferriere: Ferriere density and entropy profile'
!
!  first define reference values of pp, cs2, at midplane.  
!  pressure is set to 6 times thermal pressure, this factor roughly
!  allowing for other sources, as modelled by Ferriere.
!
      call getmu(mu)
      kpc = 3.086D21 / unit_length
      rhoscale = m_H * mu * unit_length**3
!      print *,'rhoscale: ',rhoscale, mu
      T_c=T_c_cgs/unit_temperature
      T_w=T_w_cgs/unit_temperature
      T_i=T_i_cgs/unit_temperature
      T_h=T_h_cgs/unit_temperature

!      pp0=6.0*k_B*(rho0/1.38) *                                               &
!       (1.09*0.340*T_c + 1.09*0.226*T_w + 2.09*0.025*T_i + 2.27*0.00048*T_h)
!      pp0=k_B*unit_length**3*                                               &
!       (1.09*0.340*T_c + 1.09*0.226*T_w + 2.09*0.025*T_i + 2.27*0.00048*T_h)
!      cs20=gamma*pp0/rho0
!      cs0=sqrt(cs20)
!      ss0=alog(gamma*pp0/cs20/rho0)/gamma   !ss0=zero  (not needed)
!
      do n=n1,n2            ! nb: don't need to set ghost-zones here
      absz=abs(z(n))
      do m=m1,m2 
!  cold gas profile n_c (eq 6)
        n_c=0.340*(0.859*exp(-(z(n)*kpc/0.127)**2) +         &
                   0.047*exp(-(z(n)*kpc/0.318)**2) +         &
                   0.094*exp(-absz*kpc/0.403))     
!  warm gas profile n_w (eq 7)
        n_w=0.226*(0.456*exp(-(z(n)*kpc/0.127)**2) +  &
                   0.403*exp(-(z(n)*kpc/0.318)**2) +  &
                   0.141*exp(-absz*kpc/0.403))
!  ionized gas profile n_i (eq 9)
        n_i=0.0237*exp(-absz*kpc) + 0.0013* exp(-absz*kpc/0.150)
!  hot gas profile n_h (eq 13)
        n_h=0.00048*exp(-absz*kpc/1.5)         
!  normalised s.t. rho0 gives mid-plane density directly (in 10^-24 g/cm^3)
        !rho=rho0/(0.340+0.226+0.025+0.00048)*(n_c+n_w+n_i+n_h)*rhoscale
        rho=(n_c+n_w+n_i+n_h)*rhoscale
        lnrho=alog(rho)
        f(l1:l2,m,n,ilnrho)=lnrho

!  define entropy via pressure, assuming fixed T for each component
        if(lentropy) then
!  thermal pressure (eq 15)
          pp=k_B*unit_length**3 *                                 &
             (1.09*n_c*T_c + 1.09*n_w*T_w + 2.09*n_i*T_i + 2.27*n_h*T_h)
!           
          call eoscalc(ilnrho_pp,lnrho,pp,ss=ss,yH=yH,lnTT=lnTT) 
          call pressure_gradient(lnrho(1),ss(1),cs2bot,cp1tilde)
          call pressure_gradient(lnrho(nx),ss(nx),cs2top,cp1tilde)
!
          f(l1:l2,m,n,iss)=ss
!        
          fmpi1=(/ cs2bot /)
          call mpibcast_real(fmpi1,1,0)
          cs2bot=fmpi1(1) 
          fmpi1=(/ cs2top /)
          call mpibcast_real(fmpi1,1,ncpus-1)
          cs2top=fmpi1(1)
!
        endif
       enddo
      enddo
!      
      if (lroot) print*, 'ferriere: cs2bot=',cs2bot, ' cs2top=',cs2top
!
    endsubroutine ferriere
!***********************************************************************
    subroutine shock2d(f,xx,yy,zz)
!
!  shock2d
!
! taken from clawpack:
!     =====================================================
!       subroutine ic2rp2(maxmx,maxmy,meqn,mbc,mx,my,x,y,dx,dy,q)
!     =====================================================
!
!     # Set initial conditions for q.
!
!      # Data is piecewise constant with 4 values in 4 quadrants
!      # 2D Riemann problem from Figure 4 of
!        @article{csr-col-glaz,
!          author="C. W. Schulz-Rinne and J. P. Collins and H. M. Glaz",
!          title="Numerical Solution of the {R}iemann Problem for
!                 Two-Dimensional Gas Dynamics",
!          journal="SIAM J. Sci. Comput.",
!          volume="14",
!          year="1993",
!          pages="1394-1414" }
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real, dimension (4) :: rpp,rpr,rpu,rpv
!
      if (lroot) print*,'shock2d: initial condition, gamma=',gamma
!
!      # First quadrant:
        rpp(1) = 1.5d0
        rpr(1) = 1.5d0
        rpu(1) = 0.d0
        rpv(1) = 0.d0
!
!      # Second quadrant:
        rpp(2) = 0.3d0
        rpr(2) = 0.532258064516129d0
        rpu(2) = 1.206045378311055d0
        rpv(2) = 0.0d0
!
!      # Third quadrant:
        rpp(3) = 0.029032258064516d0
        rpr(3) = 0.137992831541219d0
        rpu(3) = 1.206045378311055d0
        rpv(3) = 1.206045378311055d0
!
!      # Fourth quadrant:
        rpp(4) = 0.3d0
        rpr(4) = 0.532258064516129d0
        rpu(4) = 0.0d0
        rpv(4) = 1.206045378311055d0
!
!  s=-lnrho+alog(gamma*p)/gamma
!
        where ( (xx>=0.) .and. (yy>=0.) )
          f(:,:,:,ilnrho)=alog(rpr(1))
          f(:,:,:,iss)=alog(gamma*rpp(1))/gamma-f(:,:,:,ilnrho)
          f(:,:,:,iux)=rpu(1)
          f(:,:,:,iuy)=rpv(1)
        endwhere
        where ( (xx<0.) .and. (yy>=0.) )
          f(:,:,:,ilnrho)=alog(rpr(2))
          f(:,:,:,iss)=alog(gamma*rpp(2))/gamma-f(:,:,:,ilnrho)
          f(:,:,:,iux)=rpu(2)
          f(:,:,:,iuy)=rpv(2)
        endwhere
        where ( (xx<0.) .and. (yy<0.) )
          f(:,:,:,ilnrho)=alog(rpr(3))
          f(:,:,:,iss)=alog(gamma*rpp(3))/gamma-f(:,:,:,ilnrho)
          f(:,:,:,iux)=rpu(3)
          f(:,:,:,iuy)=rpv(3)
        endwhere
        where ( (xx>=0.) .and. (yy<0.) )
          f(:,:,:,ilnrho)=alog(rpr(4))
          f(:,:,:,iss)=alog(gamma*rpp(4))/gamma-f(:,:,:,ilnrho)
          f(:,:,:,iux)=rpu(4)
          f(:,:,:,iuy)=rpv(4)
        endwhere
!
    endsubroutine shock2d
!**********************************************************************
    subroutine dss_dt(f,df,uu,glnrho,divu,rho1,lnrho,cs2,TT1,shock,gshock,bb,bij)
!
!  Entropy auxiliary variables and diagnostics  
!
!   9-jun-04/anders: adapted from entropy
!
      use Cdata
      use Ionization
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: bij
      real, dimension (nx,3) :: uu,glnrho,gshock,bb
      real, dimension (nx) :: divu,rho1,lnrho,cs2,shock
      real, dimension (nx) :: ss,lnTT,TT1,cp1tilde,ee,rho
!
      intent(in) :: f,uu,glnrho,rho1,lnrho,shock,gshock
      intent(out) :: cs2,TT1
!
!  identify module 
!
      if (headtt.or.ldebug) print*,'dss_dt: fixed entropy'
!
!  calculate the necessary thermodynamics
!  yH and TT have already been calculated in the beginning of pencil loop
!
      ss=f(l1:l2,m,n,iss)
      call eoscalc(f,nx,ee=ee,lnTT=lnTT)
      call pressure_gradient(f,cs2,cp1tilde)
      TT1=exp(-lnTT)
!
!  phi-averages
!  Note that this does not necessarily happen with ldiagnos=.true.
!
      if (l2davgfirst) then
        call phisum_mn_name_rz(ss,i_ssmphi)
      endif
!
!  Calculate entropy related diagnostics
!
      if(ldiagnos) then
        if (i_dtc/=0) call max_mn_name(sqrt(cs2)/dxmin/cdt,i_dtc,l_dt=.true.)
        rho=exp(lnrho)
        if(i_eth/=0) call sum_mn_name(rho*ee,i_eth)
        if(i_ethtot/=0) call integrate_mn_name(rho*ee,i_ethtot)
        if(i_ethdivum/=0) call sum_mn_name(rho*ee*divu,i_ethdivum)
        if(i_ssm/=0) call sum_mn_name(ss,i_ssm)
      endif
!
    endsubroutine dss_dt
!***********************************************************************
    subroutine rprint_entropy(lreset,lwrite)
!
!  reads and registers print parameters relevant to entropy
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Cdata
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
        i_dtc=0; i_eth=0; i_ethdivum=0; i_ssm=0; i_ugradpm=0; i_ethtot=0
        i_dtchi=0
        i_ssmphi=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtc',i_dtc)
        call parse_name(iname,cname(iname),cform(iname),'dtchi',i_dtchi)
        call parse_name(iname,cname(iname),cform(iname),'ethtot',i_ethtot)
        call parse_name(iname,cname(iname),cform(iname),'ethdivum',i_ethdivum)
        call parse_name(iname,cname(iname),cform(iname),'eth',i_eth)
        call parse_name(iname,cname(iname),cform(iname),'ssm',i_ssm)
        call parse_name(iname,cname(iname),cform(iname),'ugradpm',i_ugradpm)
      enddo
!
!  check for those quantities for which we want phi-averages
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'ssmphi',i_ssmphi)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_dtc=',i_dtc
        write(3,*) 'i_dtchi=',i_dtchi
        write(3,*) 'i_ethtot=',i_ethtot
        write(3,*) 'i_ethdivum=',i_ethdivum
        write(3,*) 'i_eth=',i_eth
        write(3,*) 'i_ssm=',i_ssm
        write(3,*) 'i_ugradpm=',i_ugradpm
        write(3,*) 'nname=',nname
        write(3,*) 'iss=',iss
        write(3,*) 'i_ssmphi=',i_ssmphi
      endif
!
    endsubroutine rprint_entropy
!***********************************************************************
endmodule Entropy

