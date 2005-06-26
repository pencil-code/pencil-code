! $Id: entropy.f90,v 1.333 2005-06-26 07:00:15 bingert Exp $

!  This module takes care of entropy (initial condition
!  and time advance)

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
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
  real :: TT_int,TT_ext,cs2_int,cs2_ext,cool_int=0.,cool_ext=0.,ampl_TT=0.
  real :: chi=0.,chi_t=0.,chi_shock=0.,Kappa0=0.,Kisotr=0.
  real :: Kgperp=0.,Kgpara=0.,tdown=1,allp=2
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
  real :: heat_uniform=0.,cool_RTV=0.
  real :: deltaT_poleq=0.
  logical :: lturbulent_heat=.false.
  logical :: lcalc_heatcond_simple=.false.,lmultilayer=.true.
  logical :: lcalc_heatcond=.false.,lcalc_heatcond_constchi=.false.
  logical :: lupw_ss=.false.
  logical :: lgaspressuregradient=.true.
  character (len=labellen), dimension(ninit) :: initss='nothing'
  character (len=labellen) :: pertss='zero'
  character (len=labellen) :: cooltype='Temp',iheatcond='K-const'
  
  !
  ! Parameters for subroutine cool_RTV in SI units (from Cook et Al. 1989)
  ! 
  double precision, parameter, dimension (10) :: & 
       intlnT_1 =(/4.605, 8.959, 9.906, 10.534, 11.283, 12.434, 13.286, 14.541, 17.51, 20.723 /) 
  double precision, parameter, dimension (9) :: &
       lnH_1 = (/ -190.884,  -141.165, -80.245, -101.314, -78.748, -53.88, -80.452, -70.758, -91.182/), &
       B_1   = (/     11.7,      6.15,      0.,      2.0,      0.,    -2.,      0., -0.6667,    0.5 /)  
  !
  ! A second set of parameters for cool_RTV (from interstellar.f90)
  !
  double precision, parameter, dimension(7) ::  &
       intlnT_2 = (/ 5.704,7.601 , 8.987 , 11.513 , 17.504 , 20.723, 24.0 /) 
  double precision, parameter, dimension(6) ::  &
       lnH_2 = (/-102.811, -99.01, -111.296, -70.804, -90.934, -80.572 /),   &
       B_2   = (/    2.0,     1.5,   2.867,  -0.65,   0.5, 0.0 /)

  ! input parameters
  namelist /entropy_init_pars/ &
       initss,pertss,grads0,radius_ss,ampl_ss,widthss,epsilon_ss, &
       ss_left,ss_right,ss_const,mpoly0,mpoly1,mpoly2,isothtop, &
       khor_ss,thermal_background,thermal_peak,thermal_scaling,cs2cool, &
       center1_x, center1_y, center1_z, center2_x, center2_y, center2_z, &
       T0,ampl_TT,kx_ss

  ! run parameters
  namelist /entropy_run_pars/ &
       hcond0,hcond1,hcond2,widthss, &
       luminosity,wheat,cooltype,cool,cs2cool,rcool,wcool,Fbot, &
       chi_t,chi_shock,chi,Kappa0,Kisotr,iheatcond, &
       Kgperp,Kgpara, cool_RTV, &
       lcalc_heatcond_simple,lcalc_heatcond,lcalc_heatcond_constchi,&
       tau_ss_exterior,lmultilayer,Kbot,tau_cor,TT_cor,z_cor, &
       tauheat_buffer,TTheat_buffer,zheat_buffer,dheat_buffer1, &
       heat_uniform,lupw_ss,lcalc_cp,cool_int,cool_ext, &
       lturbulent_heat, deltaT_poleq, lgaspressuregradient, &
       tdown, allp
  
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
      lentropy = .true.
!
      iss = nvar+1             ! index to access entropy
      nvar = nvar+1
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_entropy: nvar = ', nvar
        print*, 'register_entropy: iss = ', iss
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: entropy.f90,v 1.333 2005-06-26 07:00:15 bingert Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_entropy: nvar > mvar')
      endif
!
!  Put variable name in array
!
      varname(iss) = 'ss'
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',ss $'
          if (nvar == mvar) write(4,*) ',ss'
        else
          write(4,*) ',ss $'
        endif
        write(15,*) 'ss = fltarr(mx,my,mz)*one'
      endif
!
    endsubroutine register_entropy
!***********************************************************************
    subroutine initialize_entropy(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  21-jul-2002/wolf: coded
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
        case('geo-kws','geo-benchmark')
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
      if (ip==0) print*,f,lstarting  !(to keep compiler quiet)
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
      real :: cs2int,ssint,ztop,ss_ext,pot0,pot_ext
      logical :: lnothing=.true.
!
      intent(in) :: xx,yy,zz
      intent(inout) :: f
!
      if (pretend_lnTT) f(:,:,:,iss)=f(:,:,:,iss)+(f(:,:,:,ilnrho)*gamma1-log(gamma1))/gamma
      do iinit=1,ninit
!
!  if we pretend that ss in in reality g1lnTT, we initialize the background
!  of lnTT/gamma such that it corresponds to ss=0.
!

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
        case('wave'); f(:,:,:,iss)=ss_const+ampl_ss*sin(kx_ss*xx(:,:,:) + pi)
        case('Ferriere'); call ferriere(f) 
        case('xjump'); call jump(f,iss,ss_left,ss_right,widthss,'x')
        case('yjump'); call jump(f,iss,ss_left,ss_right,widthss,'y')
        case('zjump'); call jump(f,iss,ss_left,ss_right,widthss,'z')
        case('hor-fluxtube'); call htube(ampl_ss,f,iss,iss,xx,yy,zz,radius_ss,epsilon_ss)
        case('hor-tube'); call htube2(ampl_ss,f,iss,iss,xx,yy,zz,radius_ss,epsilon_ss)

        case('sedov') 
          if (lroot) print*,'init_ss: sedov - thermal background with gaussian energy burst'
        call blob(thermal_peak,f,iss,radius_ss,center1_x,center1_y,center1_z)
      !   f(:,:,:,iss) = f(:,:,:,iss) + (log(f(:,:,:,iss) + thermal_background)+log(thermal_scaling))/gamma 

        case('sedov-dual') 
          if (lroot) print*,'init_ss: sedov - thermal background with gaussian energy burst'
        call blob(thermal_peak,f,iss,radius_ss,center1_x,center1_y,center1_z)
        call blob(thermal_peak,f,iss,radius_ss,center2_x,center2_y,center2_z)
      !   f(:,:,:,iss) = (log(f(:,:,:,iss) + thermal_background)+log(thermal_scaling))/gamma 
   
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
          ! ss0=log(-gamma1*gravz*zinfty)/gamma
          ! print*,'init_ss: isentropic stratification; ss=',ss0
          f(:,:,:,iss)=0.
          if (ampl_ss/=0.) then
            print*,'init_ss: put bubble: radius_ss,ampl_ss=',radius_ss,ampl_ss
            tmp=xx**2+yy**2+zz**2
            f(:,:,:,iss)=f(:,:,:,iss)+ampl_ss*exp(-tmp/max(radius_ss**2-tmp,1e-20))
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
              ss_ext = 0. + log(cs2cool/cs2_ext)
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

        case ('geo-benchmark')
          !
          ! radial temperature profiles for spherical shell problem
          !
          if (lroot) print*,'init_ss: benchmark temperature in spherical shell'
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
                      * log(tmp)
        ssint = ssint + (1-mpoly*gamma1)/gamma & ! ss at layer interface
                        * log(1 + beta1*(zbot-zint)/cs2int)
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
                      * log(tmp)
        ssint = ssint + (1-mpoly*gamma1)/gamma & ! ss at layer interface
                        * log(1 + beta1*(zbot**2-zint**2)/cs2int/2.)
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
      real, dimension (nx) :: lnrho,lnTT,TT,ss,pert_TT
      real :: beta1
!
      beta1 = g0/(mpoly+1)
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
!
        call calc_unitvects_sphere()
        call shell_ss_perturb(pert_TT)
!
        where (r_mn >= r_ext) TT = TT_ext
        where (r_mn < r_ext .AND. r_mn > r_int) TT = 1+beta1*(1/r_mn-1)+pert_TT
!       where (r_mn < r_ext .AND. r_mn > r_int) TT = gamma/gamma1*(1+beta1*(1/r_mn-1))
!       goes with alternate scaling in initialize_entropy
        where (r_mn <= r_int) TT = TT_int
!
        lnrho=f(l1:l2,m,n,ilnrho)
        lnTT=log(TT)
        call getentropy(lnrho,lnTT,ss)
        f(l1:l2,m,n,iss)=ss
!
      enddo 
!      
    endsubroutine shell_ss
!***********************************************************************
    subroutine shell_ss_perturb(pert_TT)
!
!  Compute perturbation to initial temperature profile
!
!  22-june-04/dave -- coded
!
      use Sub, only: calc_phiavg_general

      real, dimension (nx), intent(out) :: pert_TT
      real, dimension (nx) :: xr,cos_4phi,sin_theta4
      real :: ampl0=.885065
!
      select case(initss(1))
!
        case ('geo-kws')
          pert_TT=0.
!
        case ('geo-benchmark')
          call calc_phiavg_general()
          xr=2*r_mn-r_int-r_ext              ! radial part of perturbation
          cos_4phi=cos(4*phi_mn)             ! azimuthal part
          sin_theta4=(rcyl_mn/r_mn)**4       ! meridional part
          pert_TT=ampl0*ampl_TT*(1-3*xr**2+3*xr**4-xr**6)*sin_theta4*cos_4phi
!
      endselect
!
    endsubroutine shell_ss_perturb
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
!     double precision :: pp0 
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
!      ss0=log(gamma*pp0/cs20/rho0)/gamma   !ss0=zero  (not needed)
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
        lnrho=log(rho)
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
!  s=-lnrho+log(gamma*p)/gamma
!
        where ( (xx>=0.) .and. (yy>=0.) )
          f(:,:,:,ilnrho)=log(rpr(1))
          f(:,:,:,iss)=log(gamma*rpp(1))/gamma-f(:,:,:,ilnrho)
          f(:,:,:,iux)=rpu(1)
          f(:,:,:,iuy)=rpv(1)
        endwhere
        where ( (xx<0.) .and. (yy>=0.) )
          f(:,:,:,ilnrho)=log(rpr(2))
          f(:,:,:,iss)=log(gamma*rpp(2))/gamma-f(:,:,:,ilnrho)
          f(:,:,:,iux)=rpu(2)
          f(:,:,:,iuy)=rpv(2)
        endwhere
        where ( (xx<0.) .and. (yy<0.) )
          f(:,:,:,ilnrho)=log(rpr(3))
          f(:,:,:,iss)=log(gamma*rpp(3))/gamma-f(:,:,:,ilnrho)
          f(:,:,:,iux)=rpu(3)
          f(:,:,:,iuy)=rpv(3)
        endwhere
        where ( (xx>=0.) .and. (yy<0.) )
          f(:,:,:,ilnrho)=log(rpr(4))
          f(:,:,:,iss)=log(gamma*rpp(4))/gamma-f(:,:,:,ilnrho)
          f(:,:,:,iux)=rpu(4)
          f(:,:,:,iuy)=rpv(4)
        endwhere
!
    if (ip==0) print*,zz
    endsubroutine shock2d
!**********************************************************************
    subroutine dss_dt(f,df,uu,divu,lnrho,rho,rho1,glnrho,cs2,TT1,shock,gshock,bb, bij)
!
!  calculate right hand side of entropy equation
!  heat condution is currently disabled until old stuff,
!  which in now in calc_heatcond, has been reinstalled.
!
!  17-sep-01/axel: coded
!   9-jun-02/axel: pressure gradient added to du/dt already here
!   2-feb-03/axel: added possibility of ionization
!
      use Cdata
      use Mpicomm
      use Ionization
      use Sub
      use Global
      use Special, only: special_calc_entropy
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: bij,hss,hlnrho,hlnTT
      real, dimension (nx,3) :: uu,glnrho,gss,gshock,glnTT,bb
      real, dimension (nx) :: ugss,uglnrho,divu,rhs
      real, dimension (nx) :: lnrho,ss,rho,rho1,cs2,yH,lnTT,TT1,cp1tilde
      real, dimension (nx) :: ee,shock
      real, dimension (nx) :: vKpara,vKperp
      real :: zbot,ztop,xi,profile_cor
      integer :: j,ju
!
      intent(in)  :: f,uu,glnrho,divu,rho,rho1,lnrho,shock,gshock,bb,bij
      intent(out) :: df,cs2,TT1
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dss_dt: SOLVE dss_dt'
      if (headtt) call identify_bcs('ss',iss)
!
!  define bottom and top z positions
!
      zbot=xyz0(3)
      ztop=xyz0(3)+Lxyz(3)
!
!  entropy gradient: needed for advection and pressure gradient
!
      call grad(f,iss,gss)
!
!  calculate the necessary thermodynamics
!  yH and TT have already been calculated in the beginning of pencil loop
!
      ss=f(l1:l2,m,n,iss)
      call eoscalc(f,nx,ee=ee,lnTT=lnTT)
      call pressure_gradient(f,lnrho,cs2,cp1tilde)
      call temperature_gradient(f,glnrho,gss,glnTT)
      TT1=exp(-lnTT)
!
!  calculate cs2, TT1, and cp1tilde in a separate routine
!  With IONIZATION=noionization, assume perfect gas with const coeffs
!
      if (headtt) print*,'dss_dt: lnTT,cs2,cp1tilde=',lnTT(1),cs2(1),cp1tilde(1)
!
!  ``cs2/dx^2'' for timestep
!
      if (lfirst.and.ldt) advec_cs2=cs2*dxyz_2
      if (headtt.or.ldebug) print*,'dss_dt: max(advec_cs2) =',maxval(advec_cs2)
!
!  don't print out cs20 (in runs with ionization is in not defined)
!  and ideally it should be avoided from the entropy module altoghether
!
!     if (headtt) print*, &
!                'dss_dt: entropy (but not used with ionization): cs20=',cs20
!
      if (lhydro) then
!
!  subtract pressure gradient term in momentum equation
!  (allow suppression of pressure gradient for test purposes)
!
        if (lgaspressuregradient) then
          do j=1,3
            ju=j+iuu-1
            if (pretend_lnTT) then
              df(l1:l2,m,n,ju)=df(l1:l2,m,n,ju)-cs2*(glnrho(:,j)/gamma+cp1tilde*gss(:,j))
            else
              df(l1:l2,m,n,ju)=df(l1:l2,m,n,ju)-cs2*(glnrho(:,j)+cp1tilde*gss(:,j))
            endif
           enddo
        endif
!
!  add pressure gradient force from isothermal background density change in x
!
        if (dlnrhobdx .ne. 0) then
          df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-cs2*dlnrhobdx/gamma
        endif
!
!  velocity damping in the coronal heating zone
!
        if (tau_cor>0) then
          if (z(n)>=z_cor) then
            xi=(z(n)-z_cor)/(ztop-z_cor)
            profile_cor=xi**2*(3-2*xi)
            df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz) &
                             -profile_cor*f(l1:l2,m,n,iux:iuz)/tau_cor
          endif
        endif
!
      endif
!
!  advection term
!
      call u_dot_gradf(f,iss,gss,uu,ugss,UPWIND=lupw_ss)
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) - ugss
!
!  if pretend_lnTT=.true., we pretend that ss is actually lnTT/gamma
!
      if (pretend_lnTT) then
        df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)-divu*gamma1/gamma
      endif
!
!ajwm - lviscosity always true and there is not a noviscosity module
      if (lviscosity) call calc_viscous_heat(f,df,glnrho,divu,rho1,cs2,TT1,shock)
!
!  thermal conduction
!
!  For backward compatibility: map old logicals onto iheatcond
!
      if (lcalc_heatcond_simple) then
        if (headtt) print*,"OBSOLETE: use iheatcond='simple' instead of lcalc_heatcond_simple=T"
        iheatcond='simple'
      elseif (lcalc_heatcond_constchi) then
        if (headtt) print*,"OBSOLETE: use iheatcond='chi-const' instead of lcalc_heatcond_constchi=T"
        iheatcond='chi-const'
      elseif (lcalc_heatcond) then
        if (headtt) print*,"OBSOLETE: use iheatcond='K-const' instead of lcalc_heatcond=T"
        iheatcond='K-const'
      endif
!
!  the new scheme using iheatcond
!
      select case(iheatcond)
      case('K-const')
         call calc_heatcond(f,df,rho1,glnrho,gss)
      case('simple')
         call calc_heatcond_simple(f,df,rho1,glnrho,gss)
      case('chi-const')
         call calc_heatcond_constchi(f,df,rho1,glnrho,gss)
      case ('tensor-diffusion')
        call g2ij(f,iss,hss)        
        if (pretend_lnTT) then
          glnTT=gss
          hlnTT=hss
        else
          call g2ij(f,ilnrho,hlnrho)
          call temperature_hessian(f,hlnrho,hss,hlnTT)
        endif
        vKpara(:) = Kgpara
        vKperp(:) = Kgperp
        call tensor_diffusion_coef(glnTT,hlnTT,bij,bb,vKperp,vKpara,rhs,llog=.true.)
        df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)+rhs*rho1
      case ('spitzer')
        call calc_heatcond_spitzer(f,df,rho1,glnrho,gss,bb,bij,cs2)
      case ('corona')
        call calc_heatcond_spitzer(f,df,rho1,glnrho,gss,bb,bij,cs2)
        call newton_cool(f,df,lnTT,lnrho)
        call calc_heat_cool_RTV(f,df,lnrho,lnTT)     
      case ('shock')
         call calc_heatcond_shock(f,df,rho1,glnrho,glnTT,gss,shock,gshock)
      case default
         if (lroot) then
            print*,'dss_dt: No such value iheatcond = ', trim(iheatcond)
         endif
         call stop_it("")
      endselect
!
!  heating/cooling
!
      if ((luminosity /= 0) .or. &
          (cool /= 0) .or. &
          (tau_cor /= 0) .or. &
          (tauheat_buffer /= 0) .or. &
          (heat_uniform /= 0) .or. &
          (cool_ext /= 0 .AND. cool_int /= 0) .or. &
          (lturbulent_heat)) &
        call calc_heat_cool(f,df,rho1,cs2,cp1tilde,ss,TT1)
!
!  interstellar radiative cooling and UV heating
!
      if (linterstellar) &
        call calc_heat_cool_interstellar(df,rho,rho1,TT1,yH)
!
!  possibility of entropy relaxation in exterior region
!
      if (tau_ss_exterior/=0.) call calc_tau_ss_exterior(f,df)
!
!  entry possibility for "personal" entries.
!  In that case you'd need to provide your own "special" routine.
!
      if (lspecial) call special_calc_entropy(f,df,uu,glnrho,divu,rho1,lnrho,cs2,TT1)
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
        if (i_dtc/=0) call max_mn_name(sqrt(advec_cs2)/cdt,i_dtc,l_dt=.true.)
        if(i_eth/=0) call sum_mn_name(rho*ee,i_eth)
        if(i_ethtot/=0) call integrate_mn_name(rho*ee,i_ethtot)
        if(i_ethdivum/=0) call sum_mn_name(rho*ee*divu,i_ethdivum)
        if(i_ssm/=0) call sum_mn_name(ss,i_ssm)
        if(i_ugradpm/=0) then
          call dot_mn(uu,glnrho,uglnrho)
          call sum_mn_name(cs2*(uglnrho+ugss),i_ugradpm)
        endif
      endif
!
    endsubroutine dss_dt
!***********************************************************************
    subroutine calc_heatcond_constchi(f,df,rho1,glnrho,gss)
!
!  Heat conduction for constant value of chi=K/(rho*cp)
!  This routine also adds in turbulent diffusion, if chi_t /= 0.
!  Ds/Dt = ... + 1/(rho*T) grad(flux), where
!  flux = chi*rho*gradT + chi_t*rho*T*grads
!  This routine is currently not correct when ionization is used.
!
!  29-sep-02/axel: adapted from calc_heatcond_simple
!
      use Cdata
      use Mpicomm
      use Sub
      use Gravity
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: glnrho,gss,glnT,glnP
      real, dimension (nx) :: rho1
      real, dimension (nx) :: thdiff,del2ss,del2lnrho,g2
!
      intent(in) :: f,glnrho,gss
      intent(out) :: df
!
!  check that chi is ok
!
      if(headtt) print*,'calc_heatcond_constchi: chi==',chi
!
!  Heat conduction
!  Note: these routines require revision when ionization turned on
!
      call del2(f,iss,del2ss)
      call del2(f,ilnrho,del2lnrho)
      glnT = gamma*gss + gamma1*glnrho
      glnP = gamma*gss + gamma*glnrho
      call dot_mn(glnP,glnT,g2)
      thdiff = chi * (gamma*del2ss+gamma1*del2lnrho + g2)
      if(chi_t/=0.) then
        call dot_mn(glnP,gss,g2)
        thdiff = thdiff + chi_t*(del2ss+g2)
      endif
!
!  add heat conduction to entropy equation
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
      if (headtt) print*,'calc_heatcond_constchi: added thdiff'
!
!  check maximum diffusion from thermal diffusion
!  With heat conduction, the second-order term for entropy is
!  gamma*chi*del2ss
!
      if (lfirst.and.ldt) then
        diffus_chi=max(diffus_chi,(gamma*chi+chi_t)*dxyz_2)
        if (ldiagnos.and.i_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,i_dtchi,l_dt=.true.)
        endif
      endif
!
      if(ip==0) print*,rho1 !(to keep compiler quiet)
    endsubroutine calc_heatcond_constchi
!***********************************************************************
    subroutine calc_heatcond_shock(f,df,rho1,glnrho,glnTT,gss,shock,gshock)
!
!  Adds in shock entropy diffusion. There is potential for
!  recycling some quantities from previous calculations.
!  Ds/Dt = ... + 1/(rho*T) grad(flux), where
!  flux = chi_shock*rho*T*grads
!  (in comments we say chi_shock, but in the code this is "chi_shock*shock")
!  This routine should be ok with ionization.
!
!  20-jul-03/axel: adapted from calc_heatcond_constchi
!  19-nov-03/axel: added chi_t also here.
!
      use Cdata
      use Mpicomm
      use Sub
      use Gravity
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: glnrho,gss,gshock,glnTT
      real, dimension (nx) :: rho1,shock
      real, dimension (nx) :: thdiff,del2ss,g2,gshockgss
!
      intent(in) :: f,glnrho,gss,shock,gshock
      intent(out) :: df
!
!  check that chi is ok
!
      if(headtt) print*,'calc_heatcond_shock: chi_t,chi_shock=',chi_t,chi_shock
!
!  calculate terms for shock diffusion
!  Ds/Dt = ... + chi_shock*[del2ss + (glnchi_shock+glnpp).gss]
!
      call del2(f,iss,del2ss)
      call dot_mn(gshock,gss,gshockgss)
      call dot_mn(glnTT+glnrho,gss,g2)
!
!  shock entropy diffusivity
!  Write: chi_shock = chi_shock0*shock, and gshock=grad(shock), so
!  Ds/Dt = ... + chi_shock0*[shock*(del2ss+glnpp.gss) + gshock.gss]
!
      if(headtt) print*,'calc_heatcond_shock: use shock diffusion'
      thdiff=(chi_shock*shock+chi_t)*(del2ss+g2)+chi_shock*gshockgss
!
!  add heat conduction to entropy equation
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
      if (headtt) print*,'calc_heatcond_shock: added thdiff'
!
!  check maximum diffusion from thermal diffusion
!  With heat conduction, the second-order term for entropy is
!  gamma*chi*del2ss
!
      if (lfirst.and.ldt) then
        diffus_chi=max(diffus_chi,(chi_t+chi_shock*shock)*dxyz_2)
        if (ldiagnos.and.i_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,i_dtchi,l_dt=.true.)
        endif
      endif
!
      if(ip==0) print*,rho1,glnrho !(keep compiler quiet)
    endsubroutine calc_heatcond_shock
!***********************************************************************
    subroutine calc_heatcond_simple(f,df,rho1,glnrho,gss)
!
!  heat conduction
!
!   8-jul-02/axel: adapted from Wolfgang's more complex version
!
      use Cdata
      use Mpicomm
      use Sub
      use Gravity
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: glnrho,gss,glnT,glnThcond !,glhc
      real, dimension (nx) :: rho1,chix
      real, dimension (nx) :: thdiff,del2ss,del2lnrho,g2
      real, dimension (nx) :: hcond
!
      intent(in) :: f,rho1,glnrho,gss
      intent(out) :: df
!
!  This particular version assumes a simple polytrope, so mpoly is known
!
      hcond=Kbot
      if(headtt) print*,'calc_heatcond_simple: max(hcond)=',maxval(hcond)
!
!  Heat conduction
!  Note: these routines require revision when ionization turned on
!
      call del2(f,iss,del2ss)
      call del2(f,ilnrho,del2lnrho)
      chix = rho1*hcond
      glnT = gamma*gss + gamma1*glnrho ! grad ln(T)
      glnThcond = glnT !... + glhc/spread(hcond,2,3)    ! grad ln(T*hcond)
      call dot_mn(glnT,glnThcond,g2)
      thdiff = chix * (gamma*del2ss+gamma1*del2lnrho + g2)
!
!  add heat conduction to entropy equation
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
      if (headtt) print*,'calc_heatcond_simple: added thdiff'
!
!  check maximum diffusion from thermal diffusion
!  With heat conduction, the second-order term for entropy is
!  gamma*chix*del2ss
!
      if (lfirst.and.ldt) then
        diffus_chi=max(diffus_chi,(gamma*chix+chi_t)*dxyz_2)
        if (ldiagnos.and.i_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,i_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_simple
!***********************************************************************
    subroutine calc_heatcond_spitzer(f,df,rho1,glnrho,gss,bb,bij,cs2)
!
!  Calculates heat conduction parallel and perpendicular (isotropic)
!  to magnetic field lines     
!
!  See: Solar MHD; Priest 1982
!
!  10-feb-04/bing: coded
!
      use Ionization
      use Sub
      use Io
      use Mpicomm
!       
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: hlnrho,hss,hlnTT,bij
      real, dimension (nx,3) :: glnrho,gss,glnT,bb,gvKpara,gvKperp
      real, dimension (nx,3) :: tmpv1,tmpv2,tmpv3
      real, dimension (nx) :: rho1,TT,bb2,cs2,thdiff,b1
      real, dimension (nx) :: tmps,quenchfactor,vKpara,vKperp
      integer ::i,j
!
      intent(in) :: gss,glnrho,f
      intent(out) :: df
!
      glnT = gamma*gss + gamma1*glnrho
!
      TT = cs2/gamma1   
!        
!      Kgpara/Kgperp are in SI units -> get them in pencil units:
!      unit_temp = (0.667 * gamma1 * unit_velocity**2 )/8.3144e3 /gamma
!      Kgpara = Kgpara / (unit_velocity**3 * unit_density * unit_length / unit_temp)

      call g2ij(f,ilnrho,hlnrho)
      call g2ij(f,iss,hss)
      call temperature_hessian(f,hlnrho,hss,hlnTT)
!
!     Calculate variable diffusion coefficients along pencils
!
      call dot2_mn(bb,bb2)
      b1=1./max(tiny(bb2),bb2)
!
      vKpara = Kgpara * TT**2.5 * rho1      != Kgpara* T^3.5              /(rho*T)
      vKperp = Kgperp * (b1/rho1)/sqrt(TT)  != Kgperp* rho^2 sqrt(T)/B^2  /(rho*T)
!
      quenchfactor = vKpara/(vKpara+vKperp)
      vKperp=vKperp*quenchfactor
!     do i=1,nx
!        if  (vKperp(i) .gt. vKpara(i)) vKperp(i) = vKpara(i)
!     enddo
      if (notanumber(vKpara)) print*,'vKpara'
      if (notanumber(vKperp)) print*,'vKperp'
!
!     Calculate gradient of diff. coeffs
!      
      tmps = 3.5 * vKpara 
      call multsv_mn(tmps,glnT,gvKpara)
      tmps(:) = 2.
      call multsv_mn(tmps,glnrho,tmpv1)
      tmps(:) = 0.5
      call multsv_mn(tmps, glnT,tmpv2)
      do i=1,3
         do j=1,3
            tmpv3(:,i)=2*bb(:,j)*bij(:,j,i)
         end do
      end do
!
      call multsv_mn(b1,tmpv3,tmpv3)
      tmpv1=tmpv1+tmpv2+tmpv3
      call multsv_mn(vKperp,tmpv1,gvKperp)
      gvKperp=gvKperp*spread(quenchfactor,2,3)
!
!     Calculate diffusion term
!
      call  tensor_diffusion_coef(glnT,hlnTT,bij,bb,vKperp,vKpara,thdiff,GVKPERP=gvKperp,GVKPARA=gvKpara)
!
!      if (notanumber(hlnTT)) call stop_it('calc_heatcond_spitzer: hlnTT')
!     
!      (thdiff = thdiff*rho1/TT) is included in vKperp and vKpara
! 
      if (notanumber(thdiff)) then
         print*,ipx,ipy,ipz
         call stop_it('calc_heatcond_spitzer: thdiff')
      endif
!
!      thdiff = thdiff/TT*rho1
!
      df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) + thdiff 
!
!   check maximum diffusion from thermal diffusion
!   With heat conduction, the second-order term for entropy is
!   gamma*chi*del2ss
!
      if (lfirst.and.ldt) then
         diffus_chi=max(diffus_chi,(gamma*Kgpara*rho1*TT**2.5+chi_t)*dxyz_2)
         if (ldiagnos.and.i_dtchi/=0) then
            call max_mn_name(diffus_chi/cdtv,i_dtchi,l_dt=.true.)
         endif
      endif
!      
    endsubroutine calc_heatcond_spitzer
!***********************************************************************
    subroutine calc_heatcond(f,df,rho1,glnrho,gss)
!
!  heat conduction
!
!  17-sep-01/axel: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use IO
      use Gravity
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: glnrho,gss,glnT,glnThcond,glhc
      real, dimension (nx) :: rho1,chix
      real, dimension (nx) :: thdiff,del2ss,del2lnrho,g2
      real, dimension (nx) :: hcond
      real :: z_prev=-1.23e20
!
      save :: z_prev,hcond,glhc
!
      intent(in) :: f,rho1,glnrho,gss
      intent(out) :: df
!
!  Heat conduction / entropy diffusion
!

      if ((hcond0 /= 0) .or. (chi_t /= 0)) then
        call del2(f,iss,del2ss)
      endif

      if (hcond0 /= 0) then
        if(headtt) then
          print*,'calc_heatcond: hcond0=',hcond0
          if (lgravz) print*,'calc_heatcond: Fbot,Ftop=',Fbot,Ftop
        endif
        if (lgravz) then
          if(headtt) print*,'calc_heatcond: lgravz=',lgravz
          ! For vertical geometry, we only need to calculate this for each
          ! new value of z -> speedup by about 8% at 32x32x64
          if (z_mn(1) /= z_prev) then
            call heatcond(hcond)
            call gradloghcond(glhc)
            z_prev = z_mn(1)
          endif
        else
          call heatcond(hcond)       ! returns hcond=hcond0
          call gradloghcond(glhc)    ! returns glhc=0
        endif
        call del2(f,ilnrho,del2lnrho)
        chix = rho1*hcond
        glnT = gamma*gss + gamma1*glnrho             ! grad ln(T)
        glnThcond = glnT + glhc/spread(hcond,2,3)    ! grad ln(T*hcond)
        call dot_mn(glnT,glnThcond,g2)
        thdiff = chix * (gamma*del2ss+gamma1*del2lnrho + g2)
      else
        chix = 0
        thdiff = 0
        ! not really needed, I (wd) guess -- but be sure before you
        ! remove them
        hcond = 0
        glhc = 0
      endif
!
!  "turbulent" entropy diffusion
!
      if (chi_t/=0.) then
        if (headtt) then
          print*,'calc_headcond: "turbulent" entropy diffusion: chi_t=',chi_t
          if (hcond0 /= 0) then
            print*, "calc_heatcond: WARNING ", &
                  "- hcond0 and chi_t combined don't seem to make sense"
          endif
        endif
!        df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss)+chi_t*del2ss
        thdiff = chi_t*del2ss
      endif
!
!  check for NaNs initially
!
      if (headt .and. (hcond0 /= 0)) then
        if (notanumber(glhc))      print*,'calc_heatcond: NaNs in glhc'
        if (notanumber(rho1))      print*,'calc_heatcond: NaNs in rho1'
        if (notanumber(hcond))     print*,'calc_heatcond: NaNs in hcond'
        if (notanumber(chix))      print*,'calc_heatcond: NaNs in chix'
        if (notanumber(del2ss))    print*,'calc_heatcond: NaNs in del2ss'
        if (notanumber(del2lnrho)) print*,'calc_heatcond: NaNs in del2lnrho'
        if (notanumber(glhc))      print*,'calc_heatcond: NaNs in glhc'
        if (notanumber(1/hcond))   print*,'calc_heatcond: NaNs in 1/hcond'
        if (notanumber(glnT))      print*,'calc_heatcond: NaNs in glnT'
        if (notanumber(glnThcond)) print*,'calc_heatcond: NaNs in glnThcond'
        if (notanumber(g2))        print*,'calc_heatcond: NaNs in g2'
        if (notanumber(thdiff))    print*,'calc_heatcond: NaNs in thdiff'
        !
        !  most of these should trigger the following trap
        !
        if (notanumber(thdiff)) then
          print*, 'calc_heatcond: m,n,y(m),z(n)=',m,n,y(m),z(n)
          call stop_it('calc_heatcond: NaNs in thdiff')
        endif
      endif

      if (headt .and. lfirst .and. ip<=9) then
        call output_pencil(trim(directory)//'/chi.dat',chix,1)
        call output_pencil(trim(directory)//'/hcond.dat',hcond,1)
        call output_pencil(trim(directory)//'/glhc.dat',glhc,3)
      endif
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
!
      if (headtt) print*,'calc_heatcond: added thdiff'
!
!  check maximum diffusion from thermal diffusion
!  NB: With heat conduction, the second-order term for entropy is
!    gamma*chix*del2ss
!
      if (lfirst.and.ldt) then
        diffus_chi=max(diffus_chi,(gamma*chix+chi_t)*dxyz_2)
        if (ldiagnos.and.i_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,i_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond
!***********************************************************************
    subroutine calc_heat_cool(f,df,rho1,cs2,cp1tilde,ss,TT1)
!
!  add combined heating and cooling
!
!  02-jul-02/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use Gravity
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: rho1,cs2,ss,cp1tilde,TT1
      real, dimension (nx) :: heat,prof,theta_profile
      real :: ssref,zbot,ztop,profile_buffer,xi,profile_cor
!
      intent(in) :: f,rho1,cs2
      intent(out) :: df
!
!  identifier
!
      if(headtt) print*,'calc_heat_cool: lgravz=',lgravz
!
!  define bottom and top height
!
      zbot=xyz0(3)
      ztop=xyz0(3)+Lxyz(3)
!
!  initialize
!
      heat=0.
!
!  Vertical case:
!  Heat at bottom, cool top layers
!
      if (lgravz .and. (luminosity .ne. 0. .or. cool .ne. 0.)) then
!
!  TEMPORARY: Add heat near bottom. Wrong: should be Heat/(T*rho)
!
        ! heating profile, normalised, so volume integral = 1
        prof = spread(exp(-0.5*((z(n)-zbot)/wheat)**2), 1, l2-l1+1) &
             /(sqrt(pi/2.)*wheat*Lx*Ly)
        heat = luminosity*prof
        ! smoothly switch on heating if required
        if ((ttransient > 0) .and. (t < ttransient)) then
          heat = heat * t*(2*ttransient-t)/ttransient**2
        endif
        ! cooling profile; maximum = 1
        ssref = ss0 + (-log(gamma) + log(cs20))/gamma + grads0*ztop
        prof = spread(exp(-0.5*((ztop-z(n))/wcool)**2), 1, l2-l1+1)
        heat = heat - cool*prof*(cs2-cs20)/cs20
      endif
!
!  Spherical case:
!  heat at centre, cool outer layers
!
      if (lgravr) then
        ! central heating
        ! heating profile, normalised, so volume integral = 1
        prof = exp(-0.5*(r_mn/wheat)**2) * (2*pi*wheat**2)**(-1.5)
        heat = luminosity*prof
        ! surface cooling; entropy or temperature
        ! cooling profile; maximum = 1
!        prof = 0.5*(1+tanh((r_mn-1.)/wcool))
        prof = step(r_mn,rcool,wcool)
        !
        !  pick type of cooling
        !
        select case(cooltype)
        case ('cs2', 'Temp')    ! cooling to reference temperature cs2cool
          heat = heat - cool*prof*(cs2-cs2cool)/cs2cool
        case ('cs2-rho', 'Temp-rho') ! cool to reference temperature cs2cool
                                     ! in a more time-step neutral manner
          heat = heat - cool*prof*(cs2-cs2cool)/cs2cool/rho1
        case ('entropy')        ! cooling to reference entropy (currently =0)
          heat = heat - cool*prof*(f(l1:l2,m,n,iss)-0.)
! dgm
        case ('shell')          !  heating/cooling at shell boundaries
          heat=0.                            ! default
          select case(initss(1))
            case ('geo-kws'); heat=0.        ! can add heating later based on value of initss
          endselect
          !
          !  possibility of a latitudinal heating profile
          !  T=T0-(2/3)*delT*P2(costheta), for testing Taylor-Proudman theorem
          !  Note that P2(x)=(1/2)*(3*x^2-1).
          !
          if (deltaT_poleq/=0.) then
            if(headtt) print*,'calc_heat_cool: deltaT_poleq=',deltaT_poleq
            if(headtt) print*,'rcyl_mn=',rcyl_mn
            if(headtt) print*,'z_mn=',z_mn
            theta_profile=(1./3.-(rcyl_mn/z_mn)**2)*deltaT_poleq
            prof = step(r_mn,r_ext,wcool)      ! outer heating/cooling step
            heat = heat - cool_ext*prof*(cs2-cs2_ext)/cs2_ext*theta_profile
            prof = 1 - step(r_mn,r_int,wcool)  ! inner heating/cooling step
            heat = heat - cool_int*prof*(cs2-cs2_int)/cs2_int*theta_profile
          else
            prof = step(r_mn,r_ext,wcool)      ! outer heating/cooling step
            heat = heat - cool_ext*prof*(cs2-cs2_ext)/cs2_ext
            prof = 1 - step(r_mn,r_int,wcool)  ! inner heating/cooling step
            heat = heat - cool_int*prof*(cs2-cs2_int)/cs2_int
          endif
!
        case default
          if (lroot) print*, &
               'calc_heat_cool: No such value for cooltype: ', trim(cooltype)
          call stop_it("")
        endselect
      endif
!
!  add spatially uniform heating (usually as a test)
!
      if(heat_uniform/=0.) heat=heat+heat_uniform
!
!  add "coronal" heating (to simulate a hot corona)
!  assume a linearly increasing reference profile
!  This 1/rho1 business is clumpsy, but so would be obvious alternatives...
!
      if(tau_cor>0) then
        if(z(n)>=z_cor) then
          xi=(z(n)-z_cor)/(ztop-z_cor)
          profile_cor=xi**2*(3-2*xi)
          heat=heat+profile_cor*(TT_cor-1/TT1)/(rho1*tau_cor*cp1tilde)
        endif
      endif
!
!  add heating and cooling to a reference temperature in a buffer
!  zone at the z boundaries. Only regions in |z| > zheat_buffer are affected.
!  Inverse width of the transition is given by dheat_buffer1.
!
      if(tauheat_buffer/=0.) then
        profile_buffer=0.5*(1.+tanh(dheat_buffer1*(z(n)**2-zheat_buffer**2)))
!        profile_buffer=1.+0.5*(tanh(dheat_buffer1*(z(n)-z(n1)-zheat_buffer)) + tanh(dheat_buffer1*(z(n)-z(n2)-zheat_buffer)))
        heat=heat+profile_buffer*ss*(TTheat_buffer-1/TT1)/(rho1*tauheat_buffer)
      endif
!
!  Parametrized turbulent heating
!
      if(lturbulent_heat) then
        df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + TT1*nu_turb*(qshear*Omega)**2
      endif
!
!  add to entropy equation
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + TT1*rho1*heat
!
    endsubroutine calc_heat_cool
!***********************************************************************
    subroutine calc_heat_cool_RTV(f,df,lnrho,lnTT)
!
!    calculate cool term:  C = ne*ni*Q(T) 
!    with ne*ni = 1.2*np^2 = 1.2*rho^2/(1.4*mp)^2
!    Q(T) = H*T^B is piecewice poly
!
!  15-dez-04/bing: coded
!
      use IO
!     
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: lnTT,lncool,lnrho,lnneni,rtv_cool,lnTT_SI
      integer :: i,imax
      real :: unit_temp,unit_dens,unit_power
      real :: p0,p1,p2,p3,p4
!
      intent(in) :: f,lnrho,lnTT
      intent(out) :: df
!
!     All is in SI units and has to be rescaled to PENCIL units
!
      unit_temp = (0.667 * gamma1 * unit_velocity**2 )/8.3144e3 /gamma
      unit_dens = unit_density
      unit_power = unit_density * unit_velocity**3 / unit_length
!
      lnTT_SI = lnTT + alog(unit_temp) 
!
! First set of parameters
      if (cool_RTV .gt. 0.) then
         imax = size(intlnT_1,1)       
         lncool(:)=0.0
         do i=1,imax-1
            where (( intlnT_1(i) <= lnTT_SI .or. i==1 ) .and. lnTT_SI < intlnT_1(i+1) )
               lncool=lncool + lnH_1(i) + B_1(i)*lnTT_SI
            endwhere
         enddo
         where (lnTT_SI >= intlnT_1(imax) )
            lncool = lncool + lnH_1(imax-1) + B_1(imax-1)*intlnT_1(imax)
         endwhere
      endif

! Second set of parameters      
      if (cool_RTV .lt. 0) then
         cool_RTV = cool_RTV*(-1.)
        imax = size(intlnT_2,1)       
        lncool(:)=0.0
        do i=1,imax-1
          where (( intlnT_2(i) <= lnTT_SI .or. i==1 ) .and. lnTT_SI < intlnT_2(i+1) )
            lncool=lncool + lnH_2(i) + B_2(i)*lnTT_SI
          endwhere
        enddo
        where (lnTT_SI >= intlnT_2(imax) )
          lncool = lncool + lnH_2(imax-1) + B_2(imax-1)*intlnT_2(imax)
        endwhere
     endif
!
!     calculate ln(ne*ni) :
!          mp = 1.673*1e-27
!          ln(ne*ni) = ln( 1.2*rho^2/(1.4*mp)^2)
!          lnneni   = 2*lnrho + alog(1.2) - 2*alog(1.4*1.673*1e-27)
!          lnneni = 2*lnrho + 122.82     
!      
!    rtv_cool=exp(lnneni+lncool-lnrho-lnTT)/unit_power
!    =>
     rtv_cool=exp(lncool+lnrho-lnTT+122.82)/unit_power
!     
     rtv_cool=rtv_cool * cool_RTV  ! just for adjusting by setting cool_RTV in run.in
!
!     add to entropy equation
!
     df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss)-rtv_cool 
!
    endsubroutine calc_heat_cool_RTV
!***********************************************************************
    subroutine calc_tau_ss_exterior(f,df)
!
!  entropy relaxation to zero on time scale tau_ss_exterior within
!  exterior region. For the time being this means z > zgrav.
!
!  29-jul-02/axel: coded
!
      use Cdata
      use Gravity
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: scl
!
      intent(in) :: f
      intent(out) :: df
!
      if (headtt) print*,'calc_tau_ss_exterior: tau=',tau_ss_exterior
      if(z(n)>zgrav) then
        scl=1./tau_ss_exterior
        df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)-scl*f(l1:l2,m,n,iss)
      endif
!
    endsubroutine calc_tau_ss_exterior
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
    subroutine heatcond(hcond)
!
!  calculate the heat conductivity hcond along a pencil.
!  This is an attempt to remove explicit reference to hcond[0-2] from
!  code, e.g. the boundary condition routine.
!
!  NB: if you modify this profile, you *must* adapt gradloghcond below.
!
!  23-jan-2002/wolf: coded
!  18-sep-2002/axel: added lmultilayer switch
!
      use Sub, only: step
      use Gravity
!
      real, dimension (nx) :: hcond
!
      if (lgravz) then
        if (lmultilayer) then
          hcond = 1 + (hcond1-1)*step(z_mn,z1,-widthss) &
                    + (hcond2-1)*step(z_mn,z2,widthss)
          hcond = hcond0*hcond
        else
          hcond=Kbot
        endif
      else
        hcond = hcond0
      endif
!
    endsubroutine heatcond
!***********************************************************************
    subroutine gradloghcond(glhc)
!
!  calculate grad(log hcond), where hcond is the heat conductivity
!  NB: *Must* be in sync with heatcond() above.
!  23-jan-2002/wolf: coded
!
      use Sub, only: der_step
      use Gravity
!
      real, dimension (nx,3) :: glhc
!
      if (lgravz) then
        glhc(:,1:2) = 0.
        glhc(:,3) = (hcond1-1)*der_step(z_mn,z1,-widthss) &
                    + (hcond2-1)*der_step(z_mn,z2,widthss)
        glhc(:,3) = hcond0*glhc(:,3)
      else
        glhc = 0.
      endif
!
    endsubroutine gradloghcond
!***********************************************************************
    subroutine newton_cool(f,df,lnTT,lnrho)
!
!  Keeps the temperature in the lower chromosphere and the upper corona
!  at a constant level using newton cooling
!
!  15-dez-2004/bing: coded
!

      use Sub
      use IO
      use Mpicomm

      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: lnTT,lnrho,newton
      real :: lnTTor,xil,unit_temp
      real :: p0,p1,p2,p3,p4
!    
!     Initial temperature profile is given in ln(T) over z in Mm
!     It is independent of grid and unit system
!     Since I do not change initial condition this works fine
!
      p0 = 2.47955
      p1 = 4.45524
      p2 = 1.51496
      p3 = 4.44825
      p4 = 2.89396e-03
!
!     Get the heigth in Mm
!
      xil =  z(n) * unit_length * 1e-6 
!
!     Calculate ln(T) in SI
!
      lnTTor = p0*(tanh((xil-p1)/p2) +p3 +p4*xil)
!
      unit_temp = (0.667 * gamma1 * unit_velocity**2 )/8.3144e3 /gamma
!
      lnTTor = lnTTor - alog(unit_temp) 
!
      newton = (exp(lnrho-lnrho0)) ** allp  * tdown / gamma 
!
      newton = newton * (lnTT-lnTTor)
!      
!     Add cooling term to entropy
!
      if (notanumber(newton))  call stop_it('newton cool: NaN`s in newton') 
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) - newton
!
    endsubroutine newton_cool
!***********************************************************************    
endmodule Entropy
 
