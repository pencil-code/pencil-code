! $Id: entropy.f90,v 1.206 2003-10-07 14:20:10 mee Exp $

!  This module takes care of entropy (initial condition
!  and time advance)

module Entropy

  use Cparam
  use Cdata
  use Hydro
  use Interstellar
  use Viscosity
  use Density, only: lcalc_cp
  use Ionization, only: lionization_fixed,xHe

  implicit none

  !real, dimension (nx) :: cs2,TT1
  real :: radius_ss=0.1,ampl_ss=0.,widthss=2*epsi,epsilon_ss
  real :: luminosity=0.,wheat=0.1,cs2cool=0.,cool=0.,rcool=1.,wcool=0.1
  real :: chi=0.,chi_t=0.,chi_shock=0.
  real :: ss_left,ss_right
  real :: ss0=0.,khor_ss=1.,ss_const=0.
  real :: tau_ss_exterior=0.
  !parameters for Sedov type initial condition
  real :: center1_x=0., center1_y=0., center1_z=0.
  real :: center2_x=0., center2_y=0., center2_z=0.
  real :: thermal_background=0., thermal_peak=0., thermal_scaling=1.
  real :: hcond0=0.
  real :: hcond1=impossible,hcond2=impossible
  real :: Fbot=impossible,FbotKbot=impossible,Kbot=impossible
  real :: Ftop=impossible,FtopKtop=impossible,Ktop=impossible
  real :: tauheat_coronal=0.,TTheat_coronal=0.,zheat_coronal=0.
  real :: tauheat_buffer=0.,TTheat_buffer=0.,zheat_buffer=0.,dheat_buffer1=0.
  real :: heat_uniform=0.
  logical :: lcalc_heatcond_simple=.false.,lmultilayer=.true.
  logical :: lcalc_heatcond_constchi=.false.
  logical :: lupw_ss=.false.
  character (len=labellen) :: initss='nothing',pertss='zero',cooltype='Temp'

  ! input parameters
  namelist /entropy_init_pars/ &
       initss,pertss,grads0,radius_ss,ampl_ss,widthss,epsilon_ss, &
       ss_left,ss_right,ss_const,mpoly0,mpoly1,mpoly2,isothtop, &
       khor_ss, thermal_background, thermal_peak, thermal_scaling, &
       center1_x, center1_y, center1_z, &
       center2_x, center2_y, center2_z
     

  ! run parameters
  namelist /entropy_run_pars/ &
       hcond0,hcond1,hcond2,widthss, &
       luminosity,wheat,cooltype,cool,cs2cool,rcool,wcool,Fbot, &
       chi_t,chi_shock,lcalc_heatcond_simple,tau_ss_exterior, &
       chi,lcalc_heatcond_constchi,lmultilayer,Kbot, &
       tauheat_coronal,TTheat_coronal,zheat_coronal, &
       tauheat_buffer,TTheat_buffer,zheat_buffer,dheat_buffer1, &
       heat_uniform, lupw_ss, lcalc_cp

  ! other variables (needs to be consistent with reset list below)
  integer :: i_eth=0,i_TTm=0,i_yHm=0,i_ssm=0,i_ugradpm=0, i_ethtot=0

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
           "$Id: entropy.f90,v 1.206 2003-10-07 14:20:10 mee Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_entropy: nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      if (maux == 0) then
         if (nvar < mvar) write(4,*) ',ss $'
         if (nvar == mvar) write(4,*) ',ss'
      else
        write(4,*) ',ss $'
     endif
      write(5,*) 'ss = fltarr(mx,my,mz)*one'
!
    endsubroutine register_entropy
!***********************************************************************
    subroutine initialize_entropy()
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  21-jul-2002/wolf: coded
!
      use Cdata
      use Gravity, only: gravz
!
      lneed_sij = .true.   !let Hydro module know to precalculate some things
      lneed_glnrho = .true.
!
!  radiative diffusion: initialize flux etc
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
            else
              Fbot=0.
            endif
          endif
          Kbot=gamma1/gamma*(mpoly+1.)*Fbot
          FbotKbot=gamma/gamma1/(mpoly+1.)
          if(lroot) print*,'initialize_entropy: Fbot,Kbot=',Fbot,Kbot
          !
          !  calculate Ftop if it has not been set in run.in
          !
          if (Ftop==impossible) then
            if (bcz2(iss)=='c1') then
              Ftop=-gamma/(gamma-1)*hcond0*gravz/(mpoly+1)
              if (lroot) print*, &
                      'initialize_entropy: Calculated Ftop = ', Ftop
            else
              Ftop=0.
            endif
          endif
          Ktop=gamma1/gamma*(mpoly+1.)*Ftop
          FtopKtop=gamma/gamma1/(mpoly+1.)
          if(lroot) print*,'initialize_entropy: Ftop,Ktop=',Ftop,Ktop
!
        endif
      endif
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
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: tmp,xx,yy,zz
      real :: cs2int,ssint,ztop
!
      intent(in) :: xx,yy,zz
      intent(inout) :: f
!
      select case(initss)
        case('nothing'); if(lroot) print*,'init_ss: nothing'
        case('zero', '0'); f(:,:,:,iss) = 0.
        case('const_ss'); f(:,:,:,iss) = ss_const
        case('blob'); call blob(ampl_ss,f,iss,radius_ss,0.,0.,0.)
        case('isothermal'); call isothermal_entropy(f)
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

        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'init_ss: No such value for initss: ', trim(initss)
          call stop_it("")

      endselect
!
!  if ss_const/=0, add this constant to entropy
!  (ss_const is already taken care of)
!
!     if (ss_const/=0) f(:,:,:,iss)=f(:,:,:,iss)+ss_const
!
!  no entropy initialization when lgravr=.true.
!  why?
!
      if (lgravr) then
        f(:,:,:,iss) = -0.
      endif

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
    subroutine isothermal_entropy(f)
!
!  Isothermal stratification (for lnrho and ss)
!  This routine should be independent of the gravity module used.
!  When entropy is present, this module also initializes entropy.
!
!  Sound speed (and hence Temperature), is
!  initialised to the reference value:
!           sound speed: cs^2_0            from start.in  
!           density: rho0 = exp(lnrho0)
!
!  11-jun-03/tony: extracted from isothermal routine in Density module
!                  to allow isothermal condition for arbitrary density
!
!
      use Mpicomm, only: stop_it
      use Gravity
      use Ionization

      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if (lionization.or.lionization_fixed) &
       call stop_it("isothermal_entropy: NOT IMPLEMENTED FOR IONIZATION CASES")
      do n=n1,n2
      do m=m1,m2
          f(l1:l2,m,n,iss)= -gamma1*(f(l1:l2,m,n,ilnrho)-lnrho0)/gamma
                  ! + other terms for sound speed not equal to cs_0
      enddo
      enddo

!
!  cs2 values at top and bottom may be needed to boundary conditions.
!  The values calculated here may be revised in the entropy module.
!
      cs2bot=cs20
      cs2top=cs20
     
!
    endsubroutine isothermal_entropy
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
      real :: absz,n_c,n_w,n_i,n_h
!  T in K, k_B s.t. pp is in code units ( = 9.59e-15 erg/cm/s^2)
!  (i.e. k_B = 1.381e-16 (erg/K) / 9.59e-15 (erg/cm/s^2) )
      real :: T_c=500.0,T_w=8.0e3,T_i=8.0e3,T_h=1.0e6 !,k_B=0.0144
      real :: rho,lnrho,pp,pp0,ss,TT,yH 
      real, dimension(2) :: fmpi2
!      real, dimension(1) :: fmpi1
!      integer :: iproctop
!
      if (lroot) print*,'ferriere: Ferriere density and entropy profile'
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
           
          if (lionization) then
            call perturb_mass(lnrho,pp,ss,TT,yH) 
            f(l1:l2,m,n,iss)=ss
            !  calculate cs2bot,top: needed for a2/c2 b.c.s (fixed T)
            if (n == n1 .and. m == m1) cs2bot=0.
            if (n == n2 .and. m == m1) cs2top=0.
          else          
            f(l1:l2,m,n,iss)=alog(gamma*pp/cs20)/gamma +                   &
                                     gamma1/gamma*lnrho0 - lnrho
            !  calculate cs2bot,top: needed for a2/c2 b.c.s (fixed T)
            if (n == n1 .and. m == m1) cs2bot=gamma*pp/rho
            if (n == n2 .and. m == m1) cs2top=gamma*pp/rho
          endif
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
!**********************************************************************
    subroutine dss_dt(f,df,uu,glnrho,divu,rho1,lnrho,cs2,TT1)
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
      real, dimension (nx,3) :: uu,glnrho,gss
      real, dimension (nx) :: ugss,uglnrho,divu
      real, dimension (nx) :: lnrho,ss,yH,TT,rho1,cs2,TT1,cp1tilde
      real, dimension (nx) :: rho,ee
      integer :: j,ju
!
      intent(in) :: f,uu,glnrho,rho1,lnrho
      intent(out) :: df,cs2,TT1
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dss_dt: SOLVE dss_dt'
      if (headtt) call identify_bcs('ss',iss)
!
!  entropy gradient: needed for advection and pressure gradient
!
      call grad(f,iss,gss)
!
!  calculate the necessary thermodynamics
!  yH and TT have already been calculated in the beginning of pencil loop
!
      ss=f(l1:l2,m,n,iss)
      call ionget(f,yH,TT)
      call thermodynamics(lnrho,ss,yH,TT,cs2=cs2,cp1tilde=cp1tilde,ee=ee)
!
!  calculate cs2, TT1, and cp1tilde in a separate routine
!  With IONIZATION=noionization, assume perfect gas with const coeffs
!
      if (headtt) print*,'dss_dt: cs2,TT,cp1tilde=',cs2(1),TT(1),cp1tilde(1)
!
!  use sound speed in Courant condition
!
      if (lfirst.and.ldt) maxadvec2=amax1(maxadvec2,cs2)
      if (ip<8.and.lroot.and.imn==1) print*, &
                        'dss_dt: maxadvec2,cs2=',maxadvec2,cs2
      if (headtt) print*, &
                 'dss_dt: entropy (but not used with ionization): cs20=',cs20
!
!  subtract pressure gradient term in momentum equation
!
      if (lhydro) then
        do j=1,3
          ju=j+iuu-1
          df(l1:l2,m,n,ju)=df(l1:l2,m,n,ju)-cs2*(glnrho(:,j)+cp1tilde*gss(:,j))
        enddo
      endif
!
!  advection term
!
      call u_dot_gradf(f,iss,gss,uu,ugss,UPWIND=lupw_ss)
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) - ugss
!
!  T is now calculated in thermodynamics, calculate 1/T (==TT1)
!  Viscous heating depends on ivisc; no visc heating if ivisc='simplified'
!  
      TT1=1./TT
!
!ajwm - lviscosity always true and there is not a noviscosity module
      if (lviscosity) call calc_viscous_heat(f,df,glnrho,divu,rho1,cs2,TT1)
!
!  thermal conduction
!
      if (lcalc_heatcond_simple) then
        call calc_heatcond_simple(f,df,rho1,glnrho,gss)
      elseif (lcalc_heatcond_constchi) then
        call calc_heatcond_constchi(f,df,rho1,glnrho,gss)
      else
        !if ((hcond0 /= 0) .or. (chi_t /= 0)) &  !!(AB: needed?)
        call calc_heatcond(f,df,rho1,glnrho,gss)
      endif
!
!  shock entropy diffusion
!
      if (chi_shock/=0.) call calc_heatcond_shock(f,df,rho1,glnrho,gss)
!
!  heating/cooling
!
      if ((luminosity /= 0) .or. &
          (cool /= 0) .or. &
          (tauheat_coronal /= 0) .or. &
          (tauheat_buffer /= 0) .or. &
          (heat_uniform /= 0)) &
        call calc_heat_cool(f,df,rho1,cs2,ss,TT,TT1)
!
!  interstellar radiative cooling and UV heating
!
      if (linterstellar) &
        call calc_heat_cool_interstellar(df,rho1,TT,TT1,yH)
        !
!  possibility of entropy relaxation in exterior region
!
      if (tau_ss_exterior/=0.) call calc_tau_ss_exterior(f,df)
!
!
!
      if (lspecial) call special_calc_entropy(f,df,uu,glnrho,divu,rho1,lnrho,cs2,TT1)
!
!  Calculate entropy related diagnostics
!
      if(ldiagnos) then
        rho=exp(lnrho)
        if(i_eth/=0) then
          call sum_mn_name(rho*ee,i_eth)
        endif
        if(i_ethtot/=0) then
          call integrate_mn_name(rho*ee,i_ethtot)
        endif
        if(i_TTm/=0) call sum_mn_name(1./TT1,i_TTm)
        if(i_yHm/=0) call sum_mn_name(yH,i_yHm)
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
!  shock entropy diffusivity
!
      if(chi_shock/=0.) then
        if(lroot.and.ip<16) print*, &
                      'calc_heatcond_constchi: use shock diffusion'
        call dot_mn(glnP,gss,g2)
        thdiff = thdiff + chi_t*(del2ss+g2)
      endif
!
!  check maximum diffusion from thermal diffusion
!  With heat conduction, the second-order term for entropy is
!  gamma*chi*del2ss
!
      if (lfirst.and.ldt) maxdiffus=amax1(maxdiffus,(gamma*chi+chi_t))
!
      if(ip==0) print*,rho1 !(to keep compiler quiet)
    endsubroutine calc_heatcond_constchi
!***********************************************************************
    subroutine calc_heatcond_shock(f,df,rho1,glnrho,gss)
!
!  Adds in shock entropy diffusion. There is potential for
!  recycling some quantities from previous calculations.
!
!  20-jul-03/axel: adapted from calc_heatcond_constchi
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
      if(headtt) print*,'calc_heatcond_shock: chi_shock==',chi_shock
!
!  Heat conduction
!  Note: these routines require revision when ionization turned on
!
      call del2(f,iss,del2ss)
      call del2(f,ilnrho,del2lnrho)
      glnT = gamma*gss + gamma1*glnrho
      glnP = gamma*gss + gamma*glnrho
      call dot_mn(glnP,glnT,g2)
!
!  shock entropy diffusivity
!XX
      if(chi_shock/=0.) then
        if(headtt) print*,'calc_heatcond_shock: use shock diffusion'
        call dot_mn(glnP,gss,g2)
        thdiff = (chi_t+chi_shock*f(l1:l2,m,n,ishock))*(del2ss+g2)
      endif
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
!  NEED TO FIX THIS
      if (lfirst.and.ldt) maxdiffus=amax1(maxdiffus,(gamma*chi+chi_t))
!
      if(ip==0) print*,rho1 !(to keep compiler quiet)
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
      if (lfirst.and.ldt) maxdiffus=amax1(maxdiffus,(gamma*chix+chi_t))
!
    endsubroutine calc_heatcond_simple
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
!  identifier
!
      if(headtt) print*,'calc_heatcond: lgravz=',lgravz
!
!  Heat conduction / entropy diffusion
!
      if(headtt) then
        print*,'calc_heatcond: hcond0=',hcond0
        if (lgravz) print*,'calc_heatcond: Fbot,Ftop=',Fbot,Ftop
      endif

      if ((hcond0 /= 0) .or. (chi_t /= 0)) then
        call del2(f,iss,del2ss)
      endif
      if (hcond0 /= 0) then
        if (lgravz) then
          ! For vertical geometry, we only need to calculate this for each
          ! new value of z -> speedup by about 8% at 32x32x64
          if (z_mn(1) /= z_prev) then
            call heatcond(x_mn,y_mn,z_mn,hcond)
            call gradloghcond(x_mn,y_mn,z_mn, glhc)
            z_prev = z_mn(1)
          endif
        endif
        if (lgravr) then
          call heatcond(x_mn,y_mn,z_mn,hcond)
          call gradloghcond(x_mn,y_mn,z_mn, glhc)
        endif
        call del2(f,ilnrho,del2lnrho)
        chix = rho1*hcond
        glnT = gamma*gss + gamma1*glnrho ! grad ln(T)
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
      if (lfirst.and.ldt) maxdiffus=amax1(maxdiffus,(gamma*chix+chi_t))
!--   if (headtt) print*,'calc_heatcond: maxdiffus=',maxdiffus
!
    endsubroutine calc_heatcond
!***********************************************************************
    subroutine calc_heat_cool(f,df,rho1,cs2,ss,TT,TT1)
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
      real, dimension (nx) :: rho1,cs2,ss,TT,TT1
      real, dimension (nx) :: heat,prof
      real :: ssref,zbot,ztop,TTref,profile_buffer
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
      if (lgravz) then
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
        ssref = ss0 + (-alog(gamma) + alog(cs20))/gamma + grads0*ztop
        prof = spread(exp(-0.5*((ztop-z(n))/wcool)**2), 1, l2-l1+1)
        heat = heat - cool*prof*rho1*(cs2-cs20)/cs20
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
        case ('cs2', 'Temp')    ! cooling to reference temperatur cs2cool
          heat = heat - cool*prof*rho1*(cs2-cs2cool)/cs2cool
        case ('entropy')        ! cooling to reference entropy (currently =0)
          heat = heat - cool*prof*(f(l1:l2,m,n,iss)-0.)
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
!  assume a linearly increasing reference profile, TTref
!  This 1/rho1 business is clumpsy, but so would be obvious alternatives...
!
      if(tauheat_coronal/=0.) then
        TTref=(z(n)-ztop)/(zheat_coronal-ztop)*TTheat_coronal
        TTref=(z(n)-zheat_coronal)/(ztop-zheat_coronal)*TTheat_coronal
        heat=heat+amax1(0.,ss*(TTref-TT)/(rho1*tauheat_coronal))
      endif
!
!  add heating and cooling to a reference temperature in a buffer
!  zone at the z boundaries. Only regions in |z| > zheat_buffer are affected.
!  Inverse width of the transition is given by dheat_buffer1.
!
      if(tauheat_buffer/=0.) then
        profile_buffer=0.5*(1.+tanh(dheat_buffer1*(abs(z(n))-zheat_buffer)))
        heat=heat+profile_buffer*ss*(TTheat_buffer-TT)/(rho1*tauheat_buffer)
      endif
!
!  add to entropy equation
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + TT1*rho1*heat
!
    endsubroutine calc_heat_cool
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
    subroutine rprint_entropy(lreset)
!
!  reads and registers print parameters relevant to entropy
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Cdata
      use Sub
!
      integer :: iname
      logical :: lreset
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_eth=0; i_TTm=0; i_yHm=0; i_ssm=0; i_ugradpm=0; i_ethtot=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ethtot',i_ethtot)
        call parse_name(iname,cname(iname),cform(iname),'eth',i_eth)
        call parse_name(iname,cname(iname),cform(iname),'TTm',i_TTm)
        call parse_name(iname,cname(iname),cform(iname),'yHm',i_yHm)
        call parse_name(iname,cname(iname),cform(iname),'ssm',i_ssm)
        call parse_name(iname,cname(iname),cform(iname),'ugradpm',i_ugradpm)
      enddo
!
!  write column where which magnetic variable is stored
!
      write(3,*) 'i_ethtot=',i_ethtot
      write(3,*) 'i_eth=',i_eth
      write(3,*) 'i_TTm=',i_TTm
      write(3,*) 'i_yHm=',i_yHm
      write(3,*) 'i_ssm=',i_ssm
      write(3,*) 'i_ugradpm=',i_ugradpm
      write(3,*) 'nname=',nname
      write(3,*) 'iss=',iss
!
    endsubroutine rprint_entropy
!***********************************************************************
    subroutine heatcond(x,y,z,hcond)
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
      use Cdata, only: nx,lgravz,lgravr
      use Sub, only: step
      use Gravity
!
      real, dimension (nx) :: x,y,z
      real, dimension (nx) :: hcond
!
      if (lgravz) then
        if (lmultilayer) then
          hcond = 1 + (hcond1-1)*step(z,z1,-widthss) &
                    + (hcond2-1)*step(z,z2,widthss)
          hcond = hcond0*hcond
        else
          hcond=Kbot
        endif
      endif

      if (lgravr) then
        hcond = hcond0
      endif
!
      if(ip==0) print*,x,y  !(to keep compiler quiet)
    endsubroutine heatcond
!***********************************************************************
    subroutine gradloghcond(x,y,z,glhc)
!
!  calculate grad(log hcond), where hcond is the heat conductivity
!  NB: *Must* be in sync with heatcond() above.
!  23-jan-2002/wolf: coded
!
      use Cdata, only: nx,lgravz,lgravr
      use Sub, only: der_step
      use Gravity
!
      real, dimension (nx) :: x,y,z
      real, dimension (nx,3) :: glhc
!
      if (lgravz) then
        glhc(:,1:2) = 0.
        glhc(:,3) = (hcond1-1)*der_step(z,z1,-widthss) &
                    + (hcond2-1)*der_step(z,z2,widthss)
        glhc(:,3) = hcond0*glhc(:,3)
      endif

      if (lgravr) then
        glhc = 0.
      endif
!
      if(ip==0) print*,x,y  !(to keep compiler quiet)
    endsubroutine gradloghcond
!***********************************************************************
endmodule Entropy

