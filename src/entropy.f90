! $Id: entropy.f90,v 1.149 2003-04-05 19:05:07 brandenb Exp $

!  This module takes care of entropy (initial condition
!  and time advance)

module Entropy

  use Cparam
  use Cdata
  use Hydro
  use Interstellar
  use Viscosity

  implicit none

  real, dimension (nx) :: cs2,TT1
  real :: radius_ss=0.1,ampl_ss=0.,widthss=2*epsi,epsilon_ss
  real :: luminosity=0.,wheat=0.1,cs2cool=0.,cool=0.,rcool=1.,wcool=0.1
  real :: ss_left,ss_right,chi=0.,chi_t=0.,ss0=0.,khor_ss=1.,ss_const=0.
  real :: tau_ss_exterior=0.
  !parameters for Sedov type initial condition
  real :: thermal_background=0., thermal_peak=0., thermal_scaling=1.
  real :: hcond0=0.
  real :: Fbot=impossible,hcond1=impossible,hcond2=impossible
  real :: FbotKbot=impossible,Kbot=impossible
  logical :: lcalc_heatcond_simple=.false.,lmultilayer=.true.
  logical :: lcalc_heatcond_constchi=.false.
  character (len=labellen) :: initss='nothing',pertss='zero',cooltype='Temp'

  ! input parameters
  namelist /entropy_init_pars/ &
       initss,pertss,grads0,radius_ss,ampl_ss,widthss,epsilon_ss, &
       ss_left,ss_right,ss_const,mpoly0,mpoly1,mpoly2,isothtop, &
       khor_ss, thermal_background, thermal_peak, thermal_scaling

  ! run parameters
  namelist /entropy_run_pars/ &
       hcond0,hcond1,hcond2,widthss, &
       luminosity,wheat,cooltype,cool,cs2cool,rcool,wcool,Fbot, &
       chi_t,lcalc_heatcond_simple,tau_ss_exterior, &
       chi,lcalc_heatcond_constchi,lmultilayer,Kbot

  ! other variables (needs to be consistent with reset list below)
  integer :: i_ssm=0,i_ugradpm=0

  contains

!***********************************************************************
    subroutine register_entropy()
!
!  initialise variables which should know that we solve an entropy
!  equation: ient, etc; increase nvar accordingly
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
      ient = nvar+1             ! index to access entropy
      nvar = nvar+1
!
      if ((ip<=8) .and. lroot) then
        print*, 'Register_ent:  nvar = ', nvar
        print*, 'ient = ', ient
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: entropy.f90,v 1.149 2003-04-05 19:05:07 brandenb Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_ent: nvar > mvar')
      endif
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
            if (bcz1(ient)=='c1') then
              Fbot=-gamma/(gamma-1)*hcond0*gravz/(mpoly0+1)
              if (lroot) print*, 'Calculated Fbot = ', Fbot
            else
              Fbot=0.
            endif
          endif
          FbotKbot=Fbot/(hcond0*hcond1)
        else
          !
          !  Wolfgang, in future we should define chiz=chi(z) or Kz=K(z) here.
          !  calculate hcond and FbotKbot=Fbot/K, where K=hcond is radiative conductivity
          !
          Kbot=gamma1/gamma*(mpoly+1.)*Fbot
          FbotKbot=gamma/gamma1/(mpoly+1.)
          if(lroot) print*,'ss_run_hook: Fbot,Kbot=',Fbot,Kbot
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
      use IO
      use Sub
      use Gravity
      use Initcond
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: tmp,xx,yy,zz
      real :: cs2int,ssint
!
      intent(in) :: xx,yy,zz
      intent(inout) :: f
!
      select case(initss)
        case('nothing'); if(lroot) print*,'init_ss: nothing'
        case('zero', '0'); f(:,:,:,ient) = 0.
        case('const_ss'); f(:,:,:,ient) = ss_const
        case('blob'); call blob(ampl_ss,f,ient,radius_ss,0.,0.,0.)
        case('isothermal'); if(lroot) print*,'init_ss: isotherm set in density'
        case('Ferriere'); if(lroot) print*,'init_ss: Ferriere set in density'
        case('xjump'); call jump(f,ient,ss_left,ss_right,widthss,'x')
        case('hor-fluxtube'); call htube(ampl_ss,f,ient,ient,xx,yy,zz,radius_ss,epsilon_ss)
        case('hor-tube'); call htube2(ampl_ss,f,ient,ient,xx,yy,zz,radius_ss,epsilon_ss)

      case('sedov') 
        if (lroot) print*,'init_ss: sedov - thermal background with gaussian energy burst'
         call blob(thermal_peak,f,ient,radius_ss,0.,0.,0.)
         f(:,:,:,ient) = (alog(f(:,:,:,ient) + thermal_background)+alog(thermal_scaling))/gamma 
   
      case('isobaric')
        !
        !  ss = - ln(rho/rho0)
        !
        if (lroot) print*,'init_ss: isobaric stratification'
        f(:,:,:,ient) = -(f(:,:,:,ilnrho)-lnrho0)

      case('isentropic', '1')
        !
        !  ss = const.
        !
        if (lroot) print*,'isentropic stratification'
        ! ss0=alog(-gamma1*gravz*zinfty)/gamma
        ! print*,'isentropic stratification; ss=',ss0
        f(:,:,:,ient)=0.
        if (ampl_ss/=0.) then
          print*,'put bubble: radius_ss,ampl_ss=',radius_ss,ampl_ss
          tmp=xx**2+yy**2+zz**2
          f(:,:,:,ient)=f(:,:,:,ient)+ampl_ss*exp(-tmp/amax1(radius_ss**2-tmp,1e-20))
          !f(:,:,:,ient)=f(:,:,:,ient)+ampl_ss*exp(-tmp/radius_ss**2)
        endif

      case('linprof', '2')
        !
        !  linear profile of ss, centered around ss=0.
        !
        if (lroot) print*,'linear entropy profile'
        f(:,:,:,ient) = grads0*zz

      case('piecew-poly', '4')
        !
        !  piecewise polytropic convection setup
        !  cs0, rho0 and ss0=0 refer to height z=zref
        !
        if (lroot) print*,'piecewise polytropic vertical stratification (ss)'
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
        f(:,:,:,ient) = 0.    ! just in case
        ! top layer
        call polytropic_ss_z(f,mpoly2,zz,tmp,zref,z2,z0+2*Lz, &
                             isothtop,cs2int,ssint)
        ! unstable layer
        call polytropic_ss_z(f,mpoly0,zz,tmp,z2,z1,z2,0,cs2int,ssint)
        ! stable layer
        call polytropic_ss_z(f,mpoly1,zz,tmp,z1,z0,z1,0,cs2int,ssint)

      case('polytropic', '5')
        !
        !  polytropic stratification
        !  cs0, rho0 and ss0=0 refer to height z=zref
        !
        if (lroot) print*,'polytropic vertical stratification (ss)'
        !
        cs20 = cs0**2
        ss0 = 0.              ! reference value ss0 is zero
        f(:,:,:,ient) = ss0   ! just in case
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
        if (lroot) print*,'No such value for initss: ', trim(initss)
        call stop_it("")

      endselect

!      endif
!
      if (lgravr) then
          f(:,:,:,ient) = -0.
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
        if (lroot) print*,'adding hexagonal perturbation to ss'
        f(:,:,:,ient) = f(:,:,:,ient) &
                        + ampl_ss*(2*cos(sqrt(3.)*0.5*khor_ss*xx) &
                                    *cos(0.5*khor_ss*yy) &
                                   + cos(khor_ss*yy) &
                                  ) * cos(pi*zz)

      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*,'No such value for pertss:', pertss
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
      real, dimension (mx,my,mz,mvar) :: f
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
      f(:,:,:,ient) = p*f(:,:,:,ient)  + (1-p)*tmp
!
    endsubroutine polytropic_ss_z
!***********************************************************************
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
      use Slices
      use IO
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: uu,glnrho,gss
      real, dimension (nx) :: ugss,uglnrho, divu
      real, dimension (nx) :: lnrho,ss,rho1,cs2,TT1,cp1tilde
      integer :: j,ju
!
      intent(in) :: f,uu,glnrho,rho1,lnrho
      intent(out) :: df,cs2,TT1
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'SOLVE dss_dt'
      if (headtt) call identify_bcs('ss',ient)
!
!  entropy gradient: needed for advection and pressure gradient
!
      call grad(f,ient,gss)
!
!  sound speed squared
!  include in maximum advection speed (for timestep)
!
      ss=f(l1:l2,m,n,ient)
!
!  calculate cs2, TT1, and cp1tilde in a separate routine
!  With IONIZATION=noionization, assume perfect gas with const coeffs
!
      call thermodynamics(lnrho,ss,cs2,TT1,cp1tilde)
      if (headtt) print*,'dss_dt: cs2,TT,cp1tilde=',cs2(1),1./TT1(1),cp1tilde(1)
!
!  use sound speed in Courant condition
!
      if (lfirst.and.ldt) maxadvec2=amax1(maxadvec2,cs2)
      if (ip<8.and.lroot.and.imn==1) print*,'maxadvec2,cs2=',maxadvec2,cs2
      if (headtt) print*,'entropy (but not used with ionization): cs20=',cs20
!
!  subtract pressure gradient term in momentum equation
!
      if (lhydro) then
        do j=1,3
          ju=j+iuu-1
          df(l1:l2,m,n,ju)=df(l1:l2,m,n,ju)-cs2*(glnrho(:,j)+cp1tilde*gss(:,j))
        enddo
      endif
if(m==4.and.n>250) print*,'glnrho(4,3),cp1tilde(4),gss(4,3)=',glnrho(4,3),cp1tilde(4),gss(4,3)
!
!  advection term
!
      call dot_mn(uu,gss,ugss)
      df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) - ugss
!
!  1/T is now calculated in thermodynamics
!  Viscous heating depends on ivisc; no visc heating if ivisc='simplified'
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
!  heating/cooling
!
      if ((luminosity /= 0) .or. (cool /= 0)) &
           call calc_heat_cool(f,df,rho1,cs2,TT1)
!
!ngrs: switch off for debug
!      if (linterstellar) &
!        call calc_heat_cool_interstellar(df,rho1,TT1)
!
!  possibility of entropy relaxation in exterior region
!
      if (tau_ss_exterior/=0.) call calc_tau_ss_exterior(f,df)
!
!  Calculate entropy related diagnostics
!
      if (ldiagnos) then
        if (i_ssm/=0) call sum_mn_name(ss,i_ssm)
        if (i_ugradpm/=0) then
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
      real, dimension (mx,my,mz,mvar) :: f,df
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
      call del2(f,ient,del2ss)
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
      df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + thdiff
      if (headtt) print*,'calc_heatcond_simple: added thdiff'
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
      real, dimension (mx,my,mz,mvar) :: f,df
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
      call del2(f,ient,del2ss)
      call del2(f,ilnrho,del2lnrho)
      chix = rho1*hcond
      glnT = gamma*gss + gamma1*glnrho ! grad ln(T)
      glnThcond = glnT !... + glhc/spread(hcond,2,3)    ! grad ln(T*hcond)
      call dot_mn(glnT,glnThcond,g2)
      thdiff = chix * (gamma*del2ss+gamma1*del2lnrho + g2)
!
!  add heat conduction to entropy equation
!
      df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + thdiff
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
      real, dimension (mx,my,mz,mvar) :: f,df
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
        if (lgravz) print*,'Fbot=',Fbot
      endif

      if ((hcond0 /= 0) .or. (chi_t /= 0)) then
        call del2(f,ient,del2ss)
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
          print*,'"turbulent" entropy diffusion: chi_t=',chi_t
          if (hcond0 /= 0) then
            print*,"WARNING: hcond0 and chi_t combined don't seem to make sense"
          endif
        endif
!        df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient)+chi_t*del2ss
        thdiff = chi_t*del2ss
      endif
!
!  check for NaNs initially
!
      if (headt .and. (hcond0 /= 0)) then
        if (notanumber(glhc))      print*,'NaNs in glhc'
        if (notanumber(rho1))      print*,'NaNs in rho1'
        if (notanumber(hcond))     print*,'NaNs in hcond'
        if (notanumber(chix))       print*,'NaNs in chix'
        if (notanumber(del2ss))    print*,'NaNs in del2ss'
        if (notanumber(del2lnrho)) print*,'NaNs in del2lnrho'
        if (notanumber(glhc))      print*,'NaNs in glhc'
        if (notanumber(1/hcond))   print*,'NaNs in 1/hcond'
        if (notanumber(glnT))      print*,'NaNs in glnT'
        if (notanumber(glnThcond)) print*,'NaNs in glnThcond'
        if (notanumber(g2))        print*,'NaNs in g2'
        if (notanumber(thdiff))    print*,'NaNs in thdiff'
        !
        !  most of these should trigger the following trap
        !
        if (notanumber(thdiff)) then

print*, 'm,n,y(m),z(n)=',m,n,y(m),z(n)
call stop_it('NaNs in thdiff')
endif
      endif

      if (headt .and. lfirst .and. ip<=9) then
        call output_pencil(trim(directory)//'/chi.dat',chix,1)
        call output_pencil(trim(directory)//'/hcond.dat',hcond,1)
        call output_pencil(trim(directory)//'/glhc.dat',glhc,3)
      endif
      df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + thdiff
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
    subroutine calc_heat_cool(f,df,rho1,cs2,TT1)
!
!  heating and cooling
!
!  02-jul-02/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use Gravity
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx) :: rho1,cs2,TT1
      real, dimension (nx) :: heat,prof
      real :: ssref,zbot,ztop
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
          heat = heat - cool*prof*(f(l1:l2,m,n,ient)-0.)
        case default
          if (lroot) print*,'No such value for cooltype: ', trim(cooltype)
          call stop_it("")
        endselect
      endif
!
      df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + TT1*rho1*heat
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
      real, dimension (mx,my,mz,mvar) :: f,df
      real :: scl
!
      intent(in) :: f
      intent(out) :: df
!
      if (headtt) print*,'calc_tau_ss_exterior: tau=',tau_ss_exterior
      if(z(n)>zgrav) then
        scl=1./tau_ss_exterior
        df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient)-scl*f(l1:l2,m,n,ient)
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
        i_ssm=0; i_ugradpm=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ssm',i_ssm)
        call parse_name(iname,cname(iname),cform(iname),'ugradpm',i_ugradpm)
      enddo
!
!  write column where which magnetic variable is stored
!
      write(3,*) 'i_ssm=',i_ssm
      write(3,*) 'i_ugradpm=',i_ugradpm
      write(3,*) 'nname=',nname
      write(3,*) 'ient=',ient
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
    subroutine bc_ss_flux(f,topbot)
!
!  constant flux boundary condition for entropy (called when bcz='c1')
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my) :: tmp_xy,cs2_xy,rho_xy
      integer :: i
!
      if(ldebug) print*,'ENTER: bc_ss, cs20,cs0=',cs20,cs0
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('strange-bot')
        if(headtt) print*,'bc_ss_flux: hcond0,hcond1=',hcond0,hcond1
        if (bcz1(ilnrho) /= "a2") &
             call stop_it("BOUNDCONDS: Inconsistent boundary conditions 1.")
        tmp_xy = gamma1/cs20 & ! 1/T_0 (i.e. 1/T at boundary)
                 * exp(-gamma*f(:,:,n1,ient) &
                       - gamma1*(f(:,:,n1,ilnrho)-lnrho0))
        tmp_xy = Fbot/(hcond0*hcond1) * tmp_xy ! F_heat/(hcond T_0)
        do i=1,nghost
          f(:,:,n1-i,ient) = &
               (2*i*dz*tmp_xy &
                + 2*gamma1*(f(:,:,n1+i,ilnrho)-f(:,:,n1,ilnrho)) &
               )/gamma &
               + f(:,:,n1+i,ient)
        enddo
!
!  bottom boundary
!  ===============
!
      case('bot')
        if(headt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0*hcond1
!       if(bcz1(ilnrho)/="a2") call stop_it("bc_ss_flux: bad lnrho bc")
!
!  calculate Fbot/(K*cs2)
!
        rho_xy=exp(f(:,:,n1,ilnrho))
        cs2_xy=cs20*exp(gamma1*(f(:,:,n1,ilnrho)-lnrho0)+gamma*f(:,:,n1,ient))
!
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy. 
!
        if(lcalc_heatcond_constchi) then
          tmp_xy=Fbot/(rho_xy*chi*cs2_xy)
        else
          tmp_xy=FbotKbot/cs2_xy
        endif
!
!  enforce ds/dz + gamma1/gamma*dlnrho/dz = - gamma1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n1-i,ient)=f(:,:,n1+i,ient)+gamma1/gamma* &
              (f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho)+2*i*dz*tmp_xy)
        enddo
!
!  top boundary
!  ============
!
      case('top')
        if(headtt) print*,'bc_ss_flux: hcond0=',hcond0
        if (bcz2(ilnrho) /= "a2") &
             call stop_it("BOUNDCONDS: Inconsistent boundary conditions 2.")
        tmp_xy = gamma1/cs20 & ! 1/T_0 (i.e. 1/T at boundary)
                 * exp(-gamma*f(:,:,n2,ient) &
                       - gamma1*(f(:,:,n2,ilnrho)-lnrho0))
        tmp_xy = FbotKbot * tmp_xy ! F_heat/(hcond T_0)
        do i=1,nghost
          f(:,:,n2+i,ient) = &
               (-2*i*dz*tmp_xy &
                + 2*gamma1*(f(:,:,n2-i,ilnrho)-f(:,:,n2,ilnrho)) &
               )/gamma &
               + f(:,:,n2-i,ient)
        enddo
      case default
        if(lroot) print*,"invalid argument for 'bc_ss_flux'"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_flux
!***********************************************************************
    subroutine bc_ss_temp_old(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my) :: tmp_xy
      integer :: i
!
      if(ldebug) print*,'ENTER: bc_ss, cs20,cs0=',cs20,cs0
!
!  Do the `c2' boundary condition (fixed temperature/sound speed) for entropy.
!  This assumes that the density is already set (ie density must register
!  first!)
!  tmp_xy = s(x,y) on the boundary.
!  gamma*s/cp = [ln(cs2/cs20)-(gamma-1)ln(rho/rho0)]
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (ldebug) print*,'set bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2bot<=0'
        if (bcz1(ilnrho) /= "a2") &
             call stop_it("BOUNDCONDS: Inconsistent boundary conditions 3.")
        tmp_xy = (-gamma1*(f(:,:,n1,ilnrho)-lnrho0) &
                 + alog(cs2bot/cs20)) / gamma
        f(:,:,n1,ient) = tmp_xy
        do i=1,nghost
          f(:,:,n1-i,ient) = 2*tmp_xy - f(:,:,n1+i,ient)
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*,'set top temperature: cs2top=',cs2top
        if (cs2top<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2top<=0'
!       if (bcz1(ilnrho) /= "a2") &
!            call stop_it("BOUNDCONDS: Inconsistent boundary conditions 4.")
        tmp_xy = (-gamma1*(f(:,:,n2,ilnrho)-lnrho0) &
                 + alog(cs2top/cs20)) / gamma
        f(:,:,n2,ient) = tmp_xy
        do i=1,nghost
          f(:,:,n2+i,ient) = 2*tmp_xy - f(:,:,n2-i,ient)
        enddo
      case default
        if(lroot) print*,"invalid argument for 'bc_ss_flux'"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_temp_old
!***********************************************************************
    subroutine bc_ss_temp_x(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  3-aug-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar) :: f
      real :: tmp
      integer :: i
!
      if(ldebug) print*,'ENTER: bc_ss_temp_x, cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (ldebug) print*,'set x bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2bot<=0'
        tmp = 2/gamma*alog(cs2bot/cs20)
        f(l1,:,:,ient) = 0.5*tmp - gamma1/gamma*(f(l1,:,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(l1-i,:,:,ient) = -f(l1+i,:,:,ient) + tmp &
               - gamma1/gamma*(f(l1+i,:,:,ilnrho)+f(l1-i,:,:,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*,'set x top temperature: cs2top=',cs2top
        if (cs2top<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
        f(l2,:,:,ient) = 0.5*tmp - gamma1/gamma*(f(l2,:,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(l2+i,:,:,ient) = -f(l2-i,:,:,ient) + tmp &
               - gamma1/gamma*(f(l2-i,:,:,ilnrho)+f(l2+i,:,:,ilnrho)-2*lnrho0)
        enddo

      case default
        if(lroot) print*,"invalid argument for 'bc_ss_temp_x'"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_temp_x
!***********************************************************************
    subroutine bc_ss_temp_y(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  3-aug-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar) :: f
      real :: tmp
      integer :: i
!
      if(ldebug) print*,'ENTER: bc_ss_temp_y, cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (ldebug) print*,'set y bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2bot<=0'
        tmp = 2/gamma*alog(cs2bot/cs20)
        f(:,m1,:,ient) = 0.5*tmp - gamma1/gamma*(f(:,m1,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,m1-i,:,ient) = -f(:,m1+i,:,ient) + tmp &
               - gamma1/gamma*(f(:,m1+i,:,ilnrho)+f(:,m1-i,:,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*,'set y top temperature: cs2top=',cs2top
        if (cs2top<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
        f(:,m2,:,ient) = 0.5*tmp - gamma1/gamma*(f(:,m2,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,m2+i,:,ient) = -f(:,m2-i,:,ient) + tmp &
               - gamma1/gamma*(f(:,m2-i,:,ilnrho)+f(:,m2+i,:,ilnrho)-2*lnrho0)
        enddo

      case default
        if(lroot) print*,"invalid argument for 'bc_ss_temp_y'"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_temp_y
!***********************************************************************
    subroutine bc_ss_temp_z(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  3-aug-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar) :: f
      real :: tmp
      integer :: i
!
      if(ldebug) print*,'ENTER: bc_ss_temp_z, cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (ldebug) print*,'set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2bot<=0'
        tmp = 2/gamma*alog(cs2bot/cs20)
        f(:,:,n1,ient) = 0.5*tmp - gamma1/gamma*(f(:,:,n1,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,:,n1-i,ient) = -f(:,:,n1+i,ient) + tmp &
               - gamma1/gamma*(f(:,:,n1+i,ilnrho)+f(:,:,n1-i,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*,'set z top temperature: cs2top=',cs2top
        if (cs2top<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
        f(:,:,n2,ient) = 0.5*tmp - gamma1/gamma*(f(:,:,n2,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,:,n2+i,ient) = -f(:,:,n2-i,ient) + tmp &
               - gamma1/gamma*(f(:,:,n2-i,ilnrho)+f(:,:,n2+i,ilnrho)-2*lnrho0)
        enddo

      case default
        if(lroot) print*,"invalid argument for 'bc_ss_temp_z'"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_temp_z
!***********************************************************************
    subroutine bc_ss_stemp_x(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!  3-aug-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar) :: f
      integer :: i
!
      if(ldebug) print*,'ENTER: bc_ss_stemp_x, cs20,cs0=',cs20,cs0
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (cs2bot<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2bot<=0'
        do i=1,nghost
          f(l1-i,:,:,ient) = f(l1+i,:,:,ient) &
               + gamma1/gamma*(f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho))
        enddo
!
!  top boundary
!
      case('top')
        if (cs2top<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2top<=0'
        do i=1,nghost
          f(l2+i,:,:,ient) = f(l2-i,:,:,ient) &
               + gamma1/gamma*(f(l2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho))
        enddo

      case default
        if(lroot) print*,"invalid argument for 'bc_ss_stemp_x'"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_stemp_x
!***********************************************************************
    subroutine bc_ss_stemp_y(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!  3-aug-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar) :: f
      integer :: i
!
      if(ldebug) print*,'ENTER: bc_ss_stemp_y, cs20,cs0=',cs20,cs0
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (cs2bot<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2bot<=0'
        do i=1,nghost
          f(:,m1-i,:,ient) = f(:,m1+i,:,ient) &
               + gamma1/gamma*(f(:,m1+i,:,ilnrho)-f(:,m1-i,:,ilnrho))
        enddo
!
!  top boundary
!
      case('top')
        if (cs2top<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2top<=0'
        do i=1,nghost
          f(:,m2+i,:,ient) = f(:,m2-i,:,ient) &
               + gamma1/gamma*(f(:,m2-i,:,ilnrho)-f(:,m2+i,:,ilnrho))
        enddo

      case default
        if(lroot) print*,"invalid argument for 'bc_ss_stemp_y'"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_stemp_y
!***********************************************************************
    subroutine bc_ss_stemp_z(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!  3-aug-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar) :: f
      integer :: i
!
      if(ldebug) print*,'ENTER: bc_ss_stemp_x, cs20,cs0=',cs20,cs0
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (cs2bot<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2bot<=0'
        do i=1,nghost
          f(:,:,n1-i,ient) = f(:,:,n1+i,ient) &
               + gamma1/gamma*(f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho))
        enddo
!
!  top boundary
!
      case('top')
        if (cs2top<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2top<=0'
        do i=1,nghost
          f(:,:,n2+i,ient) = f(:,:,n2-i,ient) &
               + gamma1/gamma*(f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho))
        enddo

      case default
        if(lroot) print*,"invalid argument for 'bc_ss_stemp_z'"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_stemp_z
!***********************************************************************
    subroutine bc_ss_energy(f,topbot)
!
!  boundary condition for entropy
!
!  may-2002/nils: coded
!  11-jul-2002/nils: moved into the entropy module
!
      use Mpicomm, only: stop_it
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my) :: cs2_2d
      integer :: i
!
!  The 'ce' boundary condition for entropy makes the energy constant at
!  the boundaries.
!  This assumes that the density is already set (ie density must register
!  first!)
!
    select case(topbot)
!
! Bottom boundary
!
    case('bot')
      !  Set cs2 (temperature) in the ghost points to the value on
      !  the boundary
      !
      cs2_2d=cs20*exp(gamma1*f(:,:,n1,ilnrho)+gamma*f(:,:,n1,ient))
      do i=1,nghost
         f(:,:,n1-i,ient)=1./gamma*(-gamma1*f(:,:,n1-i,ilnrho)-log(cs20)&
              +log(cs2_2d))
      enddo

!
! Top boundary
!
    case('top')
      !  Set cs2 (temperature) in the ghost points to the value on
      !  the boundary
      !
      cs2_2d=cs20*exp(gamma1*f(:,:,n2,ilnrho)+gamma*f(:,:,n2,ient))
      do i=1,nghost
         f(:,:,n2+i,ient)=1./gamma*(-gamma1*f(:,:,n2+i,ilnrho)-log(cs20)&
              +log(cs2_2d))
      enddo
    case default
       if(lroot) print*,"invalid argument for 'bc_ss_flux'"
        call stop_it("")
    endselect

    end subroutine bc_ss_energy
!***********************************************************************
endmodule Entropy
