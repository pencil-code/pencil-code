! $Id: entropy.f90,v 1.458 2006-11-20 03:27:23 dobler Exp $


!  This module takes care of entropy (initial condition
!  and time advance)
 
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lentropy = .true.
! CPARAM logical, parameter :: ltemperature = .false.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED ugss,Ma2,fpres
!
!***************************************************************
module Entropy

  use Cparam
  use Cdata
  use Messages
  use Interstellar
  use Viscosity
  use EquationOfState, only: gamma, gamma1, gamma11, cs20, cs2top, cs2bot, &
                         isothtop, mpoly0, mpoly1, mpoly2, cs2cool, &
                         beta_glnrho_global

  implicit none

  include 'entropy.h'

  !real, dimension (nx) :: cs2,TT1
  real :: radius_ss=0.1,ampl_ss=0.,widthss=2*epsi,epsilon_ss=0.
  real :: luminosity=0.,wheat=0.1,cool=0.,rcool=0.,wcool=0.1
  real :: TT_int,TT_ext,cs2_int,cs2_ext,cool_int=0.,cool_ext=0.,ampl_TT=0.
  real :: chi=0.,chi_t=0.,chi_shock=0.,chi_hyper3=0.
  real :: Kgperp=0.,Kgpara=0.,tdown=1.,allp=2.
  real :: ss_left=1.,ss_right=1.
  real :: ss0=0.,khor_ss=1.,ss_const=0.
  real :: pp_const=0.
  real :: tau_ss_exterior=0.,T0=1.
  real :: mixinglength_flux=0.
  !parameters for Sedov type initial condition
  real :: center1_x=0., center1_y=0., center1_z=0.
  real :: center2_x=0., center2_y=0., center2_z=0.
  real :: kx_ss=1.
  real :: thermal_background=0., thermal_peak=0., thermal_scaling=1.
  real :: hcond0=impossible
  real :: hcond1=impossible,hcond2=impossible
  real :: Fbot=impossible,FbotKbot=impossible,Kbot=impossible
  real :: Ftop=impossible,FtopKtop=impossible,Ktop=impossible
  real :: chit_prof1=1.,chit_prof2=1.
  real :: tau_cor=0.,TT_cor=0.,z_cor=0.
  real :: tauheat_buffer=0.,TTheat_buffer=0.,zheat_buffer=0.,dheat_buffer1=0.
  real :: heat_uniform=0.,cool_RTV=0.
  real :: deltaT_poleq=0.,beta_hand=1.,r_bcz=0.
  integer, parameter :: nheatc_max=4
  logical :: lturbulent_heat=.false.
  logical :: lheatc_Kconst=.false.,lheatc_simple=.false.,lheatc_chiconst=.false.
  logical :: lheatc_tensordiffusion=.false.,lheatc_spitzer=.false.
  logical :: lheatc_corona=.false.
  logical :: lheatc_shock=.false.,lheatc_hyper3ss=.false.
  logical :: lupw_ss=.false.,lmultilayer=.true.
  logical :: lpressuregradient_gas=.true.,ladvection_entropy=.true.
  logical :: lviscosity_heat=.true.
  logical :: lfreeze_sint=.false.,lfreeze_sext=.false.
  logical :: lhcond_global=.false.

  integer :: iglobal_hcond=0
  integer :: iglobal_glhc=0

  character (len=labellen), dimension(ninit) :: initss='nothing'
  character (len=labellen) :: borderss='nothing'
  character (len=labellen) :: pertss='zero'
  character (len=labellen) :: cooltype='Temp',cooling_profile='gaussian'
  character (len=labellen), dimension(nheatc_max) :: iheatcond='nothing'
  character (len=4) :: iinit_str
  
  !
  ! Parameters for subroutine cool_RTV in SI units (from Cook et al. 1989)
  ! 
  double precision, parameter, dimension (10) :: & 
       intlnT_1 =(/4.605, 8.959, 9.906, 10.534, 11.283, 12.434, 13.286, 14.541, 17.51, 20.723 /) 
  double precision, parameter, dimension (9) :: &
       lnH_1 = (/ -542.398,  -228.833, -80.245, -101.314, -78.748, -53.88, -80.452, -70.758, -91.182/), &
       B_1   = (/     50.,      15.,      0.,      2.0,      0.,    -2.,      0., -0.6667,    0.5 /)  
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
      initss,     &
!TEST
      pertss,     &
      grads0,     &
      radius_ss,  &
      ampl_ss,    &
      widthss,    &
      epsilon_ss, &
!AB: mixinglength_flux is used as flux in mixing length initial condition
      mixinglength_flux, &
      pp_const, &
      ss_left,ss_right,ss_const,mpoly0,mpoly1,mpoly2,isothtop, &
      khor_ss,thermal_background,thermal_peak,thermal_scaling,cs2cool, &
      center1_x, center1_y, center1_z, center2_x, center2_y, center2_z, &
      T0,ampl_TT,kx_ss,beta_glnrho_global,ladvection_entropy, &
      lviscosity_heat,r_bcz
  ! run parameters
  namelist /entropy_run_pars/ &
      hcond0,hcond1,hcond2,widthss,borderss, &
!AB: allow polytropic indices to be read in also during run stage.
!AB: They are used to re-calculate the radiative conductivity profile.
      mpoly0,mpoly1,mpoly2, &
      luminosity,wheat,cooling_profile,cooltype,cool,cs2cool,rcool,wcool,Fbot, &
      chi_t,chit_prof1,chit_prof2,chi_shock,chi,iheatcond, &
      Kgperp,Kgpara, cool_RTV, &
      tau_ss_exterior,lmultilayer,Kbot,tau_cor,TT_cor,z_cor, &
      tauheat_buffer,TTheat_buffer,zheat_buffer,dheat_buffer1, &
      heat_uniform,lupw_ss,cool_int,cool_ext,chi_hyper3, &
      lturbulent_heat,deltaT_poleq,lpressuregradient_gas, &
      tdown, allp,beta_glnrho_global,ladvection_entropy, &
      lviscosity_heat,r_bcz,lfreeze_sint,lfreeze_sext,lhcond_global

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_dtc=0,idiag_eth=0,idiag_ethdivum=0,idiag_ssm=0
  integer :: idiag_eem=0,idiag_ppm=0,idiag_csm=0
  integer :: idiag_ugradpm=0,idiag_ethtot=0,idiag_dtchi=0,idiag_ssmphi=0
  integer :: idiag_yHm=0,idiag_yHmax=0,idiag_TTm=0,idiag_TTmax=0,idiag_TTmin=0
  integer :: idiag_fconvz=0,idiag_dcoolz=0,idiag_fradz=0,idiag_fturbz=0
  integer :: idiag_ssmz=0,idiag_ssmy=0,idiag_ssmx=0,idiag_TTmz=0,idiag_TTp=0

  contains

!***********************************************************************
    subroutine register_entropy()
!
!  initialise variables which should know that we solve an entropy
!  equation: iss, etc; increase nvar accordingly
!
!  6-nov-01/wolf: coded
!
      use FArrayManager
      use Cdata
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call fatal_error('register_entropy','module registration called twice')
      first = .false.
!
      call farray_register_pde('ss',iss)
!
!
!  identify version number
!
      if (lroot) call cvs_id( &

           "$Id: entropy.f90,v 1.458 2006-11-20 03:27:23 dobler Exp $")
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
      use FArrayManager
      use Gravity, only: gravz,g0
      use EquationOfState, only: cs0, get_soundspeed, get_cp1, &
                                 beta_glnrho_global, beta_glnrho_scaled, &
                                 mpoly, mpoly0, mpoly1, mpoly2, &
                                 select_eos_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      real :: beta1,cp1,lnrho_dummy=0.
      real :: beta0,TT_crit
      integer :: i
      logical :: lnothing
      type (pencil_case) :: p
      real, dimension (nx) :: hcond
      real, dimension (nx,3) :: glhc
!
! Check any module dependencies
!
      if (.not. leos) then
        call fatal_error('initialize_entropy','EOS=noeos but entropy requires an EQUATION OF STATE for the fluid')
      endif
!
! Tell the equation of state that we're here and what f variable we use
!
!ajwm      if (.not. lreloading) then ! already in place when reloading
        if (pretend_lnTT) then
          call select_eos_variable('lnTT',iss)
        else 
          call select_eos_variable('ss',iss)
        endif
!ajwm      endif
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
!   04-mai-2006/boris: hcond0 given by mixinglength_flux in the MLT case
        if (mixinglength_flux /= 0.) then
          hcond0=-mixinglength_flux/(gamma/(gamma-1.)*gravz/(mpoly0+1.))
          Kbot=hcond0
          lmultilayer=.true.
        endif
      else                      ! hcond0 = possible
        if (Kbot == impossible) then
          Kbot = hcond0
        else
          call warning('initialize_entropy','You should not set Kbot and hcond0 at the same time')
        endif
      endif
!
!  freeze entroopy
!
      if (lfreeze_sint) lfreeze_varint(iss) = .true.
      if (lfreeze_sext) lfreeze_varext(iss) = .true.
!
!  make sure the top boundary condition for temperature (if cT is used)
!  knows about the cooling function or vice versa (cs2cool will take over
!  if /=0)
!
      if (lgravz .and. lrun) then
        if (cs2top/=cs2cool) then
          if (lroot) print*,'initialize_entropy: cs2top,cs2cool=',cs2top,cs2cool
          if (cs2cool /= 0.) then ! cs2cool is the value to go for
            if (lroot) print*,'initialize_entropy: now set cs2top=cs2cool'
            cs2top=cs2cool
          else                  ! cs2cool=0, so go for cs2top
            if (lroot) print*,'initialize_entropy: now set cs2cool=cs2top'
            cs2cool=cs2top
          endif
        endif
!
!  settings for fluxes
!
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
              if (lroot) print*,'initialize_entropy: Calculated Fbot,Kbot=', &
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
              if (lroot) print*,'initialize_entropy: Ftop,Ktop=',Ftop,Ktop
            ! else
            !! Don't need Ftop in this case (?)
            !  Ftop=0.
            endif
          endif
!
        endif
      endif
!
!
!     if (mixinglength_flux /= 0.) then
!        Fbot=mixinglength_flux
!        hcond0=-Fbot/(gamma/(gamma-1.)*gravz/(mpoly0+1.))
!        Kbot=gamma1/gamma*(mpoly+1.)*Fbot
!        FbotKbot=gamma/gamma1/(mpoly+1.)
!        print*,'calculated hcond0 from mixinglength_flux: ',hcond0,Fbot
!     endif
!
!   make sure all relevant parameters are set for spherical shell problems
!
      select case(initss(1))
        case('geo-kws','geo-benchmark','shell_layers')
          if (lroot) then
            print*,'initialize_entropy: set boundary temperatures for spherical shell problem'
          endif
!
!  calulate temperature gradient from polytropic index
!
          call get_cp1(cp1)
          beta1=cp1*g0/(mpoly+1)*gamma/gamma1
!
!  temperatures at shell boundaries
!
          TT_ext=T0
          if (initss(1) .eq. 'shell_layers') then
!           lmultilayer=.true.  ! this is the default...
            if (hcond1==impossible) hcond1=(mpoly1+1.)/(mpoly0+1.)
            beta0=-g0/(mpoly0+1)*gamma/gamma1
            beta1=-g0/(mpoly1+1)*gamma/gamma1
            TT_crit=TT_ext+beta0*(r_bcz-r_ext)
            TT_int=TT_crit+beta1*(r_int-r_bcz)
          else
            lmultilayer=.false.  ! to ensure that hcond=cte
            TT_int=TT_ext*(1.+beta1*(r_ext/r_int-1.))
          endif
          if (lroot) then
            print*,'initialize_entropy: g0,mpoly,beta1',g0,mpoly,beta1
            print*,'initialize_entropy: TT_int, TT_ext=',TT_int,TT_ext
          endif
!         set up cooling parameters for spherical shell in terms of
!         sound speeds
          call get_soundspeed(log(TT_ext),cs2_ext)
          call get_soundspeed(log(TT_int),cs2_int)
          cs2cool=cs2_ext
!
      endselect
!
!  For global density gradient beta=H/r*dlnrho/dlnr, calculate actual
!  gradient dlnrho/dr = beta/H
!
      if (maxval(abs(beta_glnrho_global))/=0.0) then
        beta_glnrho_scaled=beta_glnrho_global*Omega/cs0
        if (lroot) print*, 'initialize_entropy: Global density gradient '// &
            'with beta_glnrho_global=', beta_glnrho_global
      endif
!
!  Turn off pressure gradient term and advection for 0-D runs.
!
      if (nxgrid*nygrid*nzgrid==1) then
        lpressuregradient_gas=.false.
        ladvection_entropy=.false.
        print*, 'initialize_entropy: 0-D run, turned off pressure gradient term'
        print*, 'initialize_entropy: 0-D run, turned off advection of entropy'
      endif
!
!  Initialize heat conduction.
!
      lheatc_Kconst=.false.
      lheatc_simple=.false.
      lheatc_chiconst=.false.
      lheatc_tensordiffusion=.false.
      lheatc_spitzer=.false.
      lheatc_corona=.false.
      lheatc_shock=.false.
      lheatc_hyper3ss=.false.
!
      lnothing=.false.
!
!  select which radiative heating we are using
!
      if (lroot) print*,'initialize_entropy: nheatc_max,iheatcond=',nheatc_max,iheatcond(1:nheatc_max)
      do i=1,nheatc_max
        select case (iheatcond(i))
        case('K-profile', 'K-const')
          lheatc_Kconst=.true.
          if (lroot) print*, 'heat conduction: K-profile'
        case('simple')
          lheatc_simple=.true.
          if (lroot) print*, 'heat conduction: simple'
        case('chi-const')
          lheatc_chiconst=.true.
          if (lroot) print*, 'heat conduction: constant chi'
        case ('tensor-diffusion')      
          lheatc_tensordiffusion=.true.
          if (lroot) print*, 'heat conduction: tensor diffusion'
        case ('spitzer')
          lheatc_spitzer=.true.
          if (lroot) print*, 'heat conduction: spitzer'
        case ('corona')
          lheatc_corona=.true.
          if (lroot) print*, 'heat conduction: corona'
        case ('shock')
          lheatc_shock=.true.
          if (lroot) print*, 'heat conduction: shock'
        case ('hyper3_ss')
          lheatc_hyper3ss=.true.
          if (lroot) print*, 'heat conduction: hyperdiffusivity of ss'
        case ('nothing')
          if (lroot .and. (.not. lnothing)) print*,'heat conduction: nothing'
        case default
          if (lroot) then
            write(unit=errormsg,fmt=*)  &
                'No such value iheatcond = ', trim(iheatcond(i))
            call fatal_error('initialize_entropy',errormsg)
          endif
        endselect
        lnothing=.true.
      enddo

      if (lhcond_global) then
        call farray_register_global("hcond",iglobal_hcond)
        call farray_register_global("glhc",iglobal_glhc,vector=3)
        do n=n1,n2
        do m=m1,m2
          p%r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
          p%z_mn=spread(z(n),1,nx)
          call heatcond(p,hcond)
          call gradloghcond(p,glhc)
          f(l1:l2,m,n,iglobal_hcond)=hcond
          f(l1:l2,m,n,iglobal_glhc:iglobal_glhc+2)=glhc
        enddo
        enddo
      endif
!
!  A word of warning...
!
      if (lheatc_Kconst .and. hcond0==0.0) then
        call warning('initialize_entropy', 'hcond0 is zero!')
      endif
      if (lheatc_chiconst .and. chi==0.0) then
        call warning('initialize_entropy','chi is zero!')
      endif
      if (all(iheatcond=='nothing') .and. hcond0/=0.0) then
        call warning('initialize_entropy', 'No heat conduction, but hcond0 /= 0')
      endif
      if (lheatc_simple .and. Kbot==0.0) then
        call warning('initialize_entropy','Kbot is zero!')
      endif
      if ((lheatc_spitzer.or.lheatc_corona) .and. (Kgpara==0.0 .or. Kgperp==0.0) ) then
        call warning('initialize_entropy','Kgperp or Kgpara is zero!')
      endif
      if (lheatc_hyper3ss .and. chi_hyper3==0.0) then
        call warning('initialize_entropy','chi_hyper3 is zero!')
      endif
      if (lheatc_shock .and. chi_shock==0.0) then
        call warning('initialize_entropy','chi_shock is zero!')
      endif
!
      if (NO_WARN) print*,f,lstarting  !(to keep compiler quiet)
!        
      endsubroutine initialize_entropy
!***********************************************************************
    subroutine read_entropy_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=entropy_init_pars,ERR=99, IOSTAT=iostat)
      else 
        read(unit,NML=entropy_init_pars,ERR=99) 
      endif

99    return
    endsubroutine read_entropy_init_pars
!***********************************************************************
    subroutine write_entropy_init_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=entropy_init_pars)
    endsubroutine write_entropy_init_pars
!***********************************************************************
    subroutine read_entropy_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=entropy_run_pars,ERR=99, IOSTAT=iostat)
      else 
        read(unit,NML=entropy_run_pars,ERR=99) 
      endif

99    return
    endsubroutine read_entropy_run_pars
!***********************************************************************
    subroutine write_entropy_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=entropy_run_pars)
    endsubroutine write_entropy_run_pars
!!***********************************************************************
    subroutine init_ss(f,xx,yy,zz)
!
!  initialise entropy; called from start.f90
!  07-nov-2001/wolf: coded
!  24-nov-2002/tony: renamed for consistancy (i.e. init_[variable name]) 
!
      use Cdata
      use Sub
      use Gravity
      use General, only: chn
      use Initcond
      use EquationOfState,  only: mpoly, isothtop, &
                                mpoly0, mpoly1, mpoly2, cs2cool, cs0, &
                                rho0, lnrho0, isothermal_entropy, &
                                isothermal_lnrho_ss, eoscalc, ilnrho_pp, &
                                eosperturb
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz,tmp,pot
      real, dimension (nx) :: pp,lnrho,ss
      real, dimension (mx) :: ss_mx 
      real :: cs2int,ssint,ztop,ss_ext,pot0,pot_ext
      logical :: lnothing=.true., save_pretend_lnTT
!
      intent(in) :: xx,yy,zz
      intent(inout) :: f
!
!mee Rather initialise everything as if entropy then convert it all at the end!
!mee      if (pretend_lnTT) f(:,:,:,iss)=f(:,:,:,iss)+(f(:,:,:,ilnrho)*gamma1-log(gamma1))/gamma

!
! If pretend_lnTT is set then turn it off so that initial conditions are
! correctly generated in terms of entropy, but then restore it later
! when we convert the data back again.
!
      save_pretend_lnTT=pretend_lnTT
      pretend_lnTT=.false.

      do iinit=1,ninit
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
        case('blob'); call blob(ampl_ss,f,iss,radius_ss,center1_x,center1_y,center1_z)
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

        case('mixinglength')
           call mixinglength(cs2cool,mixinglength_flux,f)
           if (ampl_ss/=0.) then
            print*,'init_ss: put bubble: radius_ss,ampl_ss=',radius_ss,ampl_ss
            tmp=xx**2+yy**2+(zz-0.5*(z1+z2))**2
            f(:,:,:,iss)=f(:,:,:,iss)+ampl_ss*exp(-tmp/max(radius_ss**2-tmp,1e-20))
            hcond0=-gamma/(gamma-1)*gravz/(mpoly0+1)/mixinglength_flux
            print*,'init_ss: hcond0=',hcond0
           endif

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
          if (pp_const==0.) then
            f(:,:,:,iss) = -(f(:,:,:,ilnrho)-lnrho0)
          else
            pp=pp_const
            do n=n1,n2
            do m=m1,m2
              lnrho=f(l1:l2,m,n,ilnrho)
              call eoscalc(ilnrho_pp,lnrho,pp,ss=ss)
              f(l1:l2,m,n,iss)=ss
            enddo
            enddo
          endif

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
          ! BEWARNED: isentropic star requires initlnrho=initss=isentropic-star
               call fatal_error('isentropic-star','requires density.f90')
          if (lgravr) then
            if (lroot) print*, &
                 'init_lnrho: isentropic star with isothermal atmosphere'
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
                   call fatal_error('init_ss',"inconsistency - cs2cool can't be 0")
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

        case ('shell_layers')
          !
          ! radial temperature profiles for spherical shell problem
          !
          call information('init_ss',' two polytropic layers in a spherical shell')
          call shell_ss_layers(f)

        case ('polytropic_simple')     
          !
          ! vertical temperature profiles for convective layer problem
          !
          call layer_ss(f)

        case default
          !
          !  Catch unknown values
          !
          write(unit=errormsg,fmt=*) 'No such value for initss(' &
                           //trim(iinit_str)//'): ',trim(initss(iinit))
          call fatal_error('init_ss',errormsg)

      endselect

      if (lroot) print*,'init_ss: initss(' &
                        //trim(iinit_str)//') = ',trim(initss(iinit))

      endif

      enddo

      if (lnothing.and.lroot) print*,'init_ss: nothing'
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
        write (unit=errormsg,fmt=*) 'No such value for pertss:', pertss
        call fatal_error('init_ss',errormsg)

      endselect
!
!  replace ss by lnTT when pretend_lnTT is true
!
      if (save_pretend_lnTT) then
        pretend_lnTT=.true.
        do m=1,my
        do n=1,mz
          ss_mx=f(:,m,n,iss)
          call eosperturb(f,mx,ss=ss_mx)
        enddo; enddo
      endif
!
      if (NO_WARN) print*,xx,yy  !(to keep compiler quiet)
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
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: tmp,p,zz
      real, dimension (mz) :: stp
      real :: mpoly,zint,zbot,zblend,beta1,cs2int,ssint
      integer :: isoth
!
!  Warning: beta1 is here not dT/dz, but dcs2/dz = (gamma-1)c_p dT/dz
!
      if (isoth /= 0) then ! isothermal layer
        beta1 = 0.
        tmp = ssint - gamma1*gravz*(zz-zint)/cs2int
        ssint = -gamma1*gravz*(zbot-zint)/cs2int ! ss at layer interface
        if (headt) print*,'ssint=',ssint
      else
        beta1 = gamma*gravz/(mpoly+1)
        tmp = 1 + beta1*(zz-zint)/cs2int
        ! Abort if args of log() are negative
        if (any((tmp <= 0.) .and. (zz <= zblend))) then
          call fatal_error('polytropic_ss_z', &
              'Imaginary entropy values -- your z_inf is too low.')
        endif
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
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: tmp,p,zz
      real, dimension (mz) :: stp
      real :: mpoly,zint,zbot,zblend,beta1,cs2int,ssint, nu_epicycle2
      integer :: isoth
!
!  Warning: beta1 is here not dT/dz, but dcs2/dz = (gamma-1)c_p dT/dz
!
      nu_epicycle2 = nu_epicycle**2
      if (isoth /= 0) then ! isothermal layer
        beta1 = 0.
        tmp = ssint - gamma1*gravz*nu_epicycle2*(zz**2-zint**2)/cs2int/2.
        ssint = -gamma1*gravz*nu_epicycle2*(zbot**2-zint**2)/cs2int/2. 
              ! ss at layer interface
      else
        beta1 = gamma*gravz*nu_epicycle2/(mpoly+1)
        tmp = 1 + beta1*(zz**2-zint**2)/cs2int/2.
        ! Abort if args of log() are negative
        if (any((tmp <= 0.) .and. (zz <= zblend))) then
          call fatal_error('polytropic_ss_disc', &
              'Imaginary entropy values -- your z_inf is too low.')
        endif
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
      use EquationOfState, only: eoscalc,ilnrho_ss
      use Gravity, only: gravz

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, intent(in) :: lnrho_bot,ss_const
      real :: cs2_,lnrho,lnrho_m
!
      if (.not. lgravz) then
        call fatal_error("hydrostatic_isentropic","Currently only works for vertical gravity field")
      endif

      !
      ! In case this processor is not located at the very bottom
      ! perform integration through lower lying processors
      !
      lnrho=lnrho_bot
      do n=1,nz*ipz
        call eoscalc(ilnrho_ss,lnrho,ss_const,cs2=cs2_)
        lnrho_m=lnrho+dz*gravz/cs2_/2
        call eoscalc(ilnrho_ss,lnrho_m,ss_const,cs2=cs2_)
        lnrho=lnrho+dz*gravz/cs2_
      enddo

      !
      ! Do the integration on this processor
      !
      f(:,:,n1,ilnrho)=lnrho
      do n=n1+1,n2
        call eoscalc(ilnrho_ss,lnrho,ss_const,cs2=cs2_)
        lnrho_m=lnrho+dz*gravz/cs2_/2
        call eoscalc(ilnrho_ss,lnrho_m,ss_const,cs2=cs2_)
        lnrho=lnrho+dz*gravz/cs2_
        f(:,:,n,ilnrho)=lnrho
      enddo

      !
      ! Entropy is simply constant
      !
      f(:,:,:,iss)=ss_const

    endsubroutine hydrostatic_isentropic
!***********************************************************************
    subroutine mixinglength(cs2cool,mixinglength_flux,f)
!
!  Mixing length initial condition.
!
!  ds/dz=-HT1*(F/rho*cs3)^(2/3) in the convection zone.
!  ds/dz=-HT1*(1-F/gK) in the radiative interior, where ...
!  Solves dlnrho/dz=-ds/dz-gravz/cs2 using 2nd order Runge-Kutta.
!
!  Currently only works for vertical gravity field.
!  Starts at bottom boundary where the density has to be set in the gravity
!  module. Use mixinglength_flux as flux in convection zone (no further
!  scaling is applied, ie no further free parameter is assumed.)
!
!  12-jul-05/axel: coded
!  17-Nov-2005/dintrans: updated using strat_MLT
!
      use Cdata
      use Gravity, only: gravz, z1
      use General, only: safe_character_assign
!--   use EquationOfState, only: mpoly1
      use EquationOfState, only: gamma, gamma1, cs2top
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nzgrid) :: cs2m,lnrhom,ssm
      real :: zm,ztop,cs2cool,mixinglength_flux
      real :: zbot,rbot,rt_old,rt_new,rb_old,rb_new,crit, &
              rhotop,rhobot
      integer :: iter
      character (len=120) :: wfile
!
      if (headtt) print*,'init_ss : mixinglength stratification'
      if (.not.lgravz) then
        call fatal_error("mixinglength","works only for vertical gravity")
      endif
!
!  do the calculation on all processors, and then put the relevant piece
!  into the f array.
!  choose value zbot where rhobot should be applied and give two first
!  estimates for rhotop
!
      zbot=z1
      rbot=1.
      rt_old=.1*rbot
      rt_new=.12*rbot
!
!  need to iterate for rhobot=1.
!  produce first estimate
!
      rhotop=rt_old
      cs2top=cs2cool
      call strat_MLT (rhotop, cs2top, mixinglength_flux, lnrhom, &
                  ssm, cs2m, rhobot)
      rb_old=rhobot
!
!  next estimate
!
      rhotop=rt_new
      call strat_MLT (rhotop, cs2top, mixinglength_flux, lnrhom, &
                  ssm, cs2m, rhobot)
      rb_new=rhobot

      do 10 iter=1,10
!
!  new estimate
!
        rhotop=rt_old+(rt_new-rt_old)/(rb_new-rb_old)*(rbot-rb_old)
!
!  check convergence
!
        crit=abs(rhotop/rt_new-1.)
        if (crit.le.1e-4) goto 20
!
        call strat_MLT (rhotop, cs2top, mixinglength_flux, lnrhom, &
                    ssm, cs2m, rhobot)
!
!  update new estimates
!
        rt_old=rt_new
        rb_old=rb_new
        rt_new=rhotop
        rb_new=rhobot
   10 continue
   20 if (ipz.eq.0) print*,'- iteration completed: rhotop,crit=',rhotop,crit
!
! put density and entropy into f-array
! write the initial stratification in data/proc*/stratMLT.dat
!
      ztop=xyz0(3)+Lxyz(3)
      call safe_character_assign(wfile,trim(directory)//'/stratMLT.dat')
      open(11+ipz,file=wfile,status='unknown')
      do n=1,nz
        iz=n+ipz*nz
        zm=ztop-(iz-1)*dz
        f(:,:,n+nghost,ilnrho)=lnrhom(nzgrid-iz+1)
        f(:,:,n+nghost,iss)=ssm(nzgrid-iz+1)
        write(11+ipz,'(4(2x,1pe12.5))') zm,exp(lnrhom(iz)),ssm(iz),cs2m(iz) 
      enddo
      close(11+ipz)
      return
!
    endsubroutine mixinglength
!***********************************************************************
    subroutine shell_ss(f)
!
!  Initialize entropy based on specified radial temperature profile in
!  a spherical shell
!
!  20-oct-03/dave -- coded
!
      use Gravity, only: g0
      use EquationOfState, only: eoscalc, ilnrho_lnTT, mpoly, get_cp1

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: lnrho,lnTT,TT,ss,pert_TT,r_mn
      real :: beta1,cp1,lnrho_dummy=0.
!
!  beta1 is the temperature gradient
!  1/beta = (g/cp) 1./[(1-1/gamma)*(m+1)]
!
      call get_cp1(cp1)
      beta1=cp1*g0/(mpoly+1)*gamma/gamma1
!
!  set intial condition
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
!
        r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
        call shell_ss_perturb(pert_TT)
!
        where (r_mn >= r_ext) TT = TT_ext
        where (r_mn < r_ext .AND. r_mn > r_int) &
          TT = TT_ext*(1.+beta1*(r_ext/r_mn-1))+pert_TT
        where (r_mn <= r_int) TT = TT_int
!
        lnrho=f(l1:l2,m,n,ilnrho)
        lnTT=log(TT)
        call eoscalc(ilnrho_lnTT,lnrho,lnTT,ss=ss)
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
      real, dimension (nx), intent(out) :: pert_TT
      real, dimension (nx) :: xr,cos_4phi,sin_theta4,r_mn,rcyl_mn,phi_mn
      real :: ampl0=.885065
!
      select case(initss(1))
!
        case ('geo-kws')
          pert_TT=0.
!
        case ('geo-benchmark')
          r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
          rcyl_mn=sqrt(x(l1:l2)**2+y(m)**2)
          phi_mn=atan2(spread(y(m),1,nx),x(l1:l2))
          xr=2*r_mn-r_int-r_ext              ! radial part of perturbation
          cos_4phi=cos(4*phi_mn)             ! azimuthal part
          sin_theta4=(rcyl_mn/r_mn)**4       ! meridional part
          pert_TT=ampl0*ampl_TT*(1-3*xr**2+3*xr**4-xr**6)*sin_theta4*cos_4phi
!
      endselect
!
    endsubroutine shell_ss_perturb
!***********************************************************************
    subroutine layer_ss(f)
!
!  initial state entropy profile used for initss='polytropic_simple',
!  for `conv_slab' style runs, with a layer of polytropic gas in [z0,z1].
!  generalised for cp/=1.
!
      use Cdata
      use Gravity, only: gravz
      use EquationOfState, only: eoscalc, ilnrho_lnTT, mpoly, get_cp1
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f

      real, dimension (nx) :: lnrho,lnTT,TT,ss,z_mn
      real :: beta1,cp1,lnrho_dummy=0.
!
!  beta1 is the temperature gradient
!  1/beta = (g/cp) 1./[(1-1/gamma)*(m+1)]
!
      call get_cp1(cp1)
      beta1=cp1*gamma/gamma1*gravz/(mpoly+1)
!
!  set initial condition (first in terms of TT, and then in terms of ss)
!
      do m=m1,m2
      do n=n1,n2
        z_mn = spread(z(n),1,nx)
        TT = beta1*z_mn
        lnrho=f(l1:l2,m,n,ilnrho)
        lnTT=log(TT)
        call eoscalc(ilnrho_lnTT,lnrho,lnTT,ss=ss)
        f(l1:l2,m,n,iss)=ss
      enddo
      enddo
!
    endsubroutine layer_ss
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
      use Mpicomm, only: mpibcast_real
      use EquationOfState, only: eoscalc, ilnrho_pp, getmu,rho0, eosperturb
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(nx) :: absz
      double precision, dimension(nx) :: n_c,n_w,n_i,n_h
!  T in K, k_B s.t. pp is in code units ( = 9.59e-15 erg/cm/s^2)
!  (i.e. k_B = 1.381e-16 (erg/K) / 9.59e-15 (erg/cm/s^2) )
      real, parameter :: T_c_cgs=500.0,T_w_cgs=8.0e3,T_i_cgs=8.0e3,T_h_cgs=1.0e6 
      real :: T_c,T_w,T_i,T_h
      real, dimension(nx) :: rho,pp,lnrho,ss
!      real, dimension(nx) :: pp 
!     double precision :: pp0 
!      real, dimension(2) :: fmpi2
      real, dimension(1) :: fmpi1
      real :: kpc
      double precision ::  rhoscale
!      integer :: iproctop
!
      if (lroot) print*,'ferriere: Ferriere density and entropy profile'
!
!  first define reference values of pp, cs2, at midplane.  
!  pressure is set to 6 times thermal pressure, this factor roughly
!  allowing for other sources, as modelled by Ferriere.
!
  !    call getmu(mu)
      kpc = 3.086D21 / unit_length
      rhoscale = 1.36 * m_p * unit_length**3
      print *,'ferrier: kpc, rhoscale =',kpc,rhoscale !,mu
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
        n_c=0.340*(0.859*exp(-(z(n)/(0.127*kpc))**2) +         &
                   0.047*exp(-(z(n)/(0.318*kpc))**2) +         &
                   0.094*exp(-absz/(0.403*kpc)))     
!  warm gas profile n_w (eq 7)
        n_w=0.226*(0.456*exp(-(z(n)/(0.127*kpc))**2) +  &
                   0.403*exp(-(z(n)/(0.318*kpc))**2) +  &
                   0.141*exp(-absz/(0.403*kpc)))
!  ionized gas profile n_i (eq 9)
        n_i=0.0237*exp(-absz/kpc) + 0.0013* exp(-absz/(0.150*kpc))
!  hot gas profile n_h (eq 13)
        n_h=0.00048*exp(-absz/(1.5*kpc))         
!  normalised s.t. rho0 gives mid-plane density directly (in 10^-24 g/cm^3)
        !rho=rho0/(0.340+0.226+0.025+0.00048)*(n_c+n_w+n_i+n_h)*rhoscale
        rho=real((n_c+n_w+n_i+n_h)*rhoscale)
        lnrho=log(rho)
        f(l1:l2,m,n,ilnrho)=lnrho

!  define entropy via pressure, assuming fixed T for each component
        if (lentropy) then
!  thermal pressure (eq 15)
          pp=real(k_B*unit_length**3 * &
             (1.09*n_c*T_c + 1.09*n_w*T_w + 2.09*n_i*T_i + 2.27*n_h*T_h))
!           
          call eosperturb(f,nx,pp=pp) 
          ss=f(l1:l2,m,n,ilnrho)
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
!     -----------------------------------------------------
!       subroutine ic2rp2(maxmx,maxmy,meqn,mbc,mx,my,x,y,dx,dy,q)
!     -----------------------------------------------------
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
      real, dimension (mx,my,mz,mfarray) :: f
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
    if (NO_WARN) print*,zz
    endsubroutine shock2d
!***********************************************************************
    subroutine pencil_criteria_entropy()
! 
!  All pencils that the Entropy module depends on are specified here.
! 
!  20-11-04/anders: coded
!
      use Cdata
      use EquationOfState, only: beta_glnrho_scaled
!
      lpenc_requested(i_cp1)=.true.
      if (ldt) lpenc_requested(i_cs2)=.true.
      if (lpressuregradient_gas) lpenc_requested(i_fpres)=.true.
      if (ladvection_entropy) lpenc_requested(i_ugss)=.true.
      if (lviscosity_heat) lpenc_requested(i_visc_heat)=.true.
      if (tau_cor>0.0) then
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_rho1)=.true.
      endif
      if (tauheat_buffer>0.0) then
        lpenc_requested(i_ss)=.true.
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_rho1)=.true.
      endif
      if (cool/=0.0 .or. cool_ext/=0.0 .or. cool_int/=0.0) &
          lpenc_requested(i_cs2)=.true.
      if (lgravz .and. (luminosity/=0.0 .or. cool/=0.0)) &
          lpenc_requested(i_cs2)=.true.
      if (luminosity/=0 .or. cool/=0 .or. tau_cor/=0 .or. &
          tauheat_buffer/=0 .or. heat_uniform/=0 .or. &
          (cool_ext/=0 .and. cool_int /= 0) .or. lturbulent_heat) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_TT1)=.true.
      endif
      if (pretend_lnTT) lpenc_requested(i_divu)=.true.
      if (lgravr) lpenc_requested(i_r_mn)=.true.
      if (lheatc_simple) then 
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
      endif
      if (lheatc_Kconst) then 
        if (hcond0/=0) then
          lpenc_requested(i_rho1)=.true.
          lpenc_requested(i_glnTT)=.true.
          lpenc_requested(i_del2lnTT)=.true.
          if (lmultilayer) then
            if (lgravz) then
              lpenc_requested(i_z_mn)=.true.
            else
              lpenc_requested(i_r_mn)=.true.
            endif
          endif
        endif
        if (chi_t/=0) then
          lpenc_requested(i_del2ss)=.true.
        endif
      endif
      if (lheatc_spitzer) then 
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_bij)=.true.
      endif
      if (lheatc_corona) then
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true. 
        lpenc_requested(i_bb)=.true.       
        lpenc_requested(i_bij)=.true.
      endif
      if (lheatc_chiconst) then 
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
      if (lheatc_tensordiffusion) then 
        lpenc_requested(i_bij)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true. 
      endif
      if (lheatc_shock) then 
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_gss)=.true.
        lpenc_requested(i_del2ss)=.true.
        lpenc_requested(i_gshock)=.true.
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_glnTT)=.true.
      endif
      if (lheatc_hyper3ss) lpenc_requested(i_del6ss)=.true.
      if (cooltype=='shell') then
        lpenc_requested(i_r_mn)=.true.
        if (deltaT_poleq/=0.) then
          lpenc_requested(i_z_mn)=.true.
          lpenc_requested(i_rcyl_mn)=.true.
        endif
      endif
!
      if (maxval(abs(beta_glnrho_scaled))/=0.0) lpenc_requested(i_cs2)=.true.
!
      lpenc_diagnos2d(i_ss)=.true.
!
      if (idiag_dtchi/=0) lpenc_diagnos(i_rho1)=.true.
      if (idiag_ethdivum/=0) lpenc_diagnos(i_divu)=.true.
      if (idiag_csm/=0) lpenc_diagnos(i_cs2)=.true.
      if (idiag_eem/=0) lpenc_diagnos(i_ee)=.true.
      if (idiag_ppm/=0) lpenc_diagnos(i_pp)=.true.
      if (idiag_ssm/=0 .or. idiag_ssmz/=0 .or. idiag_ssmy/=0.or.idiag_ssmx/=0) &
           lpenc_diagnos(i_ss)=.true.
      lpenc_diagnos(i_rho)=.true.
      lpenc_diagnos(i_ee)=.true.
      if (idiag_eth/=0 .or. idiag_ethtot/=0 .or. idiag_ethdivum/=0 ) then
          lpenc_diagnos(i_rho)=.true.
          lpenc_diagnos(i_ee)=.true.
      endif
      if (idiag_fconvz/=0 .or. idiag_fturbz/=0 ) then
          lpenc_diagnos(i_rho)=.true.
          lpenc_diagnos(i_TT)=.true.  !(to be replaced by enthalpy)
      endif
      if (idiag_TTm/=0 .or. idiag_TTmz/=0 .or. idiag_TTmax/=0 &
        .or. idiag_TTmin/=0) &
          lpenc_diagnos(i_TT)=.true.
      if (idiag_yHm/=0 .or. idiag_yHmax/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_dtc/=0) lpenc_diagnos(i_cs2)=.true.
      if (idiag_TTp/=0) then 
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_cs2)=.true.
        lpenc_diagnos(i_rcyl_mn)=.true.
      endif
!
    endsubroutine pencil_criteria_entropy
!***********************************************************************
    subroutine pencil_interdep_entropy(lpencil_in)
!       
!  Interdependency among pencils from the Entropy module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_ugss)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_gss)=.true.
      endif
      if (lpencil_in(i_Ma2)) then
        lpencil_in(i_u2)=.true.
        lpencil_in(i_cs2)=.true.
      endif
      if (lpencil_in(i_fpres)) then
        lpencil_in(i_cs2)=.true.
        lpencil_in(i_glnrho)=.true.
        lpencil_in(i_gss)=.true.
        if (leos_idealgas) then
          lpencil_in(i_glnTT)=.true.
        else
          lpencil_in(i_cp1)=.true.
        endif
      endif
!
    endsubroutine pencil_interdep_entropy
!***********************************************************************
    subroutine calc_pencils_entropy(f,p)
!       
!  Calculate Entropy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      use EquationOfState
      use Sub
!      
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      integer :: j
!
      intent(in) :: f
      intent(inout) :: p
! Ma2
      if (lpencil(i_Ma2)) p%Ma2=p%u2/p%cs2
! ugss
      if (lpencil(i_ugss)) &
          call u_dot_grad(f,iss,p%gss,p%uu,p%ugss,UPWIND=lupw_ss)
! fpres
      if (lpencil(i_fpres)) then
        if (leos_idealgas) then
          do j=1,3
            p%fpres(:,j)=-p%cs2*(p%glnrho(:,j) + p%glnTT(:,j))*gamma11
          enddo
!  TH: The following would work if one uncomments the intrinsic operator
!  extensions in sub.f90. Please Test.
!          p%fpres      =-p%cs2*(p%glnrho + p%glnTT)*gamma11
        else
          do j=1,3
            p%fpres(:,j)=-p%cs2*(p%glnrho(:,j) + p%cp1tilde*p%gss(:,j))
          enddo
        endif
      endif
!
    endsubroutine calc_pencils_entropy
!**********************************************************************
    subroutine dss_dt(f,df,p)
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
      use EquationOfState, only: beta_glnrho_global, beta_glnrho_scaled, gamma11, cs0
      use Sub
      use Global
      use Special, only: special_calc_entropy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: rhs,Hmax=0.
      real, dimension (nx) :: vKpara,vKperp
      real :: zbot,ztop,xi,profile_cor
      real :: uT
      integer :: j,ju
!
      intent(inout)  :: f,p
      intent(out) :: df
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
!  calculate cs2, TT1, and cp1 in a separate routine
!  With IONIZATION=noionization, assume perfect gas with const coeffs
!
      if (headtt) print*,'dss_dt: lnTT,cs2,cp1=', &
          p%lnTT(1), p%cs2(1), p%cp1(1)
!
!  ``cs2/dx^2'' for timestep
!
      if (lfirst.and.ldt) advec_cs2=p%cs2*dxyz_2
      if (headtt.or.ldebug) print*,'dss_dt: max(advec_cs2) =',maxval(advec_cs2)
!
      if (lhydro) then
!
!  Pressure term in momentum equation (setting lpressuregradient_gas to
!  .false. allows suppressing pressure term for test purposes)
!
        if (lpressuregradient_gas) then
          do j=1,3
            ju=j+iuu-1
            df(l1:l2,m,n,ju) = df(l1:l2,m,n,ju) + p%fpres(:,j)
          enddo
        endif
!
!  velocity damping in the coronal heating zone
!
        if (tau_cor>0) then
          if (z(n)>=z_cor) then
            xi=(z(n)-z_cor)/(ztop-z_cor)
            profile_cor=xi**2*(3-2*xi)
            df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - &
                profile_cor*f(l1:l2,m,n,iux:iuz)/tau_cor
          endif
        endif
!
      endif
!
!  advection term
!
      if (ladvection_entropy) df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) - p%ugss
!
!  if pretend_lnTT=.true., we pretend that ss is actually lnTT
!
      if (pretend_lnTT) df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)-p%divu*p%ss
!
!  Add pressure force from global density gradient.
!  
      if (maxval(abs(beta_glnrho_global))/=0.0) then
        if (headtt) print*, 'dss_dt: adding global pressure gradient force'
        do j=1,3
          df(l1:l2,m,n,(iux-1)+j) = df(l1:l2,m,n,(iux-1)+j) &
              - 1/gamma*p%cs2*beta_glnrho_scaled(j)
        enddo
      endif
!
!  Calculate viscous contribution to entropy
!
      if (lviscosity) then
         if (lviscosity_heat)   call calc_viscous_heat(df,p,Hmax)
      endif
!
!  thermal conduction
!
      if (lheatc_Kconst) call calc_heatcond(f,df,p)
      if (lheatc_simple) call calc_heatcond_simple(df,p)
      if (lheatc_chiconst) call calc_heatcond_constchi(df,p)
      if (lheatc_tensordiffusion) then
        vKpara(:) = Kgpara
        vKperp(:) = Kgperp
        call tensor_diffusion_coef(p%glnTT,p%hlnTT,p%bij,p%bb,vKperp,vKpara,rhs,llog=.true.)
        df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)+rhs*p%rho1
        if (lfirst.and.ldt) then         
           dt1_max=max(dt1_max,maxval(abs(rhs*p%rho1)*gamma)/(cdts))
        endif
      endif
      if (lheatc_spitzer) call calc_heatcond_spitzer(df,p)
      if (lheatc_corona) then
        call calc_heatcond_spitzer(df,p)
        call newton_cool(df,p)
        call calc_heat_cool_RTV(df,p)     
      endif
      if (lheatc_shock) call calc_heatcond_shock(df,p)
      if (lheatc_hyper3ss) call calc_heatcond_hyper3(df,p)
!
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
          call calc_heat_cool(df,p,Hmax)
!
!  interstellar radiative cooling and UV heating
!
      if (linterstellar) &
          call calc_heat_cool_interstellar(f,df,p,Hmax)
!
!  possibility of entropy relaxation in exterior region
!
      if (tau_ss_exterior/=0.) call calc_tau_ss_exterior(df,p)

!  entry possibility for "personal" entries.
!  In that case you'd need to provide your own "special" routine.
!
      if (lspecial) call special_calc_entropy(f,df,p)
!
!  Apply border profile
!
      if (lborder_profiles) call set_border_entropy(f,df,p)
!
!  phi-averages
!  Note that this does not necessarily happen with ldiagnos=.true.
!
      if (l2davgfirst) then
        call phisum_mn_name_rz(p%ss,idiag_ssmphi)
      endif
!
!  Enforce maximum heating rate timestep constraint
!
!      if (lfirst.and.ldt) dt1_max=max(dt1_max,Hmax/ee/cdts) 
!
!  Calculate entropy related diagnostics
!
      if (ldiagnos) then
        !uT=unit_temperature !(define shorthand to avoid long lines below)
        uT=1. !(AB: for the time being; to keep compatible with auto-test
        if (idiag_TTmax/=0) call max_mn_name(p%TT*uT,idiag_TTmax)
        if (idiag_TTmin/=0) call max_mn_name(-p%TT*uT,idiag_TTmin,lneg=.true.)
        if (idiag_TTm/=0) call sum_mn_name(p%TT*uT,idiag_TTm)
        if (idiag_yHmax/=0) call max_mn_name(p%yH,idiag_yHmax)
        if (idiag_yHm/=0) call sum_mn_name(p%yH,idiag_yHm)
        if (idiag_dtc/=0) &
            call max_mn_name(sqrt(advec_cs2)/cdt,idiag_dtc,l_dt=.true.)
        if (idiag_eth/=0) call sum_mn_name(p%rho*p%ee,idiag_eth)
        if (idiag_ethtot/=0) call integrate_mn_name(p%rho*p%ee,idiag_ethtot)
        if (idiag_ethdivum/=0) &
            call sum_mn_name(p%rho*p%ee*p%divu,idiag_ethdivum)
        if (idiag_ssm/=0) call sum_mn_name(p%ss,idiag_ssm)
        if (idiag_eem/=0) call sum_mn_name(p%ee,idiag_eem)
        if (idiag_ppm/=0) call sum_mn_name(p%pp,idiag_ppm)
        if (idiag_csm/=0) call sum_mn_name(p%cs2,idiag_csm,lsqrt=.true.)
        if (idiag_ugradpm/=0) &
            call sum_mn_name(p%cs2*(p%uglnrho+p%ugss),idiag_ugradpm)
        if (idiag_TTp/=0) call sum_lim_mn_name(p%rho*p%cs2*gamma11,idiag_TTp,p)
!
!  xy averages for fluxes; doesn't need to be as frequent (check later)
!  idiag_fradz is done in the calc_headcond routine
!
        if (idiag_fconvz/=0) call xysum_mn_name_z(p%rho*p%uu(:,3)*p%TT,idiag_fconvz)
        if (idiag_ssmz/=0) call xysum_mn_name_z(p%ss,idiag_ssmz)
        if (idiag_ssmy/=0) call xzsum_mn_name_y(p%ss,idiag_ssmy)
        if (idiag_ssmx/=0) call yzsum_mn_name_x(p%ss,idiag_ssmx)
        if (idiag_TTmz/=0) call xysum_mn_name_z(p%TT,idiag_TTmz)
      endif
!
    endsubroutine dss_dt
!***********************************************************************
    subroutine set_border_entropy(f,df,p)
!
!  Calculates the driving term for the border profile
!  of the ss variable.
!
!  28-jul-06/wlad: coded
!
      use BorderProfiles, only: border_driving
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(nx) :: f_target,cs2_0
      real :: g0=1.,r0_pot=0.1
!
      select case(borderss)
!
      case('zero','0')
         f_target=0.
      case('constant')
         f_target=ss_const
      case('globaldisc')
         cs2_0=(cs20*p%rcyl_mn**2)*(g0*(p%rcyl_mn**2+r0_pot**2)**(-1.5))
         !           [  r      *                omega             ]**2     
         f_target= gamma11*log(cs2_0) !- gamma1*gamma11*lnrho
      case('nothing')
         if (lroot.and.ip<=5) &
              print*,"set_border_entropy: borderss='nothing'"
      case default
         write(unit=errormsg,fmt=*) &
              'set_border_entropy: No such value for borderss: ', &
              trim(borderss)
         call fatal_error('set_border_entropy',errormsg)
      endselect
!
      call border_driving(f,df,p,f_target,iss)
!
    endsubroutine set_border_entropy
!***********************************************************************
    subroutine calc_heatcond_constchi(df,p)
!
!  Heat conduction for constant value of chi=K/(rho*cp)
!  This routine also adds in turbulent diffusion, if chi_t /= 0.
!  Ds/Dt = ... + 1/(rho*T) grad(flux), where
!  flux = chi*rho*gradT + chi_t*rho*T*grads
!  This routine is currently not correct when ionization is used.
!
!  29-sep-02/axel: adapted from calc_heatcond_simple
!  12-mar-06/axel: used p%glnTT and p%del2lnTT, so that general cp work ok
!
      use Cdata
      use Sub
      use Gravity
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: thdiff,g2
!
      intent(out) :: df
!
!  check that chi is ok
!
      if (headtt) print*,'calc_heatcond_constchi: chi=',chi
!
!  Heat conduction
!  Note: these routines require revision when ionization turned on
!  The variable g2 is reused to calculate glnP.gss a few lines below.
!
!AB: keep the following lines for the record, but they aren't
!AB: correct for general values of cp.
!
!     glnT = gamma*p%gss + gamma1*p%glnrho
!     glnP = gamma*p%gss + gamma*p%glnrho
!     call dot(glnP,glnT,g2)
!     thdiff = chi * (gamma*p%del2ss+gamma1*p%del2lnrho + g2)
!     if (chi_t/=0.) then
!       call dot(glnP,p%gss,g2)
!       thdiff = thdiff + chi_t*(p%del2ss+g2)
!     endif
!
!  diffusion of the form:
!  rho*T*Ds/Dt = ... + nab.(rho*cp*chi*gradT)
!        Ds/Dt = ... + cp*chi*[del2lnTT+(glnrho+glnTT).glnTT]
!
!  with additional turbulent diffusion
!  rho*T*Ds/Dt = ... + nab.(rho*T*chit*grads)
!        Ds/Dt = ... + chit*[del2ss+(glnrho+glnTT).gss]
!
      if (pretend_lnTT) then
        call dot(p%glnrho+p%glnTT,p%glnTT,g2)
        thdiff=gamma*chi*(p%del2lnTT+g2)
        if (chi_t/=0.) then
          call dot(p%glnrho+p%glnTT,p%gss,g2)
          thdiff=thdiff+chi_t*(p%del2ss+g2)
        endif
      else
        call dot(p%glnrho+p%glnTT,p%glnTT,g2)
        !thdiff=cp*chi*(p%del2lnTT+g2)
!AB:  divide by p%cp1, since we don't have cp here.
        thdiff=chi*(p%del2lnTT+g2)/p%cp1
        if (chi_t/=0.) then
          call dot(p%glnrho+p%glnTT,p%gss,g2)
          thdiff=thdiff+chi_t*(p%del2ss+g2)
        endif
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
        if (leos_idealgas) then 
          diffus_chi=max(diffus_chi,(gamma*chi+chi_t)*dxyz_2)
        else
          diffus_chi=max(diffus_chi,(chi+chi_t)*dxyz_2)
        endif
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_constchi
!***********************************************************************
    subroutine calc_heatcond_hyper3(df,p)
!
!  Naive hyperdiffusivity of entropy.
!
!  17-jun-05/anders: coded
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!      
      real, dimension (nx) :: thdiff
!
      intent(out) :: df
!
!  check that chi_hyper3 is ok
!
      if (headtt) print*, 'calc_heatcond_hyper3: chi_hyper3=', chi_hyper3
!
!  Heat conduction
!
      thdiff = chi_hyper3 * p%del6ss
!
!  add heat conduction to entropy equation
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
      if (headtt) print*,'calc_heatcond_hyper3: added thdiff'
!
!  check maximum diffusion from thermal diffusion
!
      if (lfirst.and.ldt) then
        diffus_chi=max(diffus_chi,chi_hyper3*dxyz_4)
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_hyper3
!***********************************************************************
    subroutine calc_heatcond_shock(df,p)
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
      use Sub
      use Gravity
      use EquationOfState, only: gamma
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: thdiff,g2,gshockglnTT
!
      intent(in) :: p
      intent(out) :: df
!
!  check that chi is ok
!
      if(headtt) print*,'calc_heatcond_shock: chi_t,chi_shock=',chi_t,chi_shock
!
!  calculate terms for shock diffusion
!  Ds/Dt = ... + chi_shock*[del2ss + (glnchi_shock+glnpp).gss]
!
      call dot(p%gshock,p%glnTT,gshockglnTT)
      call dot(p%glnrho+p%glnTT,p%glnTT,g2)
!
!  shock entropy diffusivity
!  Write: chi_shock = chi_shock0*shock, and gshock=grad(shock), so
!  Ds/Dt = ... + chi_shock0*[shock*(del2ss+glnpp.gss) + gshock.gss]
!
      if (headtt) print*,'calc_heatcond_shock: use shock diffusion'
      if (pretend_lnTT) then
        thdiff=gamma*chi_shock*(p%shock*(p%del2lnrho+g2)+gshockglnTT)
        if (chi_t/=0.) then
          call dot(p%glnrho+p%glnTT,p%gss,g2)
          thdiff=thdiff+chi_t*(p%del2ss+g2)
        endif
      else
        thdiff=chi_shock*(p%shock*(p%del2lnrho+g2)+gshockglnTT)
        if (chi_t/=0.) then
          call dot(p%glnrho+p%glnTT,p%gss,g2)
          thdiff=thdiff+chi_t*(p%del2ss+g2)
        endif
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
      if (lfirst.and.ldt) then
        if (leos_idealgas) then
          diffus_chi=max(diffus_chi,(chi_t+gamma*chi_shock*p%shock)*dxyz_2)
        else
          diffus_chi=max(diffus_chi,(chi_t+chi_shock*p%shock)*dxyz_2)
        endif
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_shock
!***********************************************************************
    subroutine calc_heatcond_simple(df,p)
!
!  heat conduction
!
!   8-jul-02/axel: adapted from Wolfgang's more complex version
!  30-mar-06/ngrs: simplified calculations using p%glnTT and p%del2lnTT
!
      use Cdata
      use Sub
!
      type (pencil_case) :: p
      real, dimension (mx,my,mz,mvar) :: df
      !real, dimension (nx,3) :: glnT,glnThcond !,glhc
      real, dimension (nx) :: chix
      real, dimension (nx) :: thdiff,g2
      real, dimension (nx) :: hcond
!
      intent(in) :: p
      intent(out) :: df
!
!  This particular version assumes a simple polytrope, so mpoly is known
!
      hcond=Kbot
      if (headtt) then
        print*,'calc_heatcond_simple: hcond=', maxval(hcond)
      endif
!
!  Heat conduction
!  Note: these routines require revision when ionization turned on
!
      ! NB: the following left in for the record, but the version below,
      !     using del2lnTT & glnTT, is simpler
      !
      !chix = p%rho1*hcond*p%cp1                    ! chix = K/(cp rho)
      !glnT = gamma*p%gss*spread(p%cp1,2,3) + gamma1*p%glnrho ! grad ln(T)
      !glnThcond = glnT !... + glhc/spread(hcond,2,3)    ! grad ln(T*hcond)
      !call dot(glnT,glnThcond,g2)
      !thdiff =  p%rho1*hcond * (gamma*p%del2ss*p%cp1 + gamma1*p%del2lnrho + g2)  

      !  diffusion of the form:
      !  rho*T*Ds/Dt = ... + nab.(K*gradT)
      !        Ds/Dt = ... + K/rho*[del2lnTT+(glnTT)^2]
      !
      ! NB: chix = K/(cp rho) is needed for diffus_chi calculation
      chix = p%rho1*hcond*p%cp1
      call dot(p%glnTT,p%glnTT,g2)
      thdiff = p%rho1*hcond * (p%del2lnTT + g2)  
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
        diffus_chi=max(diffus_chi,gamma*chix*dxyz_2)
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_simple
!***********************************************************************
    subroutine calc_heatcond_spitzer(df,p)
!
!  Calculates heat conduction parallel and perpendicular (isotropic)
!  to magnetic field lines     
!
!  See: Solar MHD; Priest 1982
!
!  10-feb-04/bing: coded
!
      use EquationOfState
      use Sub
      use IO
!       
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: gvKpara,gvKperp,tmpv
      real, dimension (nx) :: bb2,thdiff,b1
      real, dimension (nx) :: tmps,quenchfactor,vKpara,vKperp
      integer ::i,j
      type (pencil_case) :: p
!
      intent(in) :: p
      intent(out) :: df
!
      if (pretend_lnTT) call fatal_error("calc_heatcond_spitzer","not implemented when pretend_lnTT = T")
!
!     Calculate variable diffusion coefficients along pencils
!
      call dot2_mn(p%bb,bb2)
      b1=1./max(tiny(bb2),bb2)
!
      vKpara = Kgpara * p%TT**3.5
      vKperp = Kgperp * b1*exp(2*p%lnrho+0.5*p%lnTT)
!
!     limit perpendicular diffusion 
!
      if (maxval(abs(vKpara+vKperp)) .gt. tini) then         
         quenchfactor = vKpara/(vKpara+vKperp)
         vKperp=vKperp*quenchfactor
      endif
!
!     Calculate gradient of variable diffusion coefficients
!
      tmps = 3.5 * vKpara 
      call multsv_mn(tmps,p%glnTT,gvKpara)
!
      do i=1,3
         tmpv(:,i)=0.
         do j=1,3
            tmpv(:,i)=tmpv(:,i) + p%bb(:,j)*p%bij(:,j,i)
         end do
      end do
      call multsv_mn(2*b1*b1,tmpv,tmpv)
      tmpv=2.*p%glnrho+0.5*p%glnTT-tmpv
      call multsv_mn(vKperp,tmpv,gvKperp)
      gvKperp=gvKperp*spread(quenchfactor,2,3)
!
!     Calculate diffusion term
!
      call tensor_diffusion_coef(p%glnTT,p%hlnTT,p%bij,p%bb,vKperp,vKpara,thdiff,GVKPERP=gvKperp,GVKPARA=gvKpara)
!
      thdiff = thdiff*exp(-p%lnrho-p%lnTT) 
! 
      df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) + thdiff 
!
      if (lfirst.and.ldt) then
         !
         dt1_max=max(dt1_max,maxval(abs(thdiff)*gamma)/(cdts))
      endif
!      
    endsubroutine calc_heatcond_spitzer
!***********************************************************************
    subroutine calc_heatcond(f,df,p)
!
!  heat conduction
!
!  17-sep-01/axel: coded
!  14-jul-05/axel: corrected expression for chi_t diffusion.
!  30-mar-06/ngrs: simplified calculations using p%glnTT and p%del2lnTT
!
      use Cdata
      use Sub
      use IO
      use Gravity
!     use FArrayManager
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: glnThcond,glhc,glchit_prof !,glnT
      real, dimension (nx) :: chix
      real, dimension (nx) :: thdiff,g2
      real, dimension (nx) :: hcond,chit_prof
      real :: z_prev=-1.23e20
!
      save :: z_prev, hcond, glhc, chit_prof, glchit_prof
!
      intent(in) :: p
      intent(out) :: df
!
      if (pretend_lnTT) call fatal_error("calc_heatcond","not implemented when pretend_lnTT = T")
!
!  Heat conduction / entropy diffusion
!
      if (hcond0 == 0) then
        chix = 0
        thdiff = 0
        hcond = 0
        glhc = 0
      else
        if (headtt) then
          print*,'calc_heatcond: hcond0=',hcond0
          print*,'calc_heatcond: lgravz=',lgravz
          if (lgravz) print*,'calc_heatcond: Fbot,Ftop=',Fbot,Ftop
        endif

        if (lgravz) then
          ! For vertical geometry, we only need to calculate this
          ! for each new value of z -> speedup by about 8% at 32x32x64
          if (z(n) /= z_prev) then
            if (lhcond_global) then
              hcond=f(l1:l2,m,n,iglobal_hcond)
              glhc=f(l1:l2,m,n,iglobal_glhc:iglobal_glhc+2)
            else
              call heatcond(p,hcond)
              call gradloghcond(p,glhc)
            endif
            if (chi_t/=0) then
              call chit_profile(chit_prof)
              call gradlogchit_profile(glchit_prof)
            endif
            z_prev = z(n)
          endif
        else
          if (lhcond_global) then
            hcond=f(l1:l2,m,n,iglobal_hcond)
            glhc=f(l1:l2,m,n,iglobal_glhc:iglobal_glhc+2)
          else
            call heatcond(p,hcond)
            call gradloghcond(p,glhc)
          endif
          if (chi_t/=0) then
            call chit_profile(chit_prof)
            call gradlogchit_profile(glchit_prof)
          endif
        endif
        !
        ! NB: the following left in for the record, but the version below,
        !     using del2lnTT & glnTT, is simpler
        !
        !chix = p%rho1*hcond*p%cp1
        !glnT = gamma*p%gss*spread(p%cp1,2,3) + gamma1*p%glnrho        ! grad ln(T)
        !glnThcond = glnT + glhc/spread(hcond,2,3)    ! grad ln(T*hcond)
        !call dot(glnT,glnThcond,g2)
        !thdiff = p%rho1*hcond * (gamma*p%del2ss*p%cp1 + gamma1*p%del2lnrho + g2)  

        !  diffusion of the form:
        !  rho*T*Ds/Dt = ... + nab.(K*gradT)
        !        Ds/Dt = ... + K/rho*[del2lnTT+(glnTT+glnhcond).glnTT]
        !
        ! NB: chix = K/(cp rho) is needed for diffus_chi calculation
        chix = p%rho1*hcond*p%cp1
        glnThcond = p%glnTT + glhc/spread(hcond,2,3)    ! grad ln(T*hcond)
        call dot(p%glnTT,glnThcond,g2)
        thdiff = p%rho1*hcond * (p%del2lnTT + g2)  
      endif  ! hcond0/=0
!
!  write z-profile (for post-processing)
!
      if (headt .and. lfirst .and. lgravz) call write_zprof('hcond',hcond)
!
!  Write radiative flux array
!
      if (ldiagnos) then
        if (idiag_fradz/=0) call xysum_mn_name_z(-hcond*p%TT*p%glnTT(:,3),idiag_fradz)
        if (idiag_fturbz/=0) call xysum_mn_name_z(-chi_t*p%rho*p%TT*p%gss(:,3),idiag_fturbz)
      endif
!
!  "turbulent" entropy diffusion
!  should only be present if g.gradss > 0 (unstable stratification)
!
      if (chi_t/=0.) then
        if (headtt) then
          print*,'calc_headcond: "turbulent" entropy diffusion: chi_t=',chi_t
          if (hcond0 /= 0) then
            call warning('calc_heatcond',"hcond0 and chi_t combined don't seem to make sense")
          endif
        endif
        ! NB: the cp factor was previously omitted in the following calculation
        !     of glnP (=glnrho+glnTT) -- but we might as well use glnTT 
        !     directly now, since that pencil is now calculated.
        !glnP=gamma*(p%gss/cp+p%glnrho)
        !call dot(glnP+glchit_prof,p%gss,g2)
        call dot(p%glnrho+p%glnTT+glchit_prof,p%gss,g2)
        !thdiff=thdiff+chi_t*(p%del2ss+g2)
        thdiff=thdiff+chi_t*chit_prof*(p%del2ss+g2)
      endif
!
!  check for NaNs initially
!
      if (headt .and. (hcond0 /= 0)) then
        if (notanumber(glhc))      print*,'calc_heatcond: NaNs in glhc'
        if (notanumber(p%rho1))    print*,'calc_heatcond: NaNs in rho1'
        if (notanumber(hcond))     print*,'calc_heatcond: NaNs in hcond'
        if (notanumber(chix))      print*,'calc_heatcond: NaNs in chix'
        if (notanumber(p%del2ss))    print*,'calc_heatcond: NaNs in del2ss'
        if (notanumber(p%del2lnrho)) print*,'calc_heatcond: NaNs in del2lnrho'
        if (notanumber(glhc))      print*,'calc_heatcond: NaNs in glhc'
        if (notanumber(1/hcond))   print*,'calc_heatcond: NaNs in 1/hcond'
        if (notanumber(p%glnTT))      print*,'calc_heatcond: NaNs in glnT'
        if (notanumber(glnThcond)) print*,'calc_heatcond: NaNs in glnThcond'
        if (notanumber(g2))        print*,'calc_heatcond: NaNs in g2'
        if (notanumber(thdiff))    print*,'calc_heatcond: NaNs in thdiff'
        !
        !  most of these should trigger the following trap
        !
        if (notanumber(thdiff)) then
          print*, 'calc_heatcond: m,n,y(m),z(n)=',m,n,y(m),z(n)
          call fatal_error('calc_heatcond','NaNs in thdiff')
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
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond
!***********************************************************************
    subroutine calc_heat_cool(df,p,Hmax)
!
!  add combined heating and cooling
!
!  02-jul-02/wolf: coded
!
      use Cdata
      use Sub
      use Gravity
      use IO
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: Hmax
!
      real, dimension (nx) :: heat,prof,theta_profile
      real :: ssref,zbot,ztop,profile_buffer,xi,profile_cor
!
      intent(in) :: p
      intent(out) :: df
!
      if (pretend_lnTT) call fatal_error("calc_heat_cool","not implemented when pretend_lnTT = T")
!
!  identifier
!
      if (headtt) print*,'calc_heat_cool: lgravz=',lgravz
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
!AB: Wolfgang, the last part of above comment seems wrong; 
!AB: We do divide by rho and T. But what about the heating profile?
!
! heating profile, normalised, so volume integral = 1
        prof = spread(exp(-0.5*((z(n)-zbot)/wheat)**2), 1, l2-l1+1) &
             /(sqrt(pi/2.)*wheat*Lx*Ly)
        heat = luminosity*prof
! smoothly switch on heating if required
        if ((ttransient > 0) .and. (t < ttransient)) then
          heat = heat * t*(2*ttransient-t)/ttransient**2
        endif
!AB: is ssref used anywhere?
        ssref = ss0 + (-log(gamma) + log(cs20))/gamma + grads0*ztop
!
!  allow for different cooling profile functions
!  The gaussian default is rather broad and disturbs the entire interior
!
        if (headtt) print*, 'cooling_profile: cooling_profile,z2,wcool=', &
            cooling_profile, z2, wcool
        select case(cooling_profile)
        case ('gaussian')
          prof = spread(exp(-0.5*((ztop-z(n))/wcool)**2), 1, l2-l1+1)
        case ('step')
          prof = step(spread(z(n),1,nx),z2,wcool)
        case ('cubic_step')
          prof = cubic_step(spread(z(n),1,nx),z2,wcool)
        endselect
!
!  write out (during first time step only) and apply
!
        call write_zprof('cooling_profile',prof)
        heat = heat - cool*prof*(p%cs2-cs2cool)/cs2cool
!
!  Write divergence of cooling flux
!
        if (ldiagnos) then
          if (idiag_dcoolz/=0) call xysum_mn_name_z(heat,idiag_dcoolz)
        endif
      endif
!
!  Spherical case:
!  heat at centre, cool outer layers
!
      if (lgravr) then
        ! central heating
        ! heating profile, normalised, so volume integral = 1
        prof = exp(-0.5*(p%r_mn/wheat)**2) * (2*pi*wheat**2)**(-1.5)
        heat = luminosity*prof
        if (headt .and. lfirst .and. ip<=9) &
          call output_pencil(trim(directory)//'/heat.dat',heat,1)
        ! surface cooling; entropy or temperature
        ! cooling profile; maximum = 1
!        prof = 0.5*(1+tanh((r_mn-1.)/wcool))
        if (rcool==0.) rcool=r_ext
        prof = step(p%r_mn,rcool,wcool)
        !
        !  pick type of cooling
        !
        select case(cooltype)
        case ('cs2', 'Temp')    ! cooling to reference temperature cs2cool
          heat = heat - cool*prof*(p%cs2-cs2cool)/cs2cool
        case ('cs2-rho', 'Temp-rho') ! cool to reference temperature cs2cool
                                     ! in a more time-step neutral manner
          heat = heat - cool*prof*(p%cs2-cs2cool)/cs2cool/p%rho1
        case ('entropy')        ! cooling to reference entropy (currently =0)
          heat = heat - cool*prof*(p%ss-0.)
! dgm
        case ('shell')          !  heating/cooling at shell boundaries
          heat=0.                            ! default
          select case(initss(1))
!           case ('geo-kws','shell_layers'); heat=0.        ! can add heating later based on value of initss
            case ('geo-kws'); heat=0.    ! can add heating later based on value of initss
          endselect
          !
          !  possibility of a latitudinal heating profile
          !  T=T0-(2/3)*delT*P2(costheta), for testing Taylor-Proudman theorem
          !  Note that P2(x)=(1/2)*(3*x^2-1).
          !
          if (deltaT_poleq/=0.) then
            if (headtt) print*,'calc_heat_cool: deltaT_poleq=',deltaT_poleq
            if (headtt) print*,'p%rcyl_mn=',p%rcyl_mn
            if (headtt) print*,'p%z_mn=',p%z_mn
            theta_profile=(1./3.-(p%rcyl_mn/p%z_mn)**2)*deltaT_poleq
            prof = step(p%r_mn,r_ext,wcool)      ! outer heating/cooling step
            heat = heat - cool_ext*prof*(p%cs2-cs2_ext)/cs2_ext*theta_profile
            prof = 1 - step(p%r_mn,r_int,wcool)  ! inner heating/cooling step
            heat = heat - cool_int*prof*(p%cs2-cs2_int)/cs2_int*theta_profile
          else
            prof = step(p%r_mn,r_ext,wcool)      ! outer heating/cooling step
            heat = heat - cool_ext*prof*(p%cs2-cs2_ext)/cs2_ext
            prof = 1 - step(p%r_mn,r_int,wcool)  ! inner heating/cooling step
            heat = heat - cool_int*prof*(p%cs2-cs2_int)/cs2_int
          endif
!
        case default
          write(unit=errormsg,fmt=*) &
               'calc_heat_cool: No such value for cooltype: ', trim(cooltype)
          call fatal_error('calc_heat_cool',errormsg)
        endselect
      endif
!
!  add spatially uniform heating (usually as a test)
!
      if (heat_uniform/=0.) heat=heat+heat_uniform
!
!  add "coronal" heating (to simulate a hot corona)
!  assume a linearly increasing reference profile
!  This 1/rho1 business is clumpsy, but so would be obvious alternatives...
!
      if (tau_cor>0) then
        if (z(n)>=z_cor) then
          xi=(z(n)-z_cor)/(ztop-z_cor)
          profile_cor=xi**2*(3-2*xi)
          heat=heat+profile_cor*(TT_cor-1/p%TT1)/(p%rho1*tau_cor*p%cp1)
        endif
      endif
!
!  add heating and cooling to a reference temperature in a buffer
!  zone at the z boundaries. Only regions in |z| > zheat_buffer are affected.
!  Inverse width of the transition is given by dheat_buffer1.
!
      if (tauheat_buffer/=0.) then
        profile_buffer=0.5*(1.+tanh(dheat_buffer1*(z(n)-zheat_buffer)))
        !profile_buffer=0.5*(1.+tanh(dheat_buffer1*(z(n)**2-zheat_buffer**2)))
!       profile_buffer=1.+0.5*(tanh(dheat_buffer1*(z(n)-z(n1)-zheat_buffer)) + tanh(dheat_buffer1*(z(n)-z(n2)-zheat_buffer)))
        heat=heat+profile_buffer*p%ss* &
            (TTheat_buffer-1/p%TT1)/(p%rho1*tauheat_buffer)
      endif
!
!  add to entropy equation
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + p%TT1*p%rho1*heat
      if (lfirst.and.ldt) Hmax=Hmax+heat*p%rho1
!
    endsubroutine calc_heat_cool
!***********************************************************************
    subroutine calc_heat_cool_RTV(df,p)
!
!    calculate cool term:  C = ne*ni*Q(T) 
!    with ne*ni = 1.2*np^2 = 1.2*rho^2/(1.4*mp)^2
!    Q(T) = H*T^B is piecewice poly
!    [Q] = [v]^3 [rho] [l]^5
!
!  15-dec-04/bing: coded
!
      use IO
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: lnQ,rtv_cool,lnTT_SI,lnneni
      integer :: i,imax
      real :: unit_lnQ
      type (pencil_case) :: p
!
      intent(in) :: p
      intent(out) :: df
!
      if (pretend_lnTT) call fatal_error("calc_heat_cool_RTV","not implemented when pretend_lnTT = T")
!
!     All is in SI units and has to be rescaled to PENCIL units
!
      unit_lnQ=3*alog(real(unit_velocity))+5*alog(real(unit_length))+alog(real(unit_density))
      if (unit_system .eq. 'cgs') unit_lnQ = unit_lnQ+alog(1.e13)
!
      lnTT_SI = p%lnTT + alog(real(unit_temperature)) 
!
!     calculate ln(ne*ni) :
!          ln(ne*ni) = ln( 1.2*rho^2/(1.4*mp)^2)
      lnneni = 2*p%lnrho + alog(1.2) - 2*alog(1.4*real(m_p))
!      
!     rtv_cool=exp(lnneni+lnQ-unit_lnQ-lnrho-lnTT)
!   
! First set of parameters
      rtv_cool=0.
      if (cool_RTV .gt. 0.) then
        imax = size(intlnT_1,1)       
        lnQ(:)=0.0
        do i=1,imax-1
          where (( intlnT_1(i) <= lnTT_SI .or. i==1 ) .and. lnTT_SI < intlnT_1(i+1) )
            lnQ=lnQ + lnH_1(i) + B_1(i)*lnTT_SI
          endwhere
        enddo
        where (lnTT_SI >= intlnT_1(imax) )
          lnQ = lnQ + lnH_1(imax-1) + B_1(imax-1)*intlnT_1(imax)
        endwhere
        rtv_cool=exp(lnneni+lnQ-unit_lnQ-p%lnTT-p%lnrho)
      elseif (cool_RTV .lt. 0) then
! Second set of parameters      
        cool_RTV = cool_RTV*(-1.)
        imax = size(intlnT_2,1)       
        lnQ(:)=0.0
        do i=1,imax-1
          where (( intlnT_2(i) <= lnTT_SI .or. i==1 ) .and. lnTT_SI < intlnT_2(i+1) )
            lnQ=lnQ + lnH_2(i) + B_2(i)*lnTT_SI
          endwhere
        enddo
        where (lnTT_SI >= intlnT_2(imax) )
          lnQ = lnQ + lnH_2(imax-1) + B_2(imax-1)*intlnT_2(imax)
        endwhere
        rtv_cool=exp(lnneni+lnQ-unit_lnQ-p%lnTT-p%lnrho)
      else
        rtv_cool(:)=0.
      endif
!
      rtv_cool=rtv_cool * cool_RTV  ! for adjusting by setting cool_RTV in run.in
!
!     add to entropy equation
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss)-rtv_cool 
!
      if (lfirst.and.ldt) then
         !
         dt1_max=max(dt1_max,maxval(rtv_cool*gamma)/(cdts))
      endif
!
    endsubroutine calc_heat_cool_RTV
!***********************************************************************
    subroutine calc_tau_ss_exterior(df,p)
!
!  entropy relaxation to zero on time scale tau_ss_exterior within
!  exterior region. For the time being this means z > zgrav.
!
!  29-jul-02/axel: coded
!
      use Cdata
      use Gravity
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real :: scl
!
      intent(in) :: p
      intent(out) :: df
!
      if (pretend_lnTT) call fatal_error("calc_tau_ss_exterior","not implemented when pretend_lnTT = T")
!
      if (headtt) print*,'calc_tau_ss_exterior: tau=',tau_ss_exterior
      if (z(n)>zgrav) then
        scl=1./tau_ss_exterior
        df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)-scl*p%ss
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
      integer :: iname,inamez,inamey,inamex,irz
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
        idiag_dtc=0; idiag_eth=0; idiag_ethdivum=0; idiag_ssm=0
        idiag_eem=0; idiag_ppm=0; idiag_csm=0
        idiag_ugradpm=0; idiag_ethtot=0; idiag_dtchi=0; idiag_ssmphi=0
        idiag_yHmax=0; idiag_yHm=0; idiag_TTmax=0; idiag_TTmin=0; idiag_TTm=0
        idiag_fconvz=0; idiag_dcoolz=0; idiag_fradz=0; idiag_fturbz=0
        idiag_ssmz=0; idiag_ssmy=0; idiag_ssmx=0; idiag_TTmz=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtc',idiag_dtc)
        call parse_name(iname,cname(iname),cform(iname),'dtchi',idiag_dtchi)
        call parse_name(iname,cname(iname),cform(iname),'ethtot',idiag_ethtot)
        call parse_name(iname,cname(iname),cform(iname),'ethdivum',idiag_ethdivum)
        call parse_name(iname,cname(iname),cform(iname),'eth',idiag_eth)
        call parse_name(iname,cname(iname),cform(iname),'ssm',idiag_ssm)
        call parse_name(iname,cname(iname),cform(iname),'eem',idiag_eem)
        call parse_name(iname,cname(iname),cform(iname),'ppm',idiag_ppm)
        call parse_name(iname,cname(iname),cform(iname),'csm',idiag_csm)
        call parse_name(iname,cname(iname),cform(iname),'ugradpm',idiag_ugradpm)
        call parse_name(iname,cname(iname),cform(iname),'yHm',idiag_yHm)
        call parse_name(iname,cname(iname),cform(iname),'yHmax',idiag_yHmax)
        call parse_name(iname,cname(iname),cform(iname),'TTm',idiag_TTm)
        call parse_name(iname,cname(iname),cform(iname),'TTmax',idiag_TTmax)
        call parse_name(iname,cname(iname),cform(iname),'TTmin',idiag_TTmin)
        call parse_name(iname,cname(iname),cform(iname),'TTp',idiag_TTp)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fturbz',idiag_fturbz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fconvz',idiag_fconvz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'dcoolz',idiag_dcoolz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fradz',idiag_fradz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ssmz',idiag_ssmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TTmz',idiag_TTmz)
      enddo
!
!  check for those quantities for which we want xz-averages
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'ssmy',idiag_ssmy)
      enddo
!
!  check for those quantities for which we want yz-averages
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ssmx',idiag_ssmx)
      enddo
!
!  check for those quantities for which we want phi-averages
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'ssmphi',idiag_ssmphi)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_dtc=',idiag_dtc
        write(3,*) 'i_dtchi=',idiag_dtchi
        write(3,*) 'i_ethtot=',idiag_ethtot
        write(3,*) 'i_ethdivum=',idiag_ethdivum
        write(3,*) 'i_eth=',idiag_eth
        write(3,*) 'i_ssm=',idiag_ssm
        write(3,*) 'i_eem=',idiag_eem
        write(3,*) 'i_ppm=',idiag_ppm
        write(3,*) 'i_csm=',idiag_csm
        write(3,*) 'i_ugradpm=',idiag_ugradpm
        write(3,*) 'i_ssmphi=',idiag_ssmphi
        write(3,*) 'i_fturbz=',idiag_fturbz
        write(3,*) 'i_fconvz=',idiag_fconvz
        write(3,*) 'i_dcoolz=',idiag_dcoolz
        write(3,*) 'i_fradz=',idiag_fradz
        write(3,*) 'i_ssmz=',idiag_ssmz
        write(3,*) 'i_TTmz=',idiag_TTmz
        write(3,*) 'nname=',nname
        write(3,*) 'iss=',iss
        write(3,*) 'i_yHmax=',idiag_yHmax
        write(3,*) 'i_yHm=',idiag_yHm
        write(3,*) 'i_TTmax=',idiag_TTmax
        write(3,*) 'i_TTmin=',idiag_TTmin
        write(3,*) 'i_TTm=',idiag_TTm
        write(3,*) 'i_TTp=',idiag_TTp
        write(3,*) 'iyH=',iyH
        write(3,*) 'ilnTT=',ilnTT
      endif
!
    endsubroutine rprint_entropy
!***********************************************************************
    subroutine calc_heatcond_zprof(zprof_hcond,zprof_glhc)
!
!  calculate z-profile of heat conduction for multilayer setup
!
!  12-jul-05/axel: coded
!
      use Cdata
      use Sub
      use Gravity, only: z1, z2
!
      real, dimension (nz,3) :: zprof_glhc
      real, dimension (nz) :: zprof_hcond
      real :: zpt
!
      intent(out) :: zprof_hcond,zprof_glhc
!
      do n=1,nz
        zpt=z(n+nghost)
        zprof_hcond(n) = 1 + (hcond1-1)*cubic_step(zpt,z1,-widthss) &
                           + (hcond2-1)*cubic_step(zpt,z2,+widthss)
        zprof_hcond(n) = hcond0*zprof_hcond(n)
        zprof_glhc(n,1:2) = 0.
        zprof_glhc(n,3) = (hcond1-1)*cubic_der_step(zpt,z1,-widthss) &
                        + (hcond2-1)*cubic_der_step(zpt,z2,+widthss)
        zprof_glhc(n,3) = hcond0*zprof_glhc(n,3)
      enddo
!
    endsubroutine calc_heatcond_zprof
!***********************************************************************
    subroutine heatcond(p,hcond)
!
!  calculate the heat conductivity hcond along a pencil.
!  This is an attempt to remove explicit reference to hcond[0-2] from
!  code, e.g. the boundary condition routine.
!
!  NB: if you modify this profile, you *must* adapt gradloghcond below.
!
!  23-jan-2002/wolf: coded
!  18-sep-2002/axel: added lmultilayer switch
!  09-aug-2006/dintrans: added a radial profile hcond(r)
!
      use Sub, only: step
      use Gravity
!
      real, dimension (nx) :: hcond
      type (pencil_case) :: p
!
      if (lgravz) then
        if (lmultilayer) then
          hcond = 1 + (hcond1-1)*step(p%z_mn,z1,-widthss) &
                    + (hcond2-1)*step(p%z_mn,z2,widthss)
          hcond = hcond0*hcond
        else
          hcond=Kbot
        endif
      else
        if (lmultilayer) then
          hcond = 1. + (hcond1-1.)*step(p%r_mn,r_bcz,-widthss)
          hcond = hcond0*hcond
        else
          hcond = hcond0
        endif
      endif
!
    endsubroutine heatcond
!***********************************************************************
    subroutine gradloghcond(p,glhc)
!
!  calculate grad(log hcond), where hcond is the heat conductivity
!  NB: *Must* be in sync with heatcond() above.
!
!  23-jan-2002/wolf: coded
!
      use Sub, only: der_step
      use Gravity
!
      real, dimension (nx,3) :: glhc
      real, dimension (nx)   :: dhcond
      type (pencil_case) :: p
!
      if (lgravz) then
        if (lmultilayer) then
          glhc(:,3) = (hcond1-1)*der_step(p%z_mn,z1,-widthss) &
                      + (hcond2-1)*der_step(p%z_mn,z2,widthss)
          glhc(:,3) = hcond0*glhc(:,3)
        else
          glhc = 0.
        endif
      else
        if (lmultilayer) then
          dhcond=hcond0*(hcond1-1.)*der_step(p%r_mn,r_bcz,-widthss)
          glhc(:,1) = x(l1:l2)/p%r_mn*dhcond
          glhc(:,2) = y(m)/p%r_mn*dhcond
          if (lcylindrical) then 
            glhc(:,3) = 0.
          else
            glhc(:,3) = z(n)/p%r_mn*dhcond
          endif
        else
          glhc = 0.
        endif
      endif
!
    endsubroutine gradloghcond
!***********************************************************************
    subroutine chit_profile(chit_prof)
!
!  calculate the chit_profile conductivity chit_prof along a pencil.
!  This is an attempt to remove explicit reference to chit_prof[0-2] from
!  code, e.g. the boundary condition routine.
!
!  NB: if you modify this profile, you *must* adapt gradlogchit_prof below.
!
!  23-jan-2002/wolf: coded
!  18-sep-2002/axel: added lmultilayer switch
!
      use Sub, only: step
      use Gravity
!
      real, dimension (nx) :: chit_prof,z_mn
!
      if (lgravz) then
        if (lmultilayer) then
          z_mn=spread(z(n),1,nx)
          chit_prof = 1 + (chit_prof1-1)*step(z_mn,z1,-widthss) &
                        + (chit_prof2-1)*step(z_mn,z2,widthss)
        else
          chit_prof=1.
        endif
      else
        chit_prof=1.
      endif
!
    endsubroutine chit_profile
!***********************************************************************
    subroutine gradlogchit_profile(glchit_prof)
!
!  calculate grad(log chit_prof), where chit_prof is the heat conductivity
!  NB: *Must* be in sync with heatcond() above.
!
!  23-jan-2002/wolf: coded
!
      use Sub, only: der_step
      use Gravity
!
      real, dimension (nx,3) :: glchit_prof
      real, dimension (nx) :: z_mn
!
      if (lgravz) then
        if (lmultilayer) then
          z_mn=spread(z(n),1,nx)
          glchit_prof(:,1:2) = 0.
          glchit_prof(:,3) = (chit_prof1-1)*der_step(z_mn,z1,-widthss) &
                           + (chit_prof2-1)*der_step(z_mn,z2,widthss)
        else
          glchit_prof = 0.
        endif
      else
        glchit_prof = 0.
      endif
!
    endsubroutine gradlogchit_profile
!***********************************************************************
    subroutine newton_cool(df,p)
!
!  Keeps the temperature in the lower chromosphere 
!  at a constant level using newton cooling
!
!  15-dec-2004/bing: coded
!  25-sep-2006/bing: updated, using external data
!
      use EquationOfState, only: lnrho0,gamma
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: newton
      real, dimension (150), save :: b_lnT,b_z  
      real :: lnTTor
      integer :: i,lend
      type (pencil_case) :: p
!
      intent(in) :: p
      intent(inout) :: df
!
      if (pretend_lnTT) call fatal_error("newton_cool","not implemented when pretend_lnTT = T")
!
!  Initial temperature profile is given in ln(T) in [K] over z in [Mm] 
!
      if (it .eq. 1) then
         inquire(IOLENGTH=lend) lnTTor
         open (10,file='driver/b_lnT.dat',form='unformatted',status='unknown',recl=lend*150)
         read (10) b_lnT
         read (10) b_z
         close (10) 
         !
         b_lnT = b_lnT - alog(real(unit_temperature))
         if (unit_system == 'SI') then
           b_z = b_z * 1.e6 / unit_length 
         elseif (unit_system == 'cgs') then
           b_z = b_z * 1.e8 / unit_length 
         endif
      endif
!
!  Get reference temperature
!
      if (z(n) .lt. b_z(1) ) then
        lnTTor = b_lnT(1)
      elseif (z(n) .ge. b_z(150)) then
        lnTTor = b_lnT(150)
      else
        do i=1,149
          if (z(n) .ge. b_z(i) .and. z(n) .lt. b_z(i+1)) then 
            !
            ! linear interpolation
            !
            lnTTor = (b_lnT(i)*(b_z(i+1) - z(n)) +   &
                b_lnT(i+1)*(z(n) - b_z(i)) ) / (b_z(i+1)-b_z(i))
            exit
          endif
        enddo
      endif
!
      if (dt .ne. 0 )  then 
         newton = tdown/gamma/dt*((lnTTor - p%lnTT) ) 
         newton = newton * exp(allp*(-abs(p%lnrho-lnrho0)))
!
!  Add newton cooling term to entropy
!
         df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + newton
!
         if (lfirst.and.ldt) then
            !
            dt1_max=max(dt1_max,maxval(abs(newton)*gamma)/(cdts))
         endif
      endif
!
    endsubroutine newton_cool
!***********************************************************************
    subroutine strat_MLT (rhotop, cs2top, mixinglength_flux, lnrhom, &
                     ssm, cs2m, rhobot)
!
! 04-mai-2006/boris: called by 'mixinglength' to iterate the MLT
! equations until rho=1 at the bottom of convection zone (z=z1)
!
      use Cdata
      use Gravity, only: z1,z2,gravz
      use EquationOfState, only: gamma, gamma1
!
      real, dimension (nzgrid) :: lnrhom, ssm, cs2m, zz, eem
      real :: rhotop, cs2top, zm, ztop, dlnrho, dee, &
              mixinglength_flux, lnrhobot, rhobot
      real :: del, delad, fr_frac, fc_frac, fc, polyad=3./2.
      integer :: nbot1, nbot2
!
!  inital values at the top
!
      lnrhom(1)=alog(rhotop)
      cs2m(1)=cs2top
      eem(1)=cs2top/gamma/gamma1
      ssm(1)=(alog(cs2m(1))-gamma1*lnrhom(1))/gamma
      ztop=xyz0(3)+Lxyz(3)
      zz(1)=ztop
!
      delad=1.-1./gamma
      fr_frac=delad*(mpoly0+1.)
      fc_frac=1.-fr_frac
      fc=fc_frac*mixinglength_flux
!   print*,'fr_frac, fc_frac, fc=',fr_frac,fc_frac,fc
!
      do iz=2,nzgrid
        zm=ztop-(iz-1)*dz
        zz(iz)=zm
        if (zm<=z1) then
! radiative zone=polytropic stratification
          del=1./(mpoly1+1.)
        else
          if (zm<=z2) then
! convective zone=mixing-length stratification
!         del=delad
            del=delad+1.5*(fc/ &
                (exp(lnrhom(iz-1))*cs2m(iz-1)**1.5))**.6666667
          else
! upper zone=isothermal stratification
            del=0.
          endif
        endif
        dee=-polyad*del
        dlnrho=-polyad*(1.-del)/eem(iz-1)
        eem(iz)=eem(iz-1)-dee*dz
        lnrhom(iz)=lnrhom(iz-1)-dlnrho*dz
        cs2m(iz)=gamma*gamma1*eem(iz)
        ssm(iz)=(alog(cs2m(iz))-gamma1*lnrhom(iz))/gamma
      enddo
!
!  find the value of rhobot
!
      do iz=1,nzgrid
        if (zz(iz)<z1) exit
      enddo
!     stop 'find rhobot: didnt find bottom value of z'
      nbot1=iz-1
      nbot2=iz
!
!  interpolate
!
      lnrhobot=lnrhom(nbot1)+(lnrhom(nbot2)-lnrhom(nbot1))/ &
               (zz(nbot2)-zz(nbot1))*(z1-zz(nbot1))
      rhobot=exp(lnrhobot) 
!   print*,'find rhobot=',rhobot
!
    endsubroutine strat_MLT
!***********************************************************************
    subroutine shell_ss_layers(f)
!
!  Initialize entropy in a spherical shell using two polytropic layers
!
!  09-aug-06/dintrans: coded
!
      use Gravity, only: g0
      use EquationOfState, only: eoscalc, ilnrho_lnTT, mpoly0, mpoly1, lnrho0

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: lnrho,lnTT,TT,ss,r_mn
      real :: beta0,beta1,TT_crit
      real :: lnrho_int,lnrho_ext,lnrho_crit

      if (headtt) print*,'r_bcz in entropy.f90=',r_bcz
!
!  beta is the temperature gradient
!  1/beta = -(g/cp) 1./[(1-1/gamma)*(m+1)]
!  AB: Boris, did you ignore cp below?
!
      beta0=-g0/(mpoly0+1)*gamma/gamma1
      beta1=-g0/(mpoly1+1)*gamma/gamma1
      TT_crit=TT_ext+beta0*(r_bcz-r_ext)
      lnrho_ext=lnrho0
      lnrho_crit=lnrho0+ &
            mpoly0*log(TT_ext+beta0*(r_bcz-r_ext))-mpoly0*log(TT_ext)
      lnrho_int=lnrho_crit + &
            mpoly1*log(TT_crit+beta1*(r_int-r_bcz))-mpoly1*log(TT_crit)
!
!  set initial condition
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
!
        r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
!
!  convective layer
!
        where (r_mn < r_ext .AND. r_mn > r_bcz)
          TT = TT_ext+beta0*(r_mn-r_ext)
          f(l1:l2,m,n,ilnrho)=lnrho0+mpoly0*log(TT)-mpoly0*log(TT_ext)
        endwhere
!
!  radiative layer
!
        where (r_mn <= r_bcz .AND. r_mn > r_int)
          TT = TT_crit+beta1*(r_mn-r_bcz)
          f(l1:l2,m,n,ilnrho)=lnrho_crit+mpoly1*log(TT)-mpoly1*log(TT_crit)
        endwhere
        where (r_mn >= r_ext) 
          TT = TT_ext
          f(l1:l2,m,n,ilnrho)=lnrho_ext
        endwhere
        where (r_mn <= r_int) 
          TT = TT_int
          f(l1:l2,m,n,ilnrho)=lnrho_int
        endwhere
!
        lnrho=f(l1:l2,m,n,ilnrho)
        lnTT=log(TT)
        call eoscalc(ilnrho_lnTT,lnrho,lnTT,ss=ss)
        f(l1:l2,m,n,iss)=ss
!
      enddo 
!      
    endsubroutine shell_ss_layers
!***********************************************************************
endmodule Entropy
