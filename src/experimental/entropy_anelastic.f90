! $Id$
!
!  This module takes care of energy (initial condition
!  and time advance)
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lentropy = .true.
! CPARAM logical, parameter :: ltemperature = .false.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 1
! COMMUNICATED AUXILIARIES 1
!
! PENCILS PROVIDED ugss; Ma2; fpres(3); uglnTT
! PENCILS PROVIDED ugss_b; gss(3); del2ss; glnTTb(3) ; gss_b(3)
!
!***************************************************************
module Energy
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use EquationOfState, only: gamma, gamma_m1, gamma1, cs20, cs2top, cs2bot
  use Density, only: beta_glnrho_global
  use Interstellar
  use Messages
  use Viscosity
!
  implicit none
!
  include '../energy.h'
!
  real :: radius_ss=0.1,ampl_ss=0.,widthss=2*epsi,epsilon_ss=0.
  real :: luminosity=0.,wheat=0.1,cool=0.,rcool=0.,wcool=0.1
  real :: TT_int,TT_ext,cs2_int,cs2_ext,cool_int=0.,cool_ext=0.,ampl_TT=0.
  real,target :: chi=0.
  real :: chi_t=0.,chi_shock=0.,chi_hyper3=0.
  real :: Kgperp=0.,Kgpara=0.,tdown=0.,allp=2.
  real :: ss_left=1.,ss_right=1.
  real :: ss0=0.,khor_ss=1.,ss_const=0.
  real :: pp_const=0.
  real :: tau_ss_exterior=0.,T0=1.
  real :: mixinglength_flux=0.
  real :: cp1=0.0
  !parameters for Sedov type initial condition
  real :: center1_x=0., center1_y=0., center1_z=0.
  real :: center2_x=0., center2_y=0., center2_z=0.
  real :: kx_ss=1.,ky_ss=1.,kz_ss=1.
  real :: thermal_background=0., thermal_peak=0., thermal_scaling=1.
  real, target :: cs2cool
  real :: cool_fac=1., chiB=0.
  real, dimension(3) :: chi_hyper3_aniso=0.
!
  real, target :: hcond0=impossible,hcond1=impossible
  real, target :: Fbot=impossible,FbotKbot=impossible
  real, target :: Ftop=impossible,FtopKtop=impossible
!
  real :: Kbot=impossible,Ktop=impossible
  real :: hcond2=impossible
  real :: chit_prof1=1.,chit_prof2=1.
  real :: tau_cor=0.,TT_cor=0.,z_cor=0.
  real :: tauheat_buffer=0.,TTheat_buffer=0.,zheat_buffer=0.,dheat_buffer1=0.
  real :: heat_uniform=0.,cool_RTV=0.
  real :: deltaT_poleq=0.,beta_hand=1.,r_bcz=0.
  real :: tau_cool=0.0, TTref_cool=0.0, tau_cool2=0.0
  real :: cs0hs=0.0,H0hs=0.0,rho0hs=0.0
  real :: chit_aniso=0.0,xbot=0.0
  real, pointer :: mpoly
  real, target :: mpoly0=1.5, mpoly1=1.5, mpoly2=1.5
  integer, parameter :: nheatc_max=4
  logical :: lturbulent_heat=.false.
  logical :: lheatc_Kprof=.false.,lheatc_Kconst=.false.
  logical, target :: lheatc_chiconst=.false.
  logical :: lheatc_tensordiffusion=.false.,lheatc_spitzer=.false.
  logical :: lheatc_hubeny=.false.
  logical :: lheatc_corona=.false.
  logical :: lheatc_shock=.false.,lheatc_hyper3ss=.false.
  logical :: lheatc_hyper3ss_polar=.false.,lheatc_hyper3ss_aniso=.false.
  logical :: lupw_ss=.false.
  logical, target :: lmultilayer=.true.
  logical :: ladvection_entropy=.true.
  logical :: lviscosity_heat=.true.
  logical :: lfreeze_sint=.false.,lfreeze_sext=.false.
  logical :: lhcond_global=.false.
!
  integer :: iglobal_hcond=0
  integer :: iglobal_glhc=0
!
  character (len=labellen), dimension(ninit) :: initss='nothing'
  character (len=labellen) :: borderss='nothing'
  character (len=labellen) :: pertss='zero'
  character (len=labellen) :: cooltype='Temp',cooling_profile='gaussian'
  character (len=labellen), dimension(nheatc_max) :: iheatcond='nothing'
  character (len=intstr) :: iinit_str
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
!
  ! input parameters
  namelist /entropy_init_pars/ &
      initss,     &
      pertss,     &
      grads0,     &
      radius_ss,  &
      ampl_ss,    &
      widthss,    &
      epsilon_ss, &
      mixinglength_flux, &
      chi_t, &
      pp_const, &
      ss_left,ss_right,ss_const,mpoly0,mpoly1,mpoly2, &
      khor_ss,thermal_background,thermal_peak,thermal_scaling,cs2cool, &
      center1_x, center1_y, center1_z, center2_x, center2_y, center2_z, &
      T0,ampl_TT,kx_ss,ky_ss,kz_ss,beta_glnrho_global,ladvection_entropy, &
      lviscosity_heat, &
      r_bcz,luminosity,wheat,hcond0,tau_cool,TTref_cool,lhcond_global, &
      cool_fac,cs0hs,H0hs,rho0hs
!
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
      lturbulent_heat,deltaT_poleq, &
      tdown, allp,beta_glnrho_global,ladvection_entropy, &
      lviscosity_heat,r_bcz,lfreeze_sint,lfreeze_sext,lhcond_global, &
      tau_cool,TTref_cool,mixinglength_flux,chiB,chi_hyper3_aniso, Ftop, &
      chit_aniso,xbot
!
  ! diagnostic variables (need to be consistent with reset list below)
  integer :: idiag_dtc=0        ! DIAG_DOC: $\delta t/[c_{\delta t}\,\delta_x
                                ! DIAG_DOC:   /\max c_{\rm s}]$
                                ! DIAG_DOC:   \quad(time step relative to
                                ! DIAG_DOC:   acoustic time step;
                                ! DIAG_DOC:   see \S~\ref{time-step})
  integer :: idiag_ethm=0       ! DIAG_DOC: $\left<\varrho e\right>$
                                ! DIAG_DOC:   \quad(mean thermal
                                ! DIAG_DOC:   [=internal] energy)
  integer :: idiag_ethdivum=0   ! DIAG_DOC:
  integer :: idiag_ssm=0        ! DIAG_DOC: $\left<s/c_p\right>$
                                ! DIAG_DOC:   \quad(mean entropy)
  integer :: idiag_ss2m=0       ! DIAG_DOC: $\left<(s/c_p)^2\right>$
                                ! DIAG_DOC:   \quad(mean squared entropy)
  integer :: idiag_eem=0        ! DIAG_DOC: $\left<e\right>$
  integer :: idiag_ppm=0        ! DIAG_DOC: $\left<p\right>$
  integer :: idiag_csm=0        ! DIAG_DOC: $\left<c_{\rm s}\right>$
  integer :: idiag_pdivum=0     ! DIAG_DOC: $\left<p\nabla\uv\right>$
  integer :: idiag_heatm=0      ! DIAG_DOC:
  integer :: idiag_ugradpm=0    ! DIAG_DOC:
  integer :: idiag_fradbot=0    ! DIAG_DOC: $\int F_{\rm bot}\cdot d\vec{S}$
  integer :: idiag_fradtop=0    ! DIAG_DOC: $\int F_{\rm top}\cdot d\vec{S}$
  integer :: idiag_TTtop=0      ! DIAG_DOC: $\int T_{\rm top} d\vec{S}$
  integer :: idiag_ethtot=0     ! DIAG_DOC: $\int_V\varrho e\,dV$
                                ! DIAG_DOC:   \quad(total thermal
                                ! DIAG_DOC:   [=internal] energy)
  integer :: idiag_dtchi=0      ! DIAG_DOC: $\delta t / [c_{\delta t,{\rm v}}\,
                                ! DIAG_DOC:   \delta x^2/\chi_{\rm max}]$
                                ! DIAG_DOC:   \quad(time step relative to time
                                ! DIAG_DOC:   step based on heat conductivity;
                                ! DIAG_DOC:   see \S~\ref{time-step})
  integer :: idiag_ssmphi=0     ! PHIAVG_DOC: $\left<s\right>_\varphi$
  integer :: idiag_cs2mphi=0    ! PHIAVG_DOC: $\left<c^2_s\right>_\varphi$
  integer :: idiag_yHm=0        ! DIAG_DOC:
  integer :: idiag_yHmax=0      ! DIAG_DOC:
  integer :: idiag_TTm=0        ! DIAG_DOC:
  integer :: idiag_TTmax=0      ! DIAG_DOC:
  integer :: idiag_TTmin=0      ! DIAG_DOC:
  integer :: idiag_fconvm=0     ! DIAG_DOC:
  integer :: idiag_fconvz=0     ! DIAG_DOC:
  integer :: idiag_dcoolz=0     ! DIAG_DOC:
  integer :: idiag_fradz=0      ! DIAG_DOC:
  integer :: idiag_fturbz=0     ! DIAG_DOC:
  integer :: idiag_ssmx=0       ! DIAG_DOC:
  integer :: idiag_ssmy=0       ! DIAG_DOC:
  integer :: idiag_ssmz=0       ! DIAG_DOC:
  integer :: idiag_TTp=0        ! DIAG_DOC:
  integer :: idiag_ssmr=0       ! DIAG_DOC:
  integer :: idiag_TTmx=0       ! DIAG_DOC:
  integer :: idiag_TTmy=0       ! DIAG_DOC:
  integer :: idiag_TTmz=0       ! DIAG_DOC:
  integer :: idiag_TTmxy=0      ! DIAG_DOC:
  integer :: idiag_TTmxz=0      ! DIAG_DOC:
  integer :: idiag_TTmr=0       ! DIAG_DOC:
  integer :: idiag_uxTTmz=0     ! DIAG_DOC:
  integer :: idiag_uyTTmz=0     ! DIAG_DOC:
  integer :: idiag_uzTTmz=0     ! DIAG_DOC:
  integer :: idiag_uxTTmxy=0    ! DIAG_DOC:
  integer :: idiag_ssmxy=0      ! DIAG_DOC: $\left< s \right>_{z}$
  integer :: idiag_ssmxz=0      ! DIAG_DOC: $\left< s \right>_{y}$
!
!  Auxiliary variables
!
  real, dimension(nx) :: diffus_chi, diffus_chi3
!
  contains
!
!***********************************************************************
    subroutine register_energy()
!
!  initialise variables which should know that we solve an energy
!  equation: iss, etc; increase nvar accordingly
!
!  6-nov-01/wolf: coded
!
      use FArrayManager
      use SharedVariables
      use Sub
!
      call farray_register_pde('ss',iss)
      call farray_register_auxiliary('ss_b',iss_b,communicated=.true.)
!
!  identify version number
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_energy
!***********************************************************************
    subroutine initialize_energy(f)
!
!  Called by run.f90 after reading parameters, but before the time loop.
!
!  21-jul-02/wolf: coded
!
      use BorderProfiles, only: request_border_driving
      use EquationOfState, only: cs0, get_soundspeed, get_cp1, &
                                 select_eos_variable,gamma,gamma_m1
      use Density, only: beta_glnrho_global, beta_glnrho_scaled
      use FArrayManager
      use Gravity, only: gravz,g0
      use Mpicomm, only: stop_it
      use SharedVariables, only: get_shared_variable, put_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (nx,3) :: glhc
      real, dimension (nx) :: hcond
      real :: beta1, beta0, TT_crit
      integer :: i, q
      logical :: lnothing,lcompute_grav
      type (pencil_case) :: p
!
! Check any module dependencies
!
      if (.not. leos) &
        call fatal_error('initialize_energy', &
            'EOS=noeos but energy requires an EQUATION OF STATE for the fluid')
!
! Tell the equation of state that we're here and what f variable we use
!
!ajwm      if (.not. lreloading) then ! already in place when reloading
      if (pretend_lnTT) then
        call select_eos_variable('lnTT',iss)
      else
        if (gamma_m1==0.) then
          call fatal_error('initialize_energy',&
               'Use experimental/noenergy for isothermal case')
        else
          call select_eos_variable('ss',iss)
        endif
      endif

      if (ldensity.and..not.lstratz) then
        call get_shared_variable('mpoly',mpoly)
      else
        call warning('initialize_energy','mpoly not obtained from density,'// &
                     'set impossible')
        allocate(mpoly); mpoly=impossible
      endif
!
!ajwm      endif
!
!  radiative diffusion: initialize flux etc
!
      hcond = 0.
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
          call warning('initialize_energy', &
              'You should not set Kbot and hcond0 at the same time')
        endif
      endif
!
!  hcond0 is given by mixinglength_flux in the MLT case
!  hcond1 and hcond2 follow below in the block 'if (lmultilayer) etc...'
!
      if (mixinglength_flux /= 0.) then
        Fbot=mixinglength_flux
        hcond0=-mixinglength_flux/(gamma/(gamma-1.)*gravz/(mpoly0+1.))
        hcond1 = (mpoly1+1.)/(mpoly0+1.)
        hcond2 = (mpoly2+1.)/(mpoly0+1.)
        Kbot=hcond0
        lmultilayer=.true.  ! just to be sure...
        if (lroot) print*, &
            'initialize_energy: hcond0 given by mixinglength_flux=', &
            hcond0, Fbot
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
          if (lroot) print*,'initialize_energy: cs2top,cs2cool=',cs2top,cs2cool
          if (cs2cool /= 0.) then ! cs2cool is the value to go for
            if (lroot) print*,'initialize_energy: now set cs2top=cs2cool'
            cs2top=cs2cool
          else                  ! cs2cool=0, so go for cs2top
            if (lroot) print*,'initialize_energy: now set cs2cool=cs2top'
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
            if (bcz12(iss,1)=='c1') then
              Fbot=-gamma/(gamma-1)*hcond0*gravz/(mpoly0+1)
              if (lroot) print*, &
                      'initialize_energy: Calculated Fbot = ', Fbot
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
            if (bcz12(iss,2)=='c1') then
              Ftop=-gamma/(gamma-1)*hcond0*gravz/(mpoly0+1)
              if (lroot) print*, &
                      'initialize_energy: Calculated Ftop = ',Ftop
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
          !  NOTE: in future we should define chiz=chi(z) or Kz=K(z) here.
          !  calculate hcond and FbotKbot=Fbot/K
          !  (K=hcond is radiative conductivity)
          !
          !  calculate Fbot if it has not been set in run.in
          !
          if (Fbot==impossible) then
            if (bcz12(iss,1)=='c1') then
              Fbot=-gamma/(gamma-1)*hcond0*gravz/(mpoly+1)
              if (lroot) print*, &
                   'initialize_energy: Calculated Fbot = ', Fbot
!
              Kbot=gamma_m1/gamma*(mpoly+1.)*Fbot
              FbotKbot=gamma/gamma_m1/(mpoly+1.)
              if (lroot) print*,'initialize_energy: Calculated Fbot,Kbot=', &
                   Fbot,Kbot
            ! else
            !! Don't need Fbot in this case (?)
            !  Fbot=-gamma/(gamma-1)*hcond0*gravz/(mpoly+1)
            !  if (lroot) print*, &
            !       'initialize_energy: Calculated Fbot = ', Fbot
            endif
          endif
          !
          !  calculate Ftop if it has not been set in run.in
          !
          if (Ftop==impossible) then
            if (bcz12(iss,2)=='c1') then
              Ftop=-gamma/(gamma-1)*hcond0*gravz/(mpoly+1)
              if (lroot) print*, &
                      'initialize_energy: Calculated Ftop = ', Ftop
              Ktop=gamma_m1/gamma*(mpoly+1.)*Ftop
              FtopKtop=gamma/gamma_m1/(mpoly+1.)
              if (lroot) print*,'initialize_energy: Ftop,Ktop=',Ftop,Ktop
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
      call get_cp1(cp1)
      select case (initss(1))
        case ('geo-kws','geo-benchmark','shell_layers')
          if (lroot) then
            print*,'initialize_energy: set boundary temperatures for spherical shell problem'
          endif
!
!  calulate temperature gradient from polytropic index
!
          beta1=cp1*g0/(mpoly+1)*gamma/gamma_m1
!
!  temperatures at shell boundaries
!
          if (initss(1) == 'shell_layers') then
!           lmultilayer=.true.   ! this is the default...
            if (hcond1==impossible) hcond1=(mpoly1+1.)/(mpoly0+1.)
            if (hcond2==impossible) hcond2=(mpoly2+1.)/(mpoly0+1.)
            beta0=-cp1*g0/(mpoly0+1)*gamma/gamma_m1
            beta1=-cp1*g0/(mpoly1+1)*gamma/gamma_m1
            T0=cs20/gamma_m1       ! T0 defined from cs20
            TT_ext=T0
            TT_crit=TT_ext+beta0*(r_bcz-r_ext)
            TT_int=TT_crit+beta1*(r_int-r_bcz)
            cs2top=cs20
            cs2bot=gamma_m1*TT_int
          else
            lmultilayer=.false.  ! to ensure that hcond=cte
            TT_ext=T0            ! T0 defined in start.in for geodynamo
            if (coord_system=='spherical')then
              r_ext=x(l2)
              r_int=x(l1)
            else
            endif
            TT_int=TT_ext*(1.+beta1*(r_ext/r_int-1.))
          endif
          if (lroot) then
            print*,'initialize_energy: g0,mpoly,beta1',g0,mpoly,beta1
            print*,'initialize_energy: TT_int, TT_ext=',TT_int,TT_ext
          endif
!         set up cooling parameters for spherical shell in terms of
!         sound speeds
          call get_soundspeed(log(TT_ext),cs2_ext)
          call get_soundspeed(log(TT_int),cs2_int)
          cs2cool=cs2_ext
!
        case ('star_heat')
          if (hcond1==impossible) hcond1=(mpoly1+1.)/(mpoly0+1.)
          if (hcond2==impossible) hcond2=(mpoly2+1.)/(mpoly0+1.)
          if (lroot) print*,'initialize_energy: set cs2cool=cs20'
          cs2cool=cs20
          cs2_ext=cs20
          if (rcool==0.) rcool=r_ext
          ! only compute the gravity profile
          call star_heat(f,lcompute_grav)
!
        case ('cylind_layers')
          if (bcx12(iss,1)=='c1') then
            Fbot=gamma/gamma_m1*hcond0*g0/(mpoly0+1)
            FbotKbot=gamma/gamma_m1*g0/(mpoly0+1)
          endif
          cs2cool=cs2top
!
        case ('single_polytrope')
          if (cool/=0.) cs2cool=cs0**2
          mpoly=mpoly0  ! needed to compute Fbot when bc=c1 (L383)
!
      endselect
!
!  For global density gradient beta=H/r*dlnrho/dlnr, calculate actual
!  gradient dlnrho/dr = beta/H
!
      if (maxval(abs(beta_glnrho_global))/=0.0) then
        beta_glnrho_scaled=beta_glnrho_global*Omega/cs0
        if (lroot) print*, 'initialize_energy: Global density gradient '// &
            'with beta_glnrho_global=', beta_glnrho_global
      endif
!
!  Turn off pressure gradient term and advection for 0-D runs.
!
      if (nxgrid*nygrid*nzgrid==1) then
        ladvection_entropy=.false.
        print*, 'initialize_energy: 0-D run, turned off advection of energy'
      endif
!
!  Initialize heat conduction.
!
      lheatc_Kprof=.false.
      lheatc_Kconst=.false.
      lheatc_chiconst=.false.
      lheatc_tensordiffusion=.false.
      lheatc_spitzer=.false.
      lheatc_hubeny=.false.
      lheatc_corona=.false.
      lheatc_shock=.false.
      lheatc_hyper3ss=.false.
      lheatc_hyper3ss_polar=.false.
      lheatc_hyper3ss_aniso=.false.
!
      lnothing=.false.
!
!  select which radiative heating we are using
!
      if (lroot) print*,'initialize_energy: nheatc_max,iheatcond=',nheatc_max,iheatcond(1:nheatc_max)
      do i=1,nheatc_max
        select case (iheatcond(i))
        case ('K-const')
          lheatc_Kconst=.true.
          if (lroot) print*, 'heat conduction: K=cte'
        case ('nothing')
          if (lroot .and. (.not. lnothing)) print*,'heat conduction: nothing'
        case default
          if (lroot) then
            write(unit=errormsg,fmt=*)  &
                'No such value iheatcond = ', trim(iheatcond(i))
            call fatal_error('initialize_energy',errormsg)
          endif
        endselect
        lnothing=.true.
      enddo
!
!  A word of warning...
!
      if (all(iheatcond=='nothing') .and. hcond0/=0.0) then
        call warning('initialize_energy', &
            'No heat conduction, but hcond0 /= 0')
      endif
!
! Shared variables
!
      call put_shared_variable('hcond0',hcond0,caller='initialize_energy')
      call put_shared_variable('lviscosity_heat',lviscosity_heat)
      call put_shared_variable('cs2cool',cs2cool)
      call put_shared_variable('mpoly0',mpoly0)
      call put_shared_variable('mpoly1',mpoly2)
      call put_shared_variable('mpoly2',mpoly2)
!
      call keep_compiler_quiet(f)
!
      endsubroutine initialize_energy
!***********************************************************************
    subroutine read_energy_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=entropy_init_pars, IOSTAT=iostat)
!
    endsubroutine read_energy_init_pars
!***********************************************************************
    subroutine write_energy_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=entropy_init_pars)
!
    endsubroutine write_energy_init_pars
!***********************************************************************
    subroutine read_energy_run_pars(iostat)
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=entropy_run_pars, IOSTAT=iostat)
!
    endsubroutine read_energy_run_pars
!***********************************************************************
    subroutine write_energy_run_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=entropy_run_pars)
!
    endsubroutine write_energy_run_pars
!!***********************************************************************
    subroutine init_energy(f)
!
!  initialise energy; called from start.f90
!  07-nov-2001/wolf: coded
!  24-nov-2002/tony: renamed for consistancy (i.e. init_[variable name])
!
      use Sub
      use Gravity
      use General, only: itoa
      use Initcond
      use InitialCondition, only: initial_condition_ss
      use EquationOfState,  only: cs0, rho0, lnrho0, isothermal_entropy, & 
                                  isothermal_lnrho_ss, eoscalc, ilnrho_pp, &
                                  eosperturb
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (nx) :: tmp,pot
      real, dimension (nx) :: pp,lnrho,ss,r_mn
      real, dimension (mx) :: ss_mx
      real :: cs2int,ssint,ztop,ss_ext,pot0,pot_ext
      integer :: j
      logical :: lnothing=.true., save_pretend_lnTT
!
      intent(inout) :: f
!
!  If pretend_lnTT is set then turn it off so that initial conditions are
!  correctly generated in terms of entropy, but then restore it later
!  when we convert the data back again.
!
      save_pretend_lnTT=pretend_lnTT
      pretend_lnTT=.false.
!
      do j=1,ninit
!
        if (initss(j)=='nothing') cycle
!
        lnothing=.false.
        iinit_str=itoa(j)
!
!  select different initial conditions
!
        select case (initss(j))
          case ('anelastic')
          do m=1,my
            do n=1,mz
              f(1:mx,m,n,iss)=0.0
! Base state for isothermal case
              f(1:mx,m,n,iss_b)=-z(n)*gamma_m1*gravz/cs20
            end do
          end do
         case ('polytropic_simple')
           f(:,:,:,iss)=0.0
           print *, 'init_energy: Base entropy defined in density module'
!
!
          case default
!
!  Catch unknown values
!
            write(unit=errormsg,fmt=*) 'No such value for initss(' &
                             //trim(iinit_str)//'): ',trim(initss(j))
            call fatal_error('init_energy',errormsg)
        endselect
!
        if (lroot) print*,'init_energy: initss('//trim(iinit_str)//') = ', &
            trim(initss(j))
!
      enddo
!
      if (lnothing.and.lroot) print*,'init_energy: nothing'
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
!
!  Interface fow user's own initial condition
!
      if (linitial_condition) call initial_condition_ss(f)
!
!  Replace ss by lnTT when pretend_lnTT is true.
!
      if (save_pretend_lnTT) then
        pretend_lnTT=.true.
        do m=1,my; do n=1,mz
          ss_mx=f(:,m,n,iss)
          call eosperturb(f,mx,ss=ss_mx)
        enddo; enddo
      endif
!
    endsubroutine init_energy
!***********************************************************************
    subroutine blob_radeq(ampl,f,i,radius,xblob,yblob,zblob)
!
!  add blob-like perturbation in radiative pressure equilibrium
!
!   6-may-07/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, optional :: xblob,yblob,zblob
      real :: ampl,radius,x01,y01,z01
      integer :: i
!
!  single  blob
!
      if (present(xblob)) then
        x01 = xblob
      else
        x01 = 0.0
      endif
      if (present(yblob)) then
        y01 = yblob
      else
        y01 = 0.0
      endif
      if (present(zblob)) then
        z01 = zblob
      else
        z01 = 0.0
      endif
      if (ampl==0) then
        if (lroot) print*,'ampl=0 in blob_radeq'
      else
        if (lroot.and.ip<14) print*,'blob: variable i,ampl=',i,ampl
        f(:,:,:,i)=f(:,:,:,i)+ampl*(&
           spread(spread(exp(-((x-x01)/radius)**2),2,my),3,mz)&
          *spread(spread(exp(-((y-y01)/radius)**2),1,mx),3,mz)&
          *spread(spread(exp(-((z-z01)/radius)**2),1,mx),2,my))
      endif
!
    endsubroutine blob_radeq
!***********************************************************************
    subroutine polytropic_ss_z( &
         f,mpoly,zint,zbot,zblend,isoth,cs2int,ssint)
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
      real, dimension (mz) :: stp
      real :: tmp,mpoly,zint,zbot,zblend,beta1,cs2int,ssint
      integer :: isoth
!
!  Warning: beta1 is here not dT/dz, but dcs2/dz = (gamma-1)c_p dT/dz
!
      if (headt .and. isoth/=0.0) print*,'ssint=',ssint
      stp = step(z,zblend,widthss)
!
      do n=n1,n2; do m=m1,m2
        if (isoth/=0.0) then ! isothermal layer
          beta1 = 0.0
          tmp = ssint - gamma_m1*gravz*(z(n)-zint)/cs2int
        else
          beta1 = gamma*gravz/(mpoly+1)
          tmp = 1.0 + beta1*(z(n)-zint)/cs2int
          ! Abort if args of log() are negative
          if ( (tmp<=0.0) .and. (z(n)<=zblend) ) then
            call fatal_error('polytropic_ss_z', &
                'Imaginary entropy values -- your z_inf is too low.')
          endif
          tmp = max(tmp,epsi)  ! ensure arg to log is positive
          tmp = ssint + (1-mpoly*gamma_m1)/gamma*log(tmp)
        endif
!
! smoothly blend the old value (above zblend) and the new one (below
! zblend) for the two regions:
!
        f(l1:l2,m,n,iss) = stp(n)*f(l1:l2,m,n,iss)  + (1-stp(n))*tmp
!
      enddo; enddo
!
      if (isoth/=0.0) then
        ssint = -gamma_m1*gravz*(zbot-zint)/cs2int
      else
        ssint = ssint + (1-mpoly*gamma_m1)/gamma &
                      * log(1 + beta1*(zbot-zint)/cs2int)
      endif
      cs2int = cs2int + beta1*(zbot-zint) ! cs2 at layer interface (bottom)
!
    endsubroutine polytropic_ss_z
!***********************************************************************
    subroutine polytropic_ss_disc( &
         f,mpoly,zint,zbot,zblend,isoth,cs2int,ssint)
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
      real, dimension (mz) :: stp
      real :: tmp,mpoly,zint,zbot,zblend,beta1,cs2int,ssint, nu_epicycle2
      integer :: isoth
!
!  Warning: beta1 is here not dT/dz, but dcs2/dz = (gamma-1)c_p dT/dz
!
      stp = step(z,zblend,widthss)
      nu_epicycle2 = nu_epicycle**2
      do n=n1,n2; do m=m1,m2
        if (isoth /= 0) then ! isothermal layer
          beta1 = 0.
          tmp = ssint - gamma_m1*gravz*nu_epicycle2*(z(n)**2-zint**2)/cs2int/2.
        else
          beta1 = gamma*gravz*nu_epicycle2/(mpoly+1)
          tmp = 1 + beta1*(z(n)**2-zint**2)/cs2int/2.
          ! Abort if args of log() are negative
          if ( (tmp<=0.0) .and. (z(n)<=zblend) ) then
            call fatal_error('polytropic_ss_disc', &
                'Imaginary entropy values -- your z_inf is too low.')
          endif
          tmp = max(tmp,epsi)  ! ensure arg to log is positive
          tmp = ssint + (1-mpoly*gamma_m1)/gamma*log(tmp)
        endif
!
! smoothly blend the old value (above zblend) and the new one (below
! zblend) for the two regions:
!
        f(l1:l2,m,n,iss) = stp(n)*f(l1:l2,m,n,iss)  + (1-stp(n))*tmp
      enddo; enddo
!
      if (isoth/=0.0) then
        ssint = -gamma_m1*gravz*nu_epicycle2*(zbot**2-zint**2)/cs2int/2.
      else
        ssint = ssint + (1-mpoly*gamma_m1)/gamma &
                      * log(1 + beta1*(zbot**2-zint**2)/cs2int/2.)
      endif
!
      cs2int = cs2int + beta1*(zbot**2-zint**2)/2.
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
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, intent(in) :: lnrho_bot,ss_const
      real :: cs2_,lnrho,lnrho_m
!
      if (.not. lgravz) then
        call fatal_error("hydrostatic_isentropic","Currently only works for vertical gravity field")
      endif
!
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
      !
      ! Entropy is simply constant
      !
      f(:,:,:,iss)=ss_const
!
    endsubroutine hydrostatic_isentropic
!***********************************************************************
    subroutine mixinglength(mixinglength_flux,f)
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
!  17-Nov-05/dintrans: updated using strat_MLT
!
      use Gravity, only: gravz, z1
      use General, only: safe_character_assign
      use EquationOfState, only: gamma, gamma_m1, rho0, lnrho0, cs0, &
                                 cs20, cs2top, eoscalc, ilnrho_lnTT
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nzgrid) :: tempm,lnrhom
      real :: zm,ztop,mixinglength_flux,lnrho,ss,lnTT
      real :: zbot,rbot,rt_old,rt_new,rb_old,rb_new,crit, &
              rhotop,rhobot
      integer :: iter
      character (len=120) :: wfile
!
      if (headtt) print*,'init_energy : mixinglength stratification'
      if (.not.lgravz) then
        call fatal_error("mixinglength","works only for vertical gravity")
      endif
!
!  do the calculation on all processors, and then put the relevant piece
!  into the f array.
!  choose value zbot where rhobot should be applied and give two first
!  estimates for rhotop
!
      ztop=xyz0(3)+Lxyz(3)
      zbot=z1
      rbot=1.
      rt_old=.1*rbot
      rt_new=.12*rbot
!
!  need to iterate for rhobot=1.
!  produce first estimate
!
      rhotop=rt_old
      cs2top=cs20  ! just to be sure...
      call strat_MLT (rhotop,mixinglength_flux,lnrhom,tempm,rhobot)
      rb_old=rhobot
!
!  next estimate
!
      rhotop=rt_new
      call strat_MLT (rhotop,mixinglength_flux,lnrhom,tempm,rhobot)
      rb_new=rhobot
!
      do iter=1,10
!
!  new estimate
!
        rhotop=rt_old+(rt_new-rt_old)/(rb_new-rb_old)*(rbot-rb_old)
!
!  check convergence
!
        crit=abs(rhotop/rt_new-1.)
        if (crit<=1e-4) exit
!
        call strat_MLT (rhotop,mixinglength_flux,lnrhom,tempm,rhobot)
!
!  update new estimates
!
        rt_old=rt_new
        rb_old=rb_new
        rt_new=rhotop
        rb_new=rhobot
      enddo
      if (ipz==0) print*,'- iteration completed: rhotop,crit=',rhotop,crit
!
! redefine rho0 and lnrho0 as we don't have rho0=1 at the top
! (important for eoscalc!)
! put density and entropy into f-array
! write the initial stratification in data/proc*/stratMLT.dat
!
      rho0=rhotop
      lnrho0=log(rhotop)
      print*,'new rho0 and lnrho0=',rho0,lnrho0
!
      call safe_character_assign(wfile,trim(directory)//'/stratMLT.dat')
      open(11+ipz,file=wfile,status='unknown')
      do n=1,nz
        iz=n+ipz*nz
        zm=xyz0(3)+(iz-1)*dz
        lnrho=lnrhom(nzgrid-iz+1)
        lnTT=log(tempm(nzgrid-iz+1))
        call eoscalc(ilnrho_lnTT,lnrho,lnTT,ss=ss)
        f(:,:,n+nghost,ilnrho)=lnrho
        f(:,:,n+nghost,iss)=ss
        write(11+ipz,'(4(2x,1pe12.5))') zm,exp(lnrho),ss, &
              gamma_m1*tempm(nzgrid-iz+1)
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
!  21-aug-08/dhruba: added spherical coordinates
!
      use Gravity, only: g0
      use EquationOfState, only: eoscalc, ilnrho_lnTT, mpoly, get_cp1
      use Mpicomm, only:stop_it
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: lnrho,lnTT,TT,ss,pert_TT,r_mn
      real :: beta1,cp1
!
!  beta1 is the temperature gradient
!  1/beta = (g/cp) 1./[(1-1/gamma)*(m+1)]
!
      call get_cp1(cp1)
      beta1=cp1*g0/(mpoly+1)*gamma/gamma_m1
!
!  set intial condition
!
      if (lspherical_coords)then
        do imn=1,ny*nz
          n=nn(imn)
          m=mm(imn)
          call shell_ss_perturb(pert_TT)
!
          TT(1)=TT_int
          TT(2:nx-1)=TT_ext*(1.+beta1*(x(l2)/x(l1+1:l2-1)-1))+pert_TT(2:nx-1)
          TT(nx)=TT_ext
!
          lnrho=f(l1:l2,m,n,ilnrho)
          lnTT=log(TT)
          call eoscalc(ilnrho_lnTT,lnrho,lnTT,ss=ss)
          f(l1:l2,m,n,iss)=ss
!
        enddo
      elseif (lcylindrical_coords) then
        call stop_it('shell_ss:not valid in cylindrical coordinates')
      else
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
      endif
!
    endsubroutine shell_ss
!***********************************************************************
    subroutine shell_ss_perturb(pert_TT)
!
!  Compute perturbation to initial temperature profile
!
!  22-june-04/dave -- coded
!  21-aug-08/dhruba : added spherical coords
!
      real, dimension (nx), intent(out) :: pert_TT
      real, dimension (nx) :: xr,cos_4phi,sin_theta4,r_mn,rcyl_mn,phi_mn
      real, parameter :: ampl0=.885065
!
      select case (initss(1))
!
        case ('geo-kws')
          pert_TT=0.
!
        case ('geo-benchmark')
          if (lspherical_coords)then
            xr=x(l1:l2)
            sin_theta4=sin(y(m))**4
            pert_TT=ampl0*ampl_TT*(1-3*xr**2+3*xr**4-xr**6)*&
                          (sin(y(m))**4)*cos(4*z(n))
          else
            r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
            rcyl_mn=sqrt(x(l1:l2)**2+y(m)**2)
            phi_mn=atan2(spread(y(m),1,nx),x(l1:l2))
            xr=2*r_mn-r_int-r_ext              ! radial part of perturbation
            cos_4phi=cos(4*phi_mn)             ! azimuthal part
            sin_theta4=(rcyl_mn/r_mn)**4       ! meridional part
            pert_TT=ampl0*ampl_TT*(1-3*xr**2+3*xr**4-xr**6)*sin_theta4*cos_4phi
          endif
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
      use Gravity, only: gravz, zinfty
      use EquationOfState, only: eoscalc, ilnrho_lnTT, mpoly, get_cp1
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension (nx) :: lnrho,lnTT,TT,ss,z_mn
      real :: beta1,cp1
!
!  beta1 is the temperature gradient
!  1/beta = (g/cp) 1./[(1-1/gamma)*(m+1)]
!  Also set dcs2top_ini in case one uses radiative b.c.
!  NOTE: cs2top_ini=cs20 would be wrong if zinfty is not correct in gravity
!
      call get_cp1(cp1)
      beta1=cp1*gamma/gamma_m1*gravz/(mpoly+1)
      dcs2top_ini=gamma_m1*gravz
      cs2top_ini=cs20
!
!  set initial condition (first in terms of TT, and then in terms of ss)
!
      do m=m1,m2
      do n=n1,n2
        z_mn = spread(z(n),1,nx)
        TT = beta1*(z_mn-zinfty)
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
      real :: kpc,fmpi1
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
!
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
!
!  define entropy via pressure, assuming fixed T for each component
        if (lentropy) then
!  thermal pressure (eq 15)
          pp=real(k_B*unit_length**3 * &
             (1.09*n_c*T_c + 1.09*n_w*T_w + 2.09*n_i*T_i + 2.27*n_h*T_h))
!
          call eosperturb(f,nx,pp=pp)
          ss=f(l1:l2,m,n,ilnrho)
!
          fmpi1=cs2bot
          call mpibcast_real(fmpi1,0)
          cs2bot=fmpi1
          fmpi1=cs2top
          call mpibcast_real(fmpi1,ncpus-1)
          cs2top=fmpi1
!
        endif
       enddo
      enddo
!
      if (lroot) print*, 'ferriere: cs2bot=',cs2bot, ' cs2top=',cs2top
!
    endsubroutine ferriere
!***********************************************************************
    subroutine galactic_hs(f,rho0hs,cs0hs,H0hs)
!
!   Density and isothermal entropy profile in hydrostatic equilibrium
!   with the galactic-hs-gravity profile set in gravity_simple.
!   Parameters cs0hs and H0hs need to be set consistently
!   both in grav_init_pars and in entropy_init_pars to obtain hydrostatic
!   equilibrium.
!
      use Mpicomm, only: mpibcast_real
      use EquationOfState, only: eosperturb
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(nx) :: rho,pp,lnrho,ss
      real :: rho0hs,cs0hs,H0hs,fmpi1
!
      if (lroot) print*, &
         'Galactic-hs: hydrostatic equilibrium density and entropy profiles'
!
      do n=n1,n2
      do m=m1,m2
        rho=rho0hs*exp(1 - sqrt(1 + (z(n)/H0hs)**2))
        lnrho=log(rho)
        f(l1:l2,m,n,ilnrho)=lnrho
        if (lentropy) then
!  Isothermal
          pp=rho*cs0hs**2
          call eosperturb(f,nx,pp=pp)
          ss=f(l1:l2,m,n,ilnrho)
          fmpi1=cs2bot
          call mpibcast_real(fmpi1,0)
          cs2bot=fmpi1
          fmpi1=cs2top
          call mpibcast_real(fmpi1,ncpus-1)
          cs2top=fmpi1
!
         endif
       enddo
     enddo
!
      if (lroot) print*, 'Galactic-hs: cs2bot=',cs2bot, ' cs2top=',cs2top
!
    endsubroutine galactic_hs
!***********************************************************************
    subroutine shock2d(f)
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
!
      real, dimension (4) :: rpp,rpr,rpu,rpv
      integer :: l
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
        do n=n1,n2; do m=m1,m2; do l=l1,l2
          if ( x(l)>=0.0 .and. y(m)>=0.0 ) then
            f(l,m,n,ilnrho)=log(rpr(1))
            f(l,m,n,iss)=log(gamma*rpp(1))/gamma-f(l,m,n,ilnrho)
            f(l,m,n,iux)=rpu(1)
            f(l,m,n,iuy)=rpv(1)
          endif
          if ( x(l)<0.0 .and. y(m)>=0.0 ) then
            f(l,m,n,ilnrho)=log(rpr(2))
            f(l,m,n,iss)=log(gamma*rpp(2))/gamma-f(l,m,n,ilnrho)
            f(l,m,n,iux)=rpu(2)
            f(l,m,n,iuy)=rpv(2)
          endif
          if ( x(l)<0.0 .and. y(m)<0.0 ) then
            f(l,m,n,ilnrho)=log(rpr(3))
            f(l,m,n,iss)=log(gamma*rpp(3))/gamma-f(l,m,n,ilnrho)
            f(l,m,n,iux)=rpu(3)
            f(l,m,n,iuy)=rpv(3)
          endif
          if ( x(l)>=0.0 .and. y(m)<0.0 ) then
            f(l,m,n,ilnrho)=log(rpr(4))
            f(l,m,n,iss)=log(gamma*rpp(4))/gamma-f(l,m,n,ilnrho)
            f(l,m,n,iux)=rpu(4)
            f(l,m,n,iuy)=rpv(4)
          endif
        enddo; enddo; enddo
!
    endsubroutine shock2d
!***********************************************************************
    subroutine pencil_criteria_energy()
!
!  All pencils that the Energy module depends on are specified here.
!
!  20-11-04/anders: coded
!
      if (lheatc_Kconst .or. lheatc_chiconst .or. lheatc_Kprof .or. &
          tau_cor>0) lpenc_requested(i_cp1)=.true.
      if (ladvection_entropy) then
        lpenc_requested(i_ugss)=.true.
        lpenc_requested(i_ugss_b)=.true.
      end if
      if (lviscosity.and.lviscosity_heat) then
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_visc_heat)=.true.
        if (pretend_lnTT) lpenc_requested(i_cv1)=.true.
      endif
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
      if (lgravr) then
! spherical case (cylindrical case also included)
        if (lcylindrical_coords) then
          lpenc_requested(i_rcyl_mn)=.true.
        else
          lpenc_requested(i_r_mn)=.true.
        endif
      endif
      if (lheatc_Kconst) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_glnTTb)=.true.
        lpenc_requested(i_del2ss)=.true.
      endif
!
      lpenc_diagnos2d(i_ss)=.true.
!
      if (idiag_dtchi/=0) lpenc_diagnos(i_rho1)=.true.
      if (idiag_ethdivum/=0) lpenc_diagnos(i_divu)=.true.
      if (idiag_csm/=0) lpenc_diagnos(i_cs2)=.true.
      if (idiag_eem/=0) lpenc_diagnos(i_ee)=.true.
      if (idiag_ppm/=0) lpenc_diagnos(i_pp)=.true.
      if (idiag_pdivum/=0) then
        lpenc_diagnos(i_pp)=.true.
        lpenc_diagnos(i_divu)=.true.
      endif
      if (idiag_ssm/=0 .or. idiag_ss2m/=0 .or. idiag_ssmz/=0 .or. &
          idiag_ssmy/=0.or.idiag_ssmx/=0.or.idiag_ssmr/=0) &
           lpenc_diagnos(i_ss)=.true.
      lpenc_diagnos(i_rho)=.true.
      lpenc_diagnos(i_ee)=.true.
      if (idiag_ethm/=0 .or. idiag_ethtot/=0 .or. idiag_ethdivum/=0 ) then
          lpenc_diagnos(i_rho)=.true.
          lpenc_diagnos(i_ee)=.true.
      endif
      if (idiag_fconvm/=0 .or. idiag_fconvz/=0 .or. idiag_fturbz/=0) then
          lpenc_diagnos(i_rho)=.true.
          lpenc_diagnos(i_TT)=.true.  !(to be replaced by enthalpy)
      endif
      if (idiag_TTm/=0 .or. idiag_TTmx/=0 .or. idiag_TTmy/=0 .or. &
          idiag_TTmz/=0 .or. idiag_TTmr/=0 .or. idiag_TTmax/=0 .or. &
          idiag_TTmin/=0 .or. idiag_uxTTmz/=0 .or.idiag_uyTTmz/=0 .or. &
          idiag_uzTTmz/=0) &
          lpenc_diagnos(i_TT)=.true.
      if (idiag_TTmxy/=0 .or. idiag_TTmxz/=0) lpenc_diagnos2d(i_TT)=.true.
      if (idiag_yHm/=0 .or. idiag_yHmax/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_dtc/=0) lpenc_diagnos(i_cs2)=.true.
      if (idiag_TTp/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_cs2)=.true.
        lpenc_diagnos(i_rcyl_mn)=.true.
      endif
      if (idiag_cs2mphi/=0) lpenc_diagnos(i_cs2)=.true.
!
    endsubroutine pencil_criteria_energy
!***********************************************************************
    subroutine pencil_interdep_energy(lpencil_in)
!
!  Interdependency among pencils from the Energy module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_ugss)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_gss)=.true.
      endif
!
      if (lpencil_in(i_ugss_b)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_gss_b)=.true.
      endif
!
      if (lpencil_in(i_glnTTb)) then
        lpencil_in(i_grho)=.true.
        lpencil_in(i_gss_b)=.true.
      endif
!
      if (lpencil_in(i_uglnTT)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_glnTT)=.true.
      endif
      if (lpencil_in(i_Ma2)) then
        lpencil_in(i_u2)=.true.
        lpencil_in(i_cs2)=.true.
      endif
      if (lpencil_in(i_fpres)) then
        lpencil_in(i_cs2)=.true.
        lpencil_in(i_glnrho)=.true.
        if (leos_idealgas) then
          lpencil_in(i_glnTT)=.true.
        else
          lpencil_in(i_cp1tilde)=.true.
          lpencil_in(i_gss)=.true.
        endif
      endif
!
    endsubroutine pencil_interdep_energy
!***********************************************************************
    subroutine calc_pencils_energy(f,p)
!
!  Calculate all energy pencils except the pressure gradient
!
!  08-dec-2009/piyali:adapted
!
      use EquationOfState, only: gamma,gamma_m1,cs20,lnrho0
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
! ss
!
      if (lpencil(i_ss)) p%ss = f(l1:l2,m,n,iss)
      if (lpencil(i_gss)) call grad(f,iss,p%gss)
      if (lpencil(i_gss_b)) call grad(f,iss_b,p%gss_b)
      if (lpencil(i_del2ss)) call del2(f,iss,p%del2ss)
!
! Ma2
!
      if (lpencil(i_Ma2)) p%Ma2=p%u2/p%cs2
!
! ugss & ugss_b
!
      if (lpencil(i_ugss)) &
          call u_dot_grad(f,iss,p%gss,p%uu,p%ugss,UPWIND=lupw_ss)
      if (lpencil(i_ugss_b)) &
          call u_dot_grad(f,iss_b,p%gss_b,p%uu,p%ugss_b,UPWIND=lupw_ss)
      if (lpencil(i_glnTTb)) &
          p%glnTTb(:,1)=gamma*p%gss_b(:,1)*cp1+gamma_m1*p%grho(:,1)*p%rho1
          p%glnTTb(:,2)=gamma*p%gss_b(:,2)*cp1+gamma_m1*p%grho(:,2)*p%rho1
          p%glnTTb(:,3)=gamma*p%gss_b(:,3)*cp1+gamma_m1*p%grho(:,3)*p%rho1
!
!
    endsubroutine calc_pencils_energy
!**********************************************************************
    subroutine denergy_dt(f,df,p)
!
!  Calculate right hand side of energy equation,
!  ds/dt = -u.grads + [H-C + div(K*gradT) + mu0*eta*J^2 + ...]
!
!   08-dec-09/piyali: adapted from entropy.f90
      use Diagnostics
      use EquationOfState, only: gamma1, cs0
      use Special, only: special_calc_energy
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: Hmax
      real :: ztop,xi,profile_cor,uT,fradz,TTtop
      integer :: j,ju
!
      intent(inout)  :: f,p
      intent(out) :: df
!
      Hmax = 0.
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'denergy_dt: SOLVE denergy_dt'
      if (headtt) call identify_bcs('ss',iss)
      if (headtt) call identify_bcs('ss_b',iss_b)
      if (headtt) print*,'denergy_dt: lnTT,cs2,cp1=', p%lnTT(1), p%cs2(1), p%cp1(1)
!
!  ``cs2/dx^2'' for timestep
!
      if (lfirst.and.ldt) advec_cs2=p%cs2*dxyz_2
      if (headtt.or.ldebug) print*,'denergy_dt: max(advec_cs2) =',maxval(advec_cs2)
!
!
!  Advection of entropy.
!  If pretend_lnTT=.true., we pretend that ss is actually lnTT
!  Otherwise, in the regular case with entropy, s is the dimensional
!  specific entropy, i.e. it is not divided by cp.
!  NOTE: in the entropy module is it lnTT that is advanced, so
!  there are additional cv1 terms on the right hand side.
!
      if (ladvection_entropy) &
           df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) - p%ugss - p%ugss_b
!
!  Calculate viscous contribution to entropy
!
      if (lviscosity .and. lviscosity_heat) call calc_viscous_heat(df,p,Hmax)
!
!  Entry possibility for "personal" entries.
!  In that case you'd need to provide your own "special" routine.
!
      if (lspecial) call special_calc_energy(f,df,p)
!
!  Thermal conduction delegated to different subroutines.
!
      diffus_chi=0.; diffus_chi3=0.

      if (lheatc_Kconst)   call calc_heatcond_constK(df,p)
!
!  Interstellar radiative cooling and UV heating.
!
      if (linterstellar) call calc_heat_cool_interstellar(f,df,p,Hmax)
!
!  Explicit heating/cooling terms.
!
      if ((luminosity/=0.0) .or. (cool/=0.0) .or. &
          (tau_cor/=0.0) .or. (tauheat_buffer/=0.0) .or. &
          (heat_uniform/=0.0) .or. (tau_cool/=0.0) .or. &
          (cool_ext/=0.0 .and. cool_int/=0.0) .or. lturbulent_heat .or. &
          (tau_cool2 /=0)) &
          call calc_heat_cool(df,p,Hmax)
      if (tdown/=0.0) call newton_cool(df,p)
      if (cool_RTV/=0.0) call calc_heat_cool_RTV(df,p)
!
      if (lfirst.and.ldt) then
        dt1_max=max(dt1_max,Hmax/ee/cdts)
        maxdiffus=max(maxdiffus,diffus_chi)
        maxdiffus3=max(maxdiffus3,diffus_chi3)
      endif
!
!  Possibility of entropy relaxation in exterior region.
!
      if (tau_ss_exterior/=0.) call calc_tau_ss_exterior(df,p)
!
!  Apply border profile
!
      if (lborder_profiles) call set_border_entropy(f,df,p)
!
!  Phi-averages
!
      if (l2davgfirst) then
        if (idiag_ssmphi/=0)  call phisum_mn_name_rz(p%ss,idiag_ssmphi)
        if (idiag_cs2mphi/=0) call phisum_mn_name_rz(p%cs2,idiag_cs2mphi)
      endif
!
!  Enforce maximum heating rate timestep constraint
!
!      if (lfirst.and.ldt) dt1_max=max(dt1_max,Hmax/ee/cdts)
!
!  Calculate entropy related diagnostics.
!
      if (ldiagnos) then
        !uT=unit_temperature !(define shorthand to avoid long lines below)
        uT=1. !(AB: for the time being; to keep compatible with auto-test
        if (idiag_TTmax/=0)  call max_mn_name(p%TT*uT,idiag_TTmax)
        if (idiag_TTmin/=0)  call max_mn_name(-p%TT*uT,idiag_TTmin,lneg=.true.)
        if (idiag_TTm/=0)    call sum_mn_name(p%TT*uT,idiag_TTm)
        if (idiag_pdivum/=0) call sum_mn_name(p%pp*p%divu,idiag_pdivum)
        if (idiag_yHmax/=0)  call max_mn_name(p%yH,idiag_yHmax)
        if (idiag_yHm/=0)    call sum_mn_name(p%yH,idiag_yHm)
        if (idiag_dtc/=0) &
            call max_mn_name(sqrt(advec_cs2)/cdt,idiag_dtc,l_dt=.true.)
        if (idiag_ethm/=0)    call sum_mn_name(p%rho*p%ee,idiag_ethm)
        if (idiag_ethtot/=0) call integrate_mn_name(p%rho*p%ee,idiag_ethtot)
        if (idiag_ethdivum/=0) &
            call sum_mn_name(p%rho*p%ee*p%divu,idiag_ethdivum)
        if (idiag_ssm/=0) call sum_mn_name(p%ss,idiag_ssm)
        if (idiag_ss2m/=0) call sum_mn_name(p%ss**2,idiag_ss2m)
        if (idiag_eem/=0) call sum_mn_name(p%ee,idiag_eem)
        if (idiag_ppm/=0) call sum_mn_name(p%pp,idiag_ppm)
        if (idiag_csm/=0) call sum_mn_name(p%cs2,idiag_csm,lsqrt=.true.)
        if (idiag_ugradpm/=0) &
            call sum_mn_name(p%cs2*(p%uglnrho+p%ugss),idiag_ugradpm)
        if (idiag_fconvm/=0) &
            call sum_mn_name(p%rho*p%uu(:,3)*p%TT,idiag_fconvm)
!
!  radiative heat flux at the bottom (assume here that hcond=hcond0=const)
!
        if (idiag_fradbot/=0) then
          if (ipz==0       .and.n==n1) then
            fradz=sum(-hcond0*p%TT*p%glnTT(:,3)*dsurfxy)
          else
            fradz=0.
          endif
          call surf_mn_name(fradz,idiag_fradbot)
        endif
!
!  radiative heat flux at the top (assume here that hcond=hcond0=const)
!
        if (idiag_fradtop/=0) then
          if (llast_proc_z.and.n==n2) then
            fradz=sum(-hcond0*p%TT*p%glnTT(:,3)*dsurfxy)
          else
            fradz=0.
          endif
          call surf_mn_name(fradz,idiag_fradtop)
        endif
!
!  mean temperature at the top
!
        if (idiag_TTtop/=0) then
          if (llast_proc_z.and.n==n2) then
            TTtop=sum(p%TT*dsurfxy)
          else
            TTtop=0.
          endif
          call surf_mn_name(TTtop,idiag_TTtop)
        endif
!
!  calculate integrated temperature in in limited radial range
!
        if (idiag_TTp/=0) call sum_lim_mn_name(p%rho*p%cs2*gamma1,idiag_TTp,p)
      endif
!
!  1-D averages.
!
      if (l1davgfirst) then
        if (idiag_fradz/=0) call xysum_mn_name_z(-hcond0*p%TT*p%glnTT(:,3),idiag_fradz)
        if (idiag_fconvz/=0) &
            call xysum_mn_name_z(p%rho*p%uu(:,3)*p%TT,idiag_fconvz)
        if (idiag_ssmz/=0)  call xysum_mn_name_z(p%ss,idiag_ssmz)
        if (idiag_ssmy/=0)  call xzsum_mn_name_y(p%ss,idiag_ssmy)
        if (idiag_ssmx/=0)  call yzsum_mn_name_x(p%ss,idiag_ssmx)
        if (idiag_TTmx/=0)  call yzsum_mn_name_x(p%TT,idiag_TTmx)
        if (idiag_TTmy/=0)  call xzsum_mn_name_y(p%TT,idiag_TTmy)
        if (idiag_TTmz/=0)  call xysum_mn_name_z(p%TT,idiag_TTmz)
        if (idiag_ssmr/=0)  call phizsum_mn_name_r(p%ss,idiag_ssmr)
        if (idiag_TTmr/=0)  call phizsum_mn_name_r(p%TT,idiag_TTmr)
        if (idiag_uxTTmz/=0) &
            call xysum_mn_name_z(p%uu(:,1)*p%TT,idiag_uxTTmz)
        if (idiag_uyTTmz/=0) &
            call xysum_mn_name_z(p%uu(:,2)*p%TT,idiag_uyTTmz)
        if (idiag_uzTTmz/=0) &
            call xysum_mn_name_z(p%uu(:,3)*p%TT,idiag_uzTTmz)
      endif
!
!  2-D averages.
!
      if (l2davgfirst) then
        if (idiag_TTmxy/=0) call zsum_mn_name_xy(p%TT,idiag_TTmxy)
        if (idiag_TTmxz/=0) call ysum_mn_name_xz(p%TT,idiag_TTmxz)
        if (idiag_ssmxy/=0) call zsum_mn_name_xy(p%ss,idiag_ssmxy)
        if (idiag_ssmxz/=0) call ysum_mn_name_xz(p%ss,idiag_ssmxz)
        if (idiag_uxTTmxy/=0) call zsum_mn_name_xy(p%uu(:,1)*p%TT,idiag_uxTTmxy)
      endif
!
    endsubroutine denergy_dt
!***********************************************************************
    subroutine energy_after_boundary(f)
!
!  Dummy routine.
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine energy_after_boundary
!***********************************************************************
    subroutine set_border_entropy(f,df,p)
!
!  Calculates the driving term for the border profile
!  of the ss variable.
!
!  28-jul-06/wlad: coded
!
      use BorderProfiles, only: border_driving
      use EquationOfState, only: cs20,get_ptlaw,get_cp1,lnrho0
      use Sub, only: power_law
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(nx) :: f_target,cs2
      real :: ptlaw,cp1
!
      select case (borderss)
!
      case ('zero','0')
         f_target=0.
      case ('constant')
         f_target=ss_const
      case ('power-law')
        call get_ptlaw(ptlaw)
        call get_cp1(cp1)
        call power_law(cs20,p%rcyl_mn,ptlaw,cs2,r_ref)
        f_target=1./(gamma*cp1)*(log(cs2/cs20)-gamma_m1*lnrho0)
         !f_target= gamma1*log(cs2_0) !- gamma_m1*gamma1*lnrho
      case ('nothing')
         if (lroot.and.ip<=5) &
              print*,"set_border_entropy: borderss='nothing'"
      case default
         write(unit=errormsg,fmt=*) &
              'set_border_entropy: No such value for borderss: ', &
              trim(borderss)
         call fatal_error('set_border_entropy',errormsg)
      endselect
!
      if (borderss /= 'nothing') then
        call border_driving(f,df,p,f_target,iss)
      endif
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
      use Diagnostics
      use Gravity
      use Sub
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
!
!  Provisional expression for magnetic chi_t quenching;
!  (Derivatives of B are still missing.
!
          if (chiB==0.) then
            thdiff=thdiff+chi_t*(p%del2ss+g2)
          else
            thdiff=thdiff+chi_t*(p%del2ss+g2)/(1.+chiB*p%b2)
          endif
        endif
      endif
!
!  add heat conduction to energy equation
!
      df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)+thdiff
      if (headtt) print*,'calc_heatcond_constchi: added thdiff'
!
!  check maximum diffusion from thermal diffusion
!  With heat conduction, the second-order term for entropy is
!  gamma*chi*del2ss
!
      if (lfirst.and.ldt) then
        if (leos_idealgas) then
          diffus_chi=diffus_chi+(gamma*chi+chi_t)*dxyz_2
        else
          diffus_chi=diffus_chi+(chi+chi_t)*dxyz_2
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
!  add heat conduction to energy equation
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
      if (headtt) print*,'calc_heatcond_hyper3: added thdiff'
!
!  check maximum diffusion from thermal diffusion
!
      if (lfirst.and.ldt) diffus_chi3=diffus_chi3+chi_hyper3*dxyz_6
!
    endsubroutine calc_heatcond_hyper3
!***********************************************************************
    subroutine calc_heatcond_hyper3_aniso(f,df)
!
!  Naive anisotropic hyperdiffusivity of entropy.
!
!  11-may-09/wlad: coded
!
      use Sub, only: del6fj
!
      real, dimension (mx,my,mz,mfarray),intent(in) :: f
      real, dimension (mx,my,mz,mvar),intent(out) :: df
!
      real, dimension (nx) :: thdiff,tmp
!
!  check that chi_hyper3_aniso is ok
!
      if (headtt) print*, &
           'calc_heatcond_hyper3_aniso: chi_hyper3_aniso=', &
           chi_hyper3_aniso
!
!  Heat conduction
!
      call del6fj(f,chi_hyper3_aniso,iss,tmp)
      thdiff = tmp
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
!
      if (lfirst.and.ldt) diffus_chi3=diffus_chi3+ &
                                      (chi_hyper3_aniso(1)*dline_1(:,1)**6 + &
                                       chi_hyper3_aniso(2)*dline_1(:,2)**6 + &
                                       chi_hyper3_aniso(3)*dline_1(:,3)**6)
!
    endsubroutine calc_heatcond_hyper3_aniso
!***********************************************************************
    subroutine calc_heatcond_hyper3_polar(f,df)
!
!  Naive hyperdiffusivity of entropy in polar coordinates
!
!  03-aug-08/wlad: coded
!
      use Deriv, only: der6
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: tmp
      integer :: j
!
      real, dimension (nx) :: thdiff
!
      intent(in)  :: f
      intent(out) :: df
!
      if (headtt) print*, 'calc_heatcond_hyper3: chi_hyper3=', chi_hyper3
!
      do j=1,3
        call der6(f,iss,tmp,j,IGNOREDX=.true.)
        thdiff = thdiff + chi_hyper3*pi4_1*tmp*dline_1(:,j)**2
      enddo
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
!
      if (headtt) print*,'calc_heatcond_hyper3: added thdiff'
!
      if (lfirst.and.ldt) &
           diffus_chi3=diffus_chi3+diffus_chi3*pi4_1*dxyz_2
!
    endsubroutine calc_heatcond_hyper3_polar
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
      use Diagnostics
      use EquationOfState, only: gamma
      use Gravity
      use Sub
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
      if (headtt) print*,'calc_heatcond_shock: chi_t,chi_shock=',chi_t,chi_shock
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
        thdiff=chi_shock*(p%shock*(p%del2lnTT+g2)+gshockglnTT)
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
          diffus_chi=diffus_chi+(chi_t+gamma*chi_shock*p%shock)*dxyz_2
        else
          diffus_chi=diffus_chi+(chi_t+chi_shock*p%shock)*dxyz_2
        endif
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_shock
!***********************************************************************
    subroutine calc_heatcond_constK(df,p)
!
!  heat conduction
!
!   8-jul-02/axel: adapted from Wolfgang's more complex version
!  30-mar-06/ngrs: simplified calculations using p%glnTT and p%del2lnTT
!
      use Diagnostics
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
      chix = p%rho1*hcond0*p%cp1
!
      call dot(p%gss,p%glnTTb,g2)
      thdiff = chix*(p%del2ss + g2)
!
!  add heat conduction to entropy equation
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
      if (headtt) print*,'calc_heatcond_constK: added thdiff'
!
!  check maximum diffusion from thermal diffusion
!  With heat conduction, the second-order term for entropy is
!  gamma*chix*del2ss
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chix*dxyz_2
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_constK
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
      use EquationOfState, only: gamma,gamma_m1
      use Sub
      use IO, only: output_pencil
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: gvKpara,gvKperp,tmpv,tmpv2
      real, dimension (nx) :: bb2,thdiff,b1
      real, dimension (nx) :: tmps,quenchfactor,vKpara,vKperp
!
      integer ::i,j
      type (pencil_case) :: p
!
      intent(in) :: p
      intent(out) :: df
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
      if (maxval(abs(vKpara+vKperp)) > tini) then
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
      call multsv_mn(2*b1*b1,tmpv,tmpv2)
      tmpv=2.*p%glnrho+0.5*p%glnTT-tmpv2
      call multsv_mn(vKperp,tmpv,gvKperp)
      gvKperp=gvKperp*spread(quenchfactor,2,3)
!
!     Calculate diffusion term
!
      call tensor_diffusion_coef(p%glnTT,p%hlnTT,p%bij,p%bb,vKperp,vKpara,thdiff,GVKPERP=gvKperp,GVKPARA=gvKpara)
!
      if (lfirst .and. ip == 13) then
         call output_pencil('spitzer.dat',thdiff,1)
         call output_pencil('viscous.dat',p%visc_heat*exp(p%lnrho),1)
      endif
!
      thdiff = thdiff*exp(-p%lnrho-p%lnTT)
!
      df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) + thdiff
!
      if (lfirst.and.ldt) then
!
         dt1_max=max(dt1_max,maxval(abs(thdiff)*gamma)/(cdts))
         diffus_chi=diffus_chi+gamma*Kgpara*exp(2.5*p%lnTT-p%lnrho)/p%cp*dxyz_2
      endif
!
    endsubroutine calc_heatcond_spitzer
!***********************************************************************
    subroutine calc_heatcond_tensor(df,p)
!
!  Calculates heat conduction parallel and perpendicular (isotropic)
!  to magnetic field lines
!
!
!  24-aug-09/bing: moved from denergy_dt to here
!
      use Diagnostics
      use Sub, only: tensor_diffusion_coef,dot,dot2
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: cosbgT,gT2,b2,rhs
      real, dimension (nx) :: vKpara,vKperp
!
      type (pencil_case) :: p
!
      vKpara(:) = Kgpara
      vKperp(:) = Kgperp
!
      call tensor_diffusion_coef(p%glnTT,p%hlnTT,p%bij,p%bb,vKperp,vKpara,rhs,llog=.true.)
!
      df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)+rhs*p%rho1
!
      call dot(p%bb,p%glnTT,cosbgT)
      call dot2(p%glnTT,gT2)
      call dot2(p%bb,b2)
!
      where ((gT2<=tini).or.(b2<=tini))
         cosbgT=0.
      elsewhere
         cosbgT=cosbgT/sqrt(gT2*b2)
      endwhere
!
      if (lfirst.and.ldt) then
         diffus_chi=diffus_chi+ cosbgT*gamma*Kgpara*exp(-p%lnrho)/p%cp*dxyz_2
         if (ldiagnos.and.idiag_dtchi/=0) then
            call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
         endif
      endif
!
    endsubroutine calc_heatcond_tensor
!***********************************************************************
    subroutine calc_heatcond_hubeny(df,p)
!
!  Vertically integrated heat flux from a thin globaldisc.
!  Taken from D'Angelo et al. 2003, ApJ, 599, 548, based on
!  the grey analytical model of Hubeny 1990, ApJ, 351, 632
!
!  07-feb-07/wlad+heidar : coded
!
      use EquationOfState, only: gamma,gamma_m1
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: tau,cooling,kappa,a1,a3
      real :: a2,kappa0,kappa0_cgs
!
      type (pencil_case) :: p
!
      intent(in) :: p
      intent(out) :: df
!
      if (headtt) print*,'enter heatcond hubeny'
!
      if (pretend_lnTT) call fatal_error("calc_heatcond_hubeny","not implemented when pretend_lnTT = T")
!
      kappa0_cgs=2e-6  !cm2/g
      kappa0=kappa0_cgs*unit_density/unit_length
      kappa=kappa0*p%TT**2
!
!  Optical Depth tau=kappa*rho*H
!  If we are using 2D, the pencil value p%rho is actually
!   sigma, the column density, sigma=rho*2*H
!
      if (nzgrid==1) then
         tau = .5*kappa*p%rho
      else
         call fatal_error("calc_heat_hubeny","opacity not yet implemented for 3D")
      endif
!
! Analytical gray description of Hubeny (1990)
! a1 is the optically thick contribution,
! a3 the optically thin one.
!
      a1=0.375*tau ; a2=0.433013 ; a3=0.25/tau
!
      cooling = 2*sigmaSB*p%TT**4/(a1+a2+a3)
!
      df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) - cool_fac*cooling
!
    endsubroutine calc_heatcond_hubeny
!***********************************************************************
    subroutine calc_heat_cool(df,p,Hmax)
!
!  Add combined heating and cooling to energy equation.
!
!  02-jul-02/wolf: coded
!
      use Diagnostics
      use Gravity
      use IO
      use HDF5IO, only: output_profile
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: Hmax
!
      real, dimension (nx) :: heat,prof,theta_profile
      real :: zbot,ztop,profile_buffer,xi,profile_cor
!
      intent(in) :: p
      intent(out) :: df
!
      if (pretend_lnTT) call fatal_error( &
          'calc_heat_cool','not implemented when pretend_lnTT = T')
!
!  Vertical gravity determines some heat/cool models.
!
      if (headtt) print*, 'calc_heat_cool: lgravz, lgravr= lgravx= lspherical_coords=', lgravz, lgravr,lgravx,lspherical_coords
!
!  Initialize heating/cooling term.
!
      heat=0.
!
!  Vertical gravity case: Heat at bottom, cool top layers
!
      if (lgravz .and. ( (luminosity/=0.) .or. (cool/=0.) ) ) then
        zbot=xyz0(3)
        ztop=xyz0(3)+Lxyz(3)
!
!  Add heat near bottom
!
!  Heating profile, normalised, so volume integral = 1
        prof = spread(exp(-0.5*((z(n)-zbot)/wheat)**2), 1, l2-l1+1) &
             /(sqrt(pi/2.)*wheat*Lx*Ly)
        heat = luminosity*prof
!  Smoothly switch on heating if required.
        if ((ttransient>0) .and. (t<ttransient)) then
          heat = heat * t*(2*ttransient-t)/ttransient**2
        endif
!
!  Allow for different cooling profile functions.
!  The gaussian default is rather broad and disturbs the entire interior.
!
        if (headtt) print*, 'cooling_profile: cooling_profile,z2,wcool=', &
            cooling_profile, z2, wcool
        select case (cooling_profile)
        case ('gaussian')
          prof = spread(exp(-0.5*((ztop-z(n))/wcool)**2),1,nx)
        case ('step')
          prof = spread(step(z(n),z2,wcool),1,nx)
        case ('cubic_step')
          prof = spread(cubic_step(z(n),z2,wcool),1,nx)
        endselect
        heat = heat - cool*prof*(p%cs2-cs2cool)/cs2cool
!
!  Write out cooling profile (during first time step only) and apply.
!  MR: later to be moved to initialization!
!
        if (m==m1) call output_profile('cooling_profile',z(n:n),prof(1:1),'z',lsave_name=(n==n1))
!
!  Write divergence of cooling flux.
!
        if (l1davgfirst) then
          if (idiag_dcoolz/=0) call xysum_mn_name_z(heat,idiag_dcoolz)
        endif
      endif
!
!  Spherical gravity case: heat at centre, cool outer layers.
!
      if (lgravr.and.(wheat/=0)) then
!  normalised central heating profile so volume integral = 1
        if (nzgrid == 1) then
          prof = exp(-0.5*(p%r_mn/wheat)**2) * (2*pi*wheat**2)**(-1.)  ! 2-D heating profile
        else
          prof = exp(-0.5*(p%r_mn/wheat)**2) * (2*pi*wheat**2)**(-1.5) ! 3-D one
        endif
        heat = luminosity*prof
        if (headt .and. lfirst .and. ip<=9) &
          call output_pencil('heat.dat',heat,1)
!
!  surface cooling: entropy or temperature
!  cooling profile; maximum = 1
!
!       prof = 0.5*(1+tanh((r_mn-1.)/wcool))
        if (rcool==0.) rcool=r_ext
        if (lcylindrical_coords) then
          prof = step(p%rcyl_mn,rcool,wcool)
        else
          prof = step(p%r_mn,rcool,wcool)
        endif
!
!  pick type of cooling
!
        select case (cooltype)
        case ('cs2', 'Temp')    ! cooling to reference temperature cs2cool
          heat = heat - cool*prof*(p%cs2-cs2cool)/cs2cool
        case ('cs2-rho', 'Temp-rho') ! cool to reference temperature cs2cool
                                     ! in a more time-step neutral manner
          heat = heat - cool*prof*(p%cs2-cs2cool)/cs2cool/p%rho1
        case ('entropy')        ! cooling to reference entropy (currently =0)
          heat = heat - cool*prof*(p%ss-0.)
        case ('shell')          !  heating/cooling at shell boundaries
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
            prof = step(p%r_mn,r_ext,wcool)     ! outer heating/cooling step
            heat = heat - cool_ext*prof*(p%cs2-cs2_ext)/cs2_ext
            prof = 1 - step(p%r_mn,r_int,wcool) ! inner heating/cooling step
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
!  Spherical gravity in spherical coordinate case:
!           heat at centre, cool outer layers.
!
      if (lgravx.and.lspherical_coords) then
        r_ext=x(l1)
        r_int=x(l2)
!  normalised central heating profile so volume integral = 1
        if (nzgrid == 1) then
          prof = exp(-0.5*(x(l1:l2)/wheat)**2) * (2*pi*wheat**2)**(-1.)  ! 2-D heating profile
        else
          prof = exp(-0.5*(x(l1:l2)/wheat)**2) * (2*pi*wheat**2)**(-1.5) ! 3-D one
        endif
        heat = luminosity*prof
        if (headt .and. lfirst .and. ip<=9) &
          call output_pencil('heat.dat',heat,1)
!
!  surface cooling: entropy or temperature
!  cooling profile; maximum = 1
!
!       prof = 0.5*(1+tanh((r_mn-1.)/wcool))
        if (rcool==0.) rcool=r_ext
        prof = step(x(l1:l2),rcool,wcool)
!
!
!  pick type of cooling
!
        select case (cooltype)
        case ('shell')          !  heating/cooling at shell boundaries
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
            prof = step(x(l1:l2),r_ext,wcool)     ! outer heating/cooling step
            heat = heat - cool_ext*prof*(p%cs2-cs2_ext)/cs2_ext
            prof = 1 - step(x(l1:l2),r_int,wcool) ! inner heating/cooling step
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
!  Add spatially uniform heating.
!
      if (heat_uniform/=0.) heat=heat+heat_uniform
!
!  Add cooling with constant time-scale to TTref_cool.
!
      if (tau_cool/=0.0) &
          heat=heat-p%rho*p%cp*gamma1*(p%TT-TTref_cool)/tau_cool
!
!  Add "coronal" heating (to simulate a hot corona).
!  Assume a linearly increasing reference profile.
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
!  Add heating and cooling to a reference temperature in a buffer
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
!  Add heating/cooling to entropy equation.
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + p%TT1*p%rho1*heat
      if (lfirst.and.ldt) Hmax=Hmax+heat*p%rho1
!
!  Heating/cooling related diagnostics.
!
      if (ldiagnos) then
        if (idiag_heatm/=0) call sum_mn_name(heat,idiag_heatm)
      endif
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
      use IO, only: output_pencil
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
      if (unit_system == 'cgs') unit_lnQ = unit_lnQ+alog(1.e13)
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
      if (cool_RTV > 0.) then
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
      elseif (cool_RTV < 0) then
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
      if (lfirst .and. ip == 13) &
           call output_pencil('rtv.dat',rtv_cool*exp(p%lnrho+p%lnTT),1)
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
    subroutine rprint_energy(lreset,lwrite)
!
!  reads and registers print parameters relevant to entropy
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Diagnostics
      use FArrayManager, only: farray_index_append
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      integer :: iname,inamez,inamey,inamex,inamexy,inamexz,irz,inamer
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtc=0; idiag_ethm=0; idiag_ethdivum=0; idiag_ssm=0; idiag_ss2m=0
        idiag_eem=0; idiag_ppm=0; idiag_csm=0; idiag_pdivum=0; idiag_heatm=0
        idiag_ugradpm=0; idiag_ethtot=0; idiag_dtchi=0; idiag_ssmphi=0
        idiag_fradbot=0; idiag_fradtop=0; idiag_TTtop=0
        idiag_yHmax=0; idiag_yHm=0; idiag_TTmax=0; idiag_TTmin=0; idiag_TTm=0
        idiag_fconvm=0; idiag_fconvz=0; idiag_dcoolz=0; idiag_fradz=0; idiag_fturbz=0
        idiag_ssmz=0; idiag_ssmy=0; idiag_ssmx=0; idiag_ssmr=0; idiag_TTmr=0
        idiag_TTmx=0; idiag_TTmy=0; idiag_TTmz=0; idiag_TTmxy=0; idiag_TTmxz=0
        idiag_uxTTmz=0; idiag_uyTTmz=0; idiag_uzTTmz=0; idiag_cs2mphi=0
        idiag_ssmxy=0; idiag_ssmxz=0
        cformv=''
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtc',idiag_dtc)
        call parse_name(iname,cname(iname),cform(iname),'dtchi',idiag_dtchi)
        call parse_name(iname,cname(iname),cform(iname),'ethtot',idiag_ethtot)
        call parse_name(iname,cname(iname),cform(iname),'ethdivum',idiag_ethdivum)
        call parse_name(iname,cname(iname),cform(iname),'ethm',idiag_ethm)
        call parse_name(iname,cname(iname),cform(iname),'ssm',idiag_ssm)
        call parse_name(iname,cname(iname),cform(iname),'ss2m',idiag_ss2m)
        call parse_name(iname,cname(iname),cform(iname),'eem',idiag_eem)
        call parse_name(iname,cname(iname),cform(iname),'ppm',idiag_ppm)
        call parse_name(iname,cname(iname),cform(iname),'pdivum',idiag_pdivum)
        call parse_name(iname,cname(iname),cform(iname),'heatm',idiag_heatm)
        call parse_name(iname,cname(iname),cform(iname),'csm',idiag_csm)
        call parse_name(iname,cname(iname),cform(iname),'ugradpm',idiag_ugradpm)
        call parse_name(iname,cname(iname),cform(iname),'fradbot',idiag_fradbot)
        call parse_name(iname,cname(iname),cform(iname),'fradtop',idiag_fradtop)
        call parse_name(iname,cname(iname),cform(iname),'TTtop',idiag_TTtop)
        call parse_name(iname,cname(iname),cform(iname),'yHm',idiag_yHm)
        call parse_name(iname,cname(iname),cform(iname),'yHmax',idiag_yHmax)
        call parse_name(iname,cname(iname),cform(iname),'TTm',idiag_TTm)
        call parse_name(iname,cname(iname),cform(iname),'TTmax',idiag_TTmax)
        call parse_name(iname,cname(iname),cform(iname),'TTmin',idiag_TTmin)
        call parse_name(iname,cname(iname),cform(iname),'TTp',idiag_TTp)
        call parse_name(iname,cname(iname),cform(iname),'fconvm',idiag_fconvm)
      enddo
!
!  check for those quantities for which we want yz-averages
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ssmx',idiag_ssmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'TTmx',idiag_TTmx)
      enddo
!
!  check for those quantities for which we want xz-averages
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'ssmy',idiag_ssmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'TTmy',idiag_TTmy)
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
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uxTTmz',idiag_uxTTmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uyTTmz',idiag_uyTTmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzTTmz',idiag_uzTTmz)
      enddo
!
      do inamer=1,nnamer
         call parse_name(inamer,cnamer(inamer),cformr(inamer),'ssmr',idiag_ssmr)
         call parse_name(inamer,cnamer(inamer),cformr(inamer),'TTmr',idiag_TTmr)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'TTmxy',idiag_TTmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'ssmxy',idiag_ssmxy)
      enddo
!
!  check for those quantities for which we want y-averages
!
      do inamexz=1,nnamexz
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'TTmxz',idiag_TTmxz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'ssmxz',idiag_ssmxz)
      enddo
!
!  check for those quantities for which we want phi-averages
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'ssmphi',idiag_ssmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'cs2mphi',idiag_cs2mphi)
      enddo
!       
!  check for those quantities for which we want video slices
!       
      if (lwrite_slices) then
        where(cnamev=='ss') cformv='DEFINED'
      endif
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        call farray_index_append('i_dtc',idiag_dtc)
        call farray_index_append('i_dtchi',idiag_dtchi)
        call farray_index_append('i_ethtot',idiag_ethtot)
        call farray_index_append('i_ethdivum',idiag_ethdivum)
        call farray_index_append('i_ethm',idiag_ethm)
        call farray_index_append('i_ssm',idiag_ssm)
        call farray_index_append('i_ss2m',idiag_ss2m)
        call farray_index_append('i_eem',idiag_eem)
        call farray_index_append('i_ppm',idiag_ppm)
        call farray_index_append('i_pdivum',idiag_pdivum)
        call farray_index_append('i_heatm',idiag_heatm)
        call farray_index_append('i_csm',idiag_csm)
        call farray_index_append('i_fconvm',idiag_fconvm)
        call farray_index_append('i_ugradpm',idiag_ugradpm)
        call farray_index_append('i_fradbot',idiag_fradbot)
        call farray_index_append('i_fradtop',idiag_fradtop)
        call farray_index_append('i_TTtop',idiag_TTtop)
        call farray_index_append('i_ssmphi',idiag_ssmphi)
        call farray_index_append('i_cs2mphi',idiag_cs2mphi)
        call farray_index_append('i_fturbz',idiag_fturbz)
        call farray_index_append('i_fconvz',idiag_fconvz)
        call farray_index_append('i_dcoolz',idiag_dcoolz)
        call farray_index_append('i_fradz',idiag_fradz)
        call farray_index_append('i_ssmz',idiag_ssmz)
        call farray_index_append('i_TTmz',idiag_TTmz)
        call farray_index_append('i_uxTTmz',idiag_uxTTmz)
        call farray_index_append('i_uyTTmz',idiag_uyTTmz)
        call farray_index_append('i_uzTTmz',idiag_uzTTmz)
        call farray_index_append('i_ssmr',idiag_ssmr)
        call farray_index_append('i_TTmr',idiag_TTmr)
        call farray_index_append('nname',nname)
        call farray_index_append('iss',iss)
        call farray_index_append('i_yHmax',idiag_yHmax)
        call farray_index_append('i_yHm',idiag_yHm)
        call farray_index_append('i_TTmax',idiag_TTmax)
        call farray_index_append('i_TTmin',idiag_TTmin)
        call farray_index_append('i_TTm',idiag_TTm)
        call farray_index_append('i_TTp',idiag_TTp)
        call farray_index_append('iyH',iyH)
        call farray_index_append('ilnTT',ilnTT)
        call farray_index_append('i_TTmxy',idiag_TTmxy)
        call farray_index_append('i_TTmxz',idiag_TTmxz)
        call farray_index_append('i_ssmxy',idiag_ssmxy)
        call farray_index_append('i_ssmxz',idiag_ssmxz)
      endif
!
    endsubroutine rprint_energy
!***********************************************************************
    subroutine get_slices_energy(f,slices)
!
!  Write slices for animation of Energy variables.
!
!  26-jul-06/tony: coded
!
      use EquationOfState, only: eoscalc, ilnrho_ss
      use Slices_methods, only: assign_slices_scal
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      real :: tmpval
      integer :: l
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Entropy.
!
        case ('ss'); call assign_slices_scal(slices,f,iss)

      endselect
!
    endsubroutine get_slices_energy
!***********************************************************************
    subroutine fill_farray_pressure(f)
!
!  18-feb-10/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine fill_farray_pressure
!***********************************************************************
    subroutine calc_heatcond_zprof(zprof_hcond,zprof_glhc)
!
!  calculate z-profile of heat conduction for multilayer setup
!
!  12-jul-05/axel: coded
!
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
    subroutine heatcond(hcond,p)
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
      use Gravity, only: z1, z2
!
      real, dimension (nx) :: hcond
      type (pencil_case)   :: p
!
      if (lgravz) then
        if (lmultilayer) then
          hcond = 1. + (hcond1-1.)*step(p%z_mn,z1,-widthss) &
                     + (hcond2-1.)*step(p%z_mn,z2,widthss)
          hcond = hcond0*hcond
        else
          hcond=Kbot
        endif
      else
        if (lmultilayer) then
          hcond = 1. + (hcond1-1.)*step(p%r_mn,r_bcz,-widthss) &
                     + (hcond2-1.)*step(p%r_mn,r_ext,widthss)
          hcond = hcond0*hcond
        else
          hcond = hcond0
        endif
      endif
!
    endsubroutine heatcond
!***********************************************************************
    subroutine gradloghcond(glhc,p)
!
!  calculate grad(log hcond), where hcond is the heat conductivity
!  NB: *Must* be in sync with heatcond() above.
!
!  23-jan-2002/wolf: coded
!
      use Sub, only: der_step
      use Gravity, only: z1, z2
!
      real, dimension (nx,3) :: glhc
      real, dimension (nx)   :: dhcond
      type (pencil_case)     :: p
!
      if (lgravz) then
        if (lmultilayer) then
          glhc(:,1:2) = 0.
          glhc(:,3) = (hcond1-1.)*der_step(p%z_mn,z1,-widthss) &
                      + (hcond2-1.)*der_step(p%z_mn,z2,widthss)
          glhc(:,3) = hcond0*glhc(:,3)
        else
          glhc = 0.
        endif
      else
        if (lmultilayer) then
          dhcond=(hcond1-1.)*der_step(p%r_mn,r_bcz,-widthss) &
                 + (hcond2-1.)*der_step(p%r_mn,r_ext,widthss)
          dhcond=hcond0*dhcond
          glhc(:,1) = x(l1:l2)/p%r_mn*dhcond
          glhc(:,2) = y(m)/p%r_mn*dhcond
          if (lcylinder_in_a_box) then
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
      endif
!
      if (lspherical_coords) then
        chit_prof = 1 + (chit_prof1-1)*step(x(l1:l2),xbot,-widthss)
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
      endif
!
      if (lspherical_coords) then
        glchit_prof(:,1) = (chit_prof1-1)*der_step(x(l1:l2),xbot,-widthss)
        glchit_prof(:,2:3) = 0.
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
      use Io, only:  output_pencil
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
      if (it == 1) then
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
      if (z(n) < b_z(1) ) then
        lnTTor = b_lnT(1)
      elseif (z(n) >= b_z(150)) then
        lnTTor = b_lnT(150)
      else
        do i=1,149
          if (z(n) >= b_z(i) .and. z(n) < b_z(i+1)) then
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
      if (dt /= 0 )  then
         newton = tdown/gamma/dt*((lnTTor - p%lnTT) )
         newton = newton * exp(allp*(-abs(p%lnrho-lnrho0)))
!
!  Add newton cooling term to entropy
!
         if (lfirst .and. ip == 13) &
              call output_pencil('newton.dat',newton*exp(p%lnrho+p%lnTT),1)
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
    subroutine strat_MLT (rhotop,mixinglength_flux,lnrhom,tempm,rhobot)
!
! 04-mai-06/dintrans: called by 'mixinglength' to iterate the MLT
! equations until rho=1 at the bottom of convection zone (z=z1)
! see Eqs. (20-21-22) in Brandenburg et al., AN, 326 (2005)
!
      use Gravity, only: z1,z2,gravz
      use EquationOfState, only: gamma, gamma_m1
!
      real, dimension (nzgrid) :: zz, lnrhom, tempm
      real :: rhotop, zm, ztop, dlnrho, dtemp, &
              mixinglength_flux, lnrhobot, rhobot
      real :: del, delad, fr_frac, fc_frac, fc, polyad
      integer :: nbot1, nbot2
!
!  inital values at the top
!
      lnrhom(1)=alog(rhotop)
      tempm(1)=cs2top/gamma_m1
      ztop=xyz0(3)+Lxyz(3)
      zz(1)=ztop
!
      polyad=1./gamma_m1
      delad=1.-1./gamma
      fr_frac=delad*(mpoly0+1.)
      fc_frac=1.-fr_frac
      fc=fc_frac*mixinglength_flux
!     print*,'fr_frac, fc_frac, fc=',fr_frac,fc_frac,fc,delad
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
            del=delad+1.5*(fc/ &
                (exp(lnrhom(iz-1))*(gamma_m1*tempm(iz-1))**1.5))**.6666667
          else
! upper zone=isothermal stratification
            del=0.
          endif
        endif
        dtemp=gamma*polyad*gravz*del
        dlnrho=gamma*polyad*gravz*(1.-del)/tempm(iz-1)
        tempm(iz)=tempm(iz-1)-dtemp*dz
        lnrhom(iz)=lnrhom(iz-1)-dlnrho*dz
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
      use EquationOfState, only: eoscalc, ilnrho_lnTT, lnrho0, get_cp1
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: lnrho,lnTT,TT,ss,r_mn
      real :: beta0,beta1,TT_crit,cp1
      real :: lnrho_int,lnrho_ext,lnrho_crit
!
      if (headtt) print*,'r_bcz in entropy.f90=',r_bcz
!
!  beta is the temperature gradient
!  1/beta = -(g/cp) /[(1-1/gamma)*(m+1)]
!
      call get_cp1(cp1)
      beta0=-cp1*g0/(mpoly0+1)*gamma/gamma_m1
      beta1=-cp1*g0/(mpoly1+1)*gamma/gamma_m1
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
    subroutine star_heat(f,lcompute_grav)
!
!  Initialize energy for two superposed polytropes with a central heating
!
!  20-dec-06/dintrans: coded
!  28-nov-07/dintrans: merged with strat_heat_grav
!
    use EquationOfState, only: gamma, gamma_m1, rho0, lnrho0, cs20, get_soundspeed,eoscalc, ilnrho_TT
    use FArrayManager
    use Sub, only: step, interp1, erfunc
!
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f
    integer, parameter   :: nr=100
    integer              :: i,l,iter
    real, dimension (nr) :: r,lnrho,temp,hcond
    real                 :: u,r_mn,lnrho_r,temp_r,cs2,ss,lumi,g
    real                 :: rhotop, rbot,rt_old,rt_new,rhobot,rb_old,rb_new,crit,r_max
! variables for the gravity profile
    logical, optional    :: lcompute_grav
    integer, pointer     :: iglobal_gg
    real, allocatable, dimension (:)   :: rr_mn,u_mn,lumi_mn,g_r
!
    if (present(lcompute_grav)) then
      print*,'only compute the gravity profile',lcompute_grav
      call farray_use_global('gg',iglobal_gg)
      allocate (rr_mn(nx),u_mn(nx),lumi_mn(nx),g_r(nx))
      do imn=1,ny*nz
        m=mm(imn)
        n=nn(imn)
        if (nzgrid == 1) then
          rr_mn=sqrt(x(l1:l2)**2+y(m)**2)
        else
          rr_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
        endif
        u_mn=rr_mn/sqrt(2.)/wheat
        if (nzgrid == 1) then
          lumi_mn=luminosity*(1.-exp(-u_mn**2))
          g_r=-lumi_mn/(2.*pi*rr_mn)*(mpoly0+1.)/hcond0*gamma_m1/gamma
          f(l1:l2,m,n,iglobal_gg+2) = 0.                 ! g_z=0
        else
          lumi_mn=luminosity*(erfunc(u_mn)-2.*u_mn/sqrt(pi)*exp(-u_mn**2))
          g_r=-lumi_mn/(4.*pi*rr_mn**2)*(mpoly0+1.)/hcond0*gamma_m1/gamma
          f(l1:l2,m,n,iglobal_gg+2) = z(n)/rr_mn*g_r     ! g_z
        endif
        f(l1:l2,m,n,iglobal_gg)   = x(l1:l2)/rr_mn*g_r   ! g_x
        f(l1:l2,m,n,iglobal_gg+1) = y(m)/rr_mn*g_r       ! g_y
     enddo
     deallocate(rr_mn,u_mn,lumi_mn,g_r)
     return
    endif
!
!  the bottom value that we want for density at r=r_bcz, actually given by rho0
    rbot=rho0
    rt_old=0.1*rbot
    rt_new=0.12*rbot
!
!  need to iterate for rhobot=1
!  produce first estimate
    rhotop=rt_old
    call strat_heat(lnrho,temp,rhotop,rhobot)
    rb_old=rhobot
!
!  next estimate
    rhotop=rt_new
    call strat_heat(lnrho,temp,rhotop,rhobot)
    rb_new=rhobot
!
    do iter=1,10
!  new estimate
      rhotop=rt_old+(rt_new-rt_old)/(rb_new-rb_old)*(rbot-rb_old)
!
      crit=abs(rhotop/rt_new-1.)
      if (crit<=1e-4) exit
      call strat_heat(lnrho,temp,rhotop,rhobot)
!
!  update new estimates
      rt_old=rt_new
      rb_old=rb_new
      rt_new=rhotop
      rb_new=rhobot
    enddo
!
! uncomment to force rhotop=rho0 and just compute the corresponding setup
!   rhotop=rho0
!   call strat_heat(lnrho,temp,rhotop,rhobot)
!
    print*,'- iteration completed: rhotop,crit=',rhotop,crit
!
!  redefine rho0 and lnrho0 (important for eoscalc!)
    rho0=rhotop
    lnrho0=log(rhotop)
    T0=cs20/gamma_m1
    print*,'final rho0, lnrho0, T0=',rho0, lnrho0, T0
!
! define the radial grid r=[0,r_max]
    if (nzgrid == 1) then
      r_max=sqrt(xyz1(1)**2+xyz1(2)**2)
    else
      r_max=sqrt(xyz1(1)**2+xyz1(2)**2+xyz1(3)**2)
    endif
    do i=1,nr
      r(i)=r_max*float(i-1)/(nr-1)
    enddo
!
    do imn=1,ny*nz
      n=nn(imn)
      m=mm(imn)
!
      do l=l1,l2
        if (nzgrid == 1) then
          r_mn=sqrt(x(l)**2+y(m)**2)
        else
          r_mn=sqrt(x(l)**2+y(m)**2+z(n)**2)
        endif
        lnrho_r=interp1(r,lnrho,nr,r_mn)
        temp_r=interp1(r,temp,nr,r_mn)
        f(l,m,n,ilnrho)=lnrho_r
        call eoscalc(ilnrho_TT,lnrho_r,temp_r,ss=ss)
        f(l,m,n,iss)=ss
      enddo
    enddo
!
    if (lroot) then
      hcond=1.+(hcond1-1.)*step(r,r_bcz,-widthss) &
              +(hcond2-1.)*step(r,r_ext,widthss)
      hcond=hcond0*hcond
      print*,'--> writing initial setup to data/proc0/setup.dat'
      open(unit=11,file=trim(directory)//'/setup.dat')
      write(11,'(a1,a5,5a14)') '#','r','rho','ss','cs2','grav','hcond'
      do i=nr,1,-1
        u=r(i)/sqrt(2.)/wheat
        if (nzgrid == 1) then
          lumi=luminosity*(1.-exp(-u**2))
        else
          lumi=luminosity*(erfunc(u)-2.*u/sqrt(pi)*exp(-u**2))
        endif
        if (r(i) /= 0.) then
          if (nzgrid == 1) then
            g=-lumi/(2.*pi*r(i))*(mpoly0+1.)/hcond0*gamma_m1/gamma
          else
            g=-lumi/(4.*pi*r(i)**2)*(mpoly0+1.)/hcond0*gamma_m1/gamma
          endif
        else
          g=0.
        endif
        call get_soundspeed(log(temp(i)),cs2)
        call eoscalc(ilnrho_TT,lnrho(i),temp(i),ss=ss)
        write(11,'(f6.3,4e14.5,1pe14.5)') r(i),exp(lnrho(i)),ss,cs2,g,hcond(i)
      enddo
      close(11)
    endif
!
    endsubroutine star_heat
!***********************************************************************
    subroutine strat_heat(lnrho,temp,rhotop,rhobot)
!
!  compute the radial stratification for two superposed polytropic
!  layers and a central heating
!
!  17-jan-07/dintrans: coded
!
    use EquationOfState, only: gamma, gamma_m1, lnrho0, cs20
    use Sub, only: step, erfunc, interp1
!
    integer, parameter   :: nr=100
    integer              :: i
    real, dimension (nr) :: r,lumi,hcond,g,lnrho,temp
    real                 :: dtemp,dlnrho,dr,u,rhotop,rhobot,lnrhobot,r_max
!
! define the radial grid r=[0,r_max]
    if (nzgrid == 1) then
      r_max=sqrt(xyz1(1)**2+xyz1(2)**2)
    else
      r_max=sqrt(xyz1(1)**2+xyz1(2)**2+xyz1(3)**2)
    endif
    do i=1,nr
      r(i)=r_max*float(i-1)/(nr-1)
    enddo
!
! luminosity and gravity radial profiles
    lumi(1)=0. ; g(i)=0.
    do i=2,nr
      u=r(i)/sqrt(2.)/wheat
      if (nzgrid == 1) then
        lumi(i)=luminosity*(1.-exp(-u**2))
        g(i)=-lumi(i)/(2.*pi*r(i))*(mpoly0+1.)/hcond0*gamma_m1/gamma
      else
        lumi(i)=luminosity*(erfunc(u)-2.*u/sqrt(pi)*exp(-u**2))
        g(i)=-lumi(i)/(4.*pi*r(i)**2)*(mpoly0+1.)/hcond0*gamma_m1/gamma
      endif
    enddo
!
! radiative conductivity profile
    hcond1=(mpoly1+1.)/(mpoly0+1.)
    hcond2=(mpoly2+1.)/(mpoly0+1.)
    hcond=1.+(hcond1-1.)*step(r,r_bcz,-widthss) &
            +(hcond2-1.)*step(r,r_ext,widthss)
    hcond=hcond0*hcond
!
! start from surface values for rho and temp
    temp(nr)=cs20/gamma_m1 ; lnrho(nr)=alog(rhotop)
    dr=r(2)
    do i=nr-1,1,-1
      if (r(i+1) > r_ext) then
! isothermal exterior
        dtemp=0.
        dlnrho=-gamma*g(i+1)/cs20
      elseif (r(i+1) > r_bcz) then
! convective zone
!       if (nzgrid == 1) then
!         dtemp=lumi(i+1)/(2.*pi*r(i+1))/hcond(i+1)
!       else
!         dtemp=lumi(i+1)/(4.*pi*r(i+1)**2)/hcond(i+1)
!       endif
!       dlnrho=mpoly0*dtemp/temp(i+1)
! force adiabatic stratification with m0=3/2 (assume cp=1)
        dtemp=-g(i+1)
        dlnrho=3./2.*dtemp/temp(i+1)
      else
! radiative zone
        if (nzgrid == 1) then
          dtemp=lumi(i+1)/(2.*pi*r(i+1))/hcond(i+1)
        else
          dtemp=lumi(i+1)/(4.*pi*r(i+1)**2)/hcond(i+1)
        endif
        dlnrho=mpoly1*dtemp/temp(i+1)
      endif
      temp(i)=temp(i+1)+dtemp*dr
      lnrho(i)=lnrho(i+1)+dlnrho*dr
    enddo
!
! find the value of rhobot at the bottom of convection zone
!
!   lnrhobot=interp1(r,lnrho,nr,r_max)
!   lnrhobot=interp1(r,lnrho,nr,r_ext)
    lnrhobot=interp1(r,lnrho,nr,r_bcz)
    rhobot=exp(lnrhobot)
    print*,'find rhobot=',rhobot
!
    endsubroutine strat_heat
!***********************************************************************
    subroutine cylind_layers(f)
!
!  17-mar-07/dintrans: coded
!  Initialise ss in a cylindrical ring using 2 superposed polytropic layers
!
      use Gravity, only: gravz, g0
      use EquationOfState, only: lnrho0,cs20,gamma,gamma_m1,cs2top,cs2bot, &
                                 get_cp1,eoscalc,ilnrho_TT
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: lnrho, TT, ss
      real :: beta0, beta1, TT_bcz, TT_ext, TT_int
      real :: cp1, lnrho_bcz
!
      if (headtt) print*,'r_bcz in cylind_layers=', r_bcz
!
!  beta is the (negative) temperature gradient
!  beta = (g/cp) 1./[(1-1/gamma)*(m+1)]
!
      call get_cp1(cp1)
      beta0=-cp1*g0/(mpoly0+1)*gamma/gamma_m1
      beta1=-cp1*g0/(mpoly1+1)*gamma/gamma_m1
      TT_ext=cs20/gamma_m1
      TT_bcz=TT_ext+beta0*(r_bcz-r_ext)
      TT_int=TT_bcz+beta1*(r_int-r_bcz)
      cs2top=cs20
      cs2bot=gamma_m1*TT_int
      lnrho_bcz=lnrho0+mpoly0*log(TT_bcz/TT_ext)
!
      do m=m1,m2
      do n=n1,n2
!
!  convective layer
!
        where (rcyl_mn <= r_ext .AND. rcyl_mn > r_bcz)
          TT=TT_ext+beta0*(rcyl_mn-r_ext)
          lnrho=lnrho0+mpoly0*log(TT/TT_ext)
        endwhere
!
!  radiative layer
!
        where (rcyl_mn <= r_bcz)
          TT=TT_bcz+beta1*(rcyl_mn-r_bcz)
          lnrho=lnrho_bcz+mpoly1*log(TT/TT_bcz)
        endwhere
!
        f(l1:l2,m,n,ilnrho)=lnrho
        call eoscalc(ilnrho_TT,lnrho,TT,ss=ss)
        f(l1:l2,m,n,iss)=ss
      enddo
      enddo
!
    endsubroutine cylind_layers
!***********************************************************************
    subroutine single_polytrope(f)
!
!  06-sep-07/dintrans: coded a single polytrope of index mpoly0
!  Note: both entropy and density are initialized there (compared to layer_ss)
!
      use Gravity, only: gravz, g0
      use EquationOfState, only: eoscalc, ilnrho_TT, get_cp1, &
                                 gamma_m1, lnrho0
      use SharedVariables, only: get_shared_variable
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: lnrho, TT, ss, z_mn
      real :: beta, cp1, zbot, ztop, TT0
      real, pointer :: gravx
!
!  beta is the (negative) temperature gradient
!  beta = (g/cp) 1./[(1-1/gamma)*(m+1)]
!
      call get_cp1(cp1)
      if (lcylindrical_coords) then
        call get_shared_variable('gravx', gravx, caller='single_polytrope')
        beta=cp1*gamma/gamma_m1*gravx/(mpoly0+1)
      else
        beta=cp1*gamma/gamma_m1*gravz/(mpoly0+1)
        ztop=xyz0(3)+Lxyz(3)
        zbot=xyz0(3)
      endif
      TT0=cs20/gamma_m1
!
!  set initial condition (first in terms of TT, and then in terms of ss)
!
      do m=m1,m2
      do n=n1,n2
        if (lcylindrical_coords) then
          TT = TT0+beta*(rcyl_mn-r_ext)
        else
          z_mn = spread(z(n),1,nx)
          TT = TT0+beta*(z_mn-ztop)
        endif
        lnrho=lnrho0+mpoly0*log(TT/TT0)
        f(l1:l2,m,n,ilnrho)=lnrho
        call eoscalc(ilnrho_TT,lnrho,TT,ss=ss)
        f(l1:l2,m,n,iss)=ss
      enddo
      enddo
      cs2top=cs20
      if (lcylindrical_coords) then
        cs2bot=gamma_m1*(TT0+beta*(r_int-r_ext))
      else
        cs2bot=gamma_m1*(TT0+beta*(zbot-ztop))
      endif
!
    endsubroutine single_polytrope
!***********************************************************************
    subroutine calc_heatcond_ADI(finit,f)
!
      implicit none
!
      real, dimension(mx,my,mz,mfarray) :: finit,f
!
      call keep_compiler_quiet(finit)
      call keep_compiler_quiet(f)
!
    endsubroutine calc_heatcond_ADI
!***********************************************************************
    subroutine calc_diagnostics_energy(f,p)

      real, dimension (mx,my,mz,mfarray) :: f
      type(pencil_case) :: p

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)

    endsubroutine calc_diagnostics_energy
!***********************************************************************
    subroutine read_hcond(hcond,glhc)
!
!  read radial profiles of hcond and glhc from an ascii-file.
!
!  11-jun-09/pjk: coded
!
      use Mpicomm, only: stop_it
!
      real, dimension(nx) :: hcond
      real, dimension(nx,3) :: glhc
      integer, parameter :: ntotal=nx*nprocx
      real :: var1,var2
      logical :: exist
      integer :: stat
!
!  read hcond and glhc and write into an array
!  if file is not found in run directory, search under trim(directory)
!
      inquire(file='hcond_glhc.dat',exist=exist)
      if (exist) then
        open(31,file='hcond_glhc.dat')
      else
        inquire(file=trim(directory)//'/hcond_glhc.ascii',exist=exist)
        if (exist) then
          open(31,file=trim(directory)//'/hcond_glhc.ascii')
        else
          call stop_it('read_hcond: *** error *** - no input file')
        endif
      endif
!
!  read profiles
!
      do n=1,ntotal
        read(31,*,iostat=stat) var1,var2
        if (stat>=0) then
          if (ip<5) print*,"hcond, glhc: ",var1,var2
          hcond(n)=var1
          glhc(n,1)=var2
          glhc(n,2)=0.
          glhc(n,3)=0.
        else
          exit
        endif
      enddo
!
      close(31)
!
    endsubroutine read_hcond
!***********************************************************************
endmodule Energy
