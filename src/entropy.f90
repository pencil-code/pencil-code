! $Id$
!
!  This module takes care of evolving the entropy.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lentropy = .true.
! CPARAM logical, parameter :: ltemperature = .false.
! CPARAM logical, parameter :: lthermal_energy = .false.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED ugss; Ma2; fpres(3); uglnTT; transprhos
!
!***************************************************************
module Entropy
!
  use Cdata
  use Cparam
  use EquationOfState, only: gamma, gamma_m1, gamma_inv, cs20, cs2top, cs2bot, &
                         isothtop, mpoly0, mpoly1, mpoly2, cs2cool, &
                         beta_glnrho_global, cs2top_ini, dcs2top_ini
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'entropy.h'
!
  real :: entropy_floor = 0.
  real :: radius_ss=0.1, ampl_ss=0.0, widthss=2*epsi, epsilon_ss=0.0
  real :: luminosity=0.0, wheat=0.1, cool=0.0, zcool=0.0, rcool=0.0, wcool=0.1
  real :: rcool1=0.0, rcool2=0.0, deltaT=0.0
  real :: TT_int, TT_ext, cs2_int, cs2_ext
  real :: cool_int=0.0, cool_ext=0.0, ampl_TT=0.0
  real, target :: chi=0.0
  real :: chi_t=0.0, chi_shock=0.0, chi_hyper3=0.0
  real :: chi_hyper3_mesh=5.0, chi_th=0.0, chi_rho=0.0
  real :: Kgperp=0.0, Kgpara=0.0, tdown=0.0, allp=2.0, TT_powerlaw=1.0
  real :: ss_left=1.0, ss_right=1.0
  real :: ss0=0.0, khor_ss=1.0, ss_const=0.0
  real :: pp_const=0.0
  real :: tau_ss_exterior=0.0, T0=0.0
  real :: mixinglength_flux=0.0
  real :: center1_x=0.0, center1_y=0.0, center1_z=0.0
  real :: center2_x=0.0, center2_y=0.0, center2_z=0.0
  real :: kx_ss=1.0, ky_ss=1.0, kz_ss=1.0
  real :: thermal_background=0.0, thermal_peak=0.0, thermal_scaling=1.0
  real :: cool_fac=1.0, chiB=0.0
  real, dimension(3) :: chi_hyper3_aniso=0.0
  real, target :: hcond0=impossible, hcond1=impossible
  real, target :: hcondxbot=impossible, hcondxtop=impossible
  real, target :: hcondzbot=impossible, hcondztop=impossible
  real, target :: Fbot=impossible, FbotKbot=impossible
  real, target :: Ftop=impossible, FtopKtop=impossible
  real :: Kbot=impossible, Ktop=impossible
  real :: hcond2=impossible
  real :: chit_prof1=1.0, chit_prof2=1.0
  real :: chit_aniso=0.0, chit_aniso_prof1=1.0, chit_aniso_prof2=1.0
  real :: tau_cor=0.0, TT_cor=0.0, z_cor=0.0
  real :: tauheat_buffer=0.0, TTheat_buffer=0.0
  real :: zheat_buffer=0.0, dheat_buffer1=0.0
  real :: heat_uniform=0.0, cool_uniform=0.0, cool_newton=0.0, cool_RTV=0.0
  real :: deltaT_poleq=0.0, beta_hand=1.0, r_bcz=0.0
  real :: tau_cool=0.0, tau_diff=0.0, TTref_cool=0.0, tau_cool2=0.0
  real :: tau_cool_ss=0.0
  real :: cs0hs=0.0, H0hs=0.0, rho0hs=0.0, rho0ts=impossible
  real :: ss_volaverage=0., ss_volaverage_tmp
  real, target :: T0hs=impossible
  real, dimension(mx), target :: zrho
  real :: rho0ts_cgs=1.67262158e-24, T0hs_cgs=7.202e3
  real :: xbot=0.0, xtop=0.0, alpha_MLT=0.0, xbot_aniso=0.0, xtop_aniso=0.0
  real :: zz1=impossible, zz2=impossible
  real, target :: hcond0_kramers=0.0, nkramers=0.0
  integer, parameter :: nheatc_max=4
  integer :: iglobal_hcond=0
  integer :: iglobal_glhc=0
  integer :: ippaux=0
  logical :: lturbulent_heat=.false.
  logical :: lheatc_Kprof=.false., lheatc_Kconst=.false.
  logical, target :: lheatc_chiconst=.false.
  logical :: lheatc_tensordiffusion=.false., lheatc_spitzer=.false.
  logical :: lheatc_hubeny=.false.,lheatc_sqrtrhochiconst=.false.
  logical, target :: lheatc_kramers=.false.
  logical :: lheatc_corona=.false.,lheatc_chitherm=.false.
  logical :: lheatc_shock=.false., lheatc_hyper3ss=.false.
  logical :: lheatc_hyper3ss_polar=.false., lheatc_hyper3ss_aniso=.false.
  logical :: lheatc_hyper3ss_mesh=.false.
  logical :: lcooling_general=.false.
  logical :: lupw_ss=.false.
  logical :: lcalc_ssmean=.false., lcalc_ss_volaverage=.false.
  logical :: lcalc_cs2mean=.false., lcalc_cs2mz_mean=.false.
  logical :: lcalc_ssmeanxy=.false.
  logical, target :: lmultilayer=.true.
  logical :: ladvection_entropy=.true.
  logical, pointer :: lpressuregradient_gas
  logical :: lviscosity_heat=.true.
  logical :: lfreeze_sint=.false.,lfreeze_sext=.false.
  logical :: lhcond_global=.false.,lchit_aniso_simplified=.false.
  logical :: lfpres_from_pressure=.false.
  logical :: lconvection_gravx=.false.
  logical :: ltau_cool_variable=.false.
  logical :: lprestellar_cool_iso=.false.
  logical, save :: lfirstcall_hcond=.true.
  character (len=labellen), dimension(ninit) :: initss='nothing'
  character (len=labellen) :: borderss='nothing'
  character (len=labellen) :: pertss='zero'
  character (len=labellen) :: cooltype='Temp',cooling_profile='gaussian'
  character (len=labellen), dimension(nheatc_max) :: iheatcond='nothing'
  character (len=intlen) :: iinit_str
  real, dimension (mz),   save :: hcond_zprof,chit_zprof
  real, dimension (mz,3), save :: gradloghcond_zprof,gradlogchit_zprof
  real, dimension (mx),   save :: hcond_xprof,chit_xprof
  real, dimension (mx,3), save :: gradloghcond_xprof
!
!  xy-averaged field
!
  real, dimension (mz) :: ssmz,cs2mz
  real, dimension (nz,3) :: gssmz
  real, dimension (nz) :: del2ssmz
  real, dimension (mx) :: cs2mx, ssmx
  real, dimension (mx,my) :: cs2mxy, ssmxy
!
!  Input parameters.
!
  namelist /entropy_init_pars/ &
      initss, pertss, grads0, radius_ss, ampl_ss, widthss, epsilon_ss, &
      mixinglength_flux, chi_t, chi_th, chi_rho, pp_const, ss_left, ss_right, &
      ss_const, mpoly0, mpoly1, mpoly2, isothtop, khor_ss, &
      thermal_background, thermal_peak, thermal_scaling, cs2cool, center1_x, &
      center1_y, center1_z, center2_x, center2_y, center2_z, T0, ampl_TT, &
      kx_ss, ky_ss, kz_ss, beta_glnrho_global, ladvection_entropy, &
      lviscosity_heat, r_bcz, luminosity, wheat, hcond0, tau_cool, &
      tau_cool_ss, &
      TTref_cool, lhcond_global, cool_fac, cs0hs, H0hs, rho0hs, tau_cool2, &
      rho0ts, T0hs, lconvection_gravx, Fbot, hcond0_kramers, nkramers, &
      alpha_MLT, lprestellar_cool_iso
!
!  Run parameters.
!
  namelist /entropy_run_pars/ &
      hcond0, hcond1, hcond2, widthss, borderss, mpoly0, mpoly1, mpoly2, &
      luminosity, wheat, cooling_profile, cooltype, cool, cs2cool, rcool, &
      rcool1, rcool2, deltaT, &
      wcool, Fbot, lcooling_general, &
      ss_const, chi_t, chi_th, chi_rho, chit_prof1, &
      chit_prof2, chi_shock, chi, iheatcond, Kgperp, Kgpara, cool_RTV, &
      tau_ss_exterior, lmultilayer, Kbot, tau_cor, TT_cor, z_cor, &
      tauheat_buffer, TTheat_buffer, zheat_buffer, dheat_buffer1, &
      heat_uniform, cool_uniform, cool_newton, lupw_ss, cool_int, cool_ext, &
      chi_hyper3, chi_hyper3_mesh, lturbulent_heat, deltaT_poleq, tdown, allp, &
      beta_glnrho_global, ladvection_entropy, lviscosity_heat, r_bcz, &
      lcalc_ss_volaverage, lcalc_ssmean, lcalc_cs2mean, lcalc_cs2mz_mean, &
      lfreeze_sint, lfreeze_sext, lhcond_global, tau_cool, &
      TTref_cool, mixinglength_flux, chiB, chi_hyper3_aniso, Ftop, xbot, &
      xtop, tau_cool2, tau_cool_ss, tau_diff, lfpres_from_pressure, &
      chit_aniso, chit_aniso_prof1, chit_aniso_prof2, &
      lchit_aniso_simplified, lconvection_gravx, &
      ltau_cool_variable, TT_powerlaw, lcalc_ssmeanxy, hcond0_kramers, &
      nkramers, xbot_aniso, xtop_aniso, entropy_floor, lprestellar_cool_iso, &
      zz1, zz2
!
!  Diagnostic variables for print.in
!  (need to be consistent with reset list below).
!
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
  integer :: idiag_yHm=0        ! DIAG_DOC: mean hydrogen ionization
  integer :: idiag_yHmax=0      ! DIAG_DOC: max of hydrogen ionization
  integer :: idiag_TTm=0        ! DIAG_DOC: $\left<T\right>$
  integer :: idiag_TTmax=0      ! DIAG_DOC: $T_{\max}$
  integer :: idiag_TTmin=0      ! DIAG_DOC: $T_{\min}$
  integer :: idiag_ssmax=0      ! DIAG_DOC: $s_{\max}$
  integer :: idiag_ssmin=0      ! DIAG_DOC: $s_{\min}$
  integer :: idiag_gTrms=0      ! DIAG_DOC: $(\nabla T)_{\rm rms}$
  integer :: idiag_gsrms=0      ! DIAG_DOC: $(\nabla s)_{\rm rms}$
  integer :: idiag_gTxgsrms=0   ! DIAG_DOC: $(\nabla T\times\nabla s)_{\rm rms}$
  integer :: idiag_fconvm=0     ! DIAG_DOC: $\left<\varrho u_z T \right>$
  integer :: idiag_TTp=0        ! DIAG_DOC:
  integer :: idiag_ssmr=0       ! DIAG_DOC:
  integer :: idiag_TTmr=0       ! DIAG_DOC:
  integer :: idiag_ufpresm=0    ! DIAG_DOC: $\left< -u/\rho\nabla p \right>$
  integer :: idiag_uduum=0
!
! xy averaged diagnostics given in xyaver.in
!
  integer :: idiag_fradz=0      ! XYAVG_DOC: $F_{\rm rad}$
  integer :: idiag_fconvz=0     ! XYAVG_DOC: $\left<\varrho u_z T \right>_{xy}$
  integer :: idiag_ssmz=0       ! XYAVG_DOC: $\left< s \right>_{xy}$
  integer :: idiag_ppmz=0       ! XYAVG_DOC: $\left< p \right>_{xy}$
  integer :: idiag_TTmz=0       ! XYAVG_DOC: $\left< T \right>_{xy}$
  integer :: idiag_TT2mz=0      ! XYAVG_DOC: $\left< T^2 \right>_{xy}$
  integer :: idiag_uxTTmz=0     ! XYAVG_DOC: $\left< u_x T \right>_{xy}$
  integer :: idiag_uyTTmz=0     ! XYAVG_DOC: $\left< u_y T \right>_{xy}$
  integer :: idiag_uzTTmz=0     ! XYAVG_DOC: $\left< u_z T \right>_{xy}$
  integer :: idiag_fradz_kramers=0 ! XYAVG_DOC: $F_{\rm rad}$ (from Kramers'
                                   ! XYAVG_DOC: opacity)
  integer :: idiag_fradz_Kprof=0 ! XYAVG_DOC: $F_{\rm rad}$ (from Kprof)
  integer :: idiag_fturbz=0     ! XYAVG_DOC: $\left<\varrho T \chi_t \nabla_z
                                ! XYAVG_DOC: s\right>_{xy}$ \quad(turbulent
                                ! XYAVG_DOC: heat flux)
  integer :: idiag_dcoolz=0     ! XYAVG_DOC: surface cooling flux
!
! xz averaged diagnostics given in xzaver.in
!
  integer :: idiag_ssmy=0       ! XZAVG_DOC: $\left< s \right>_{xz}$
  integer :: idiag_ppmy=0       ! XZAVG_DOC: $\left< p \right>_{xz}$
  integer :: idiag_TTmy=0       ! XZAVG_DOC: $\left< T \right>_{xz}$
!
! yz averaged diagnostics given in yzaver.in
!
  integer :: idiag_ssmx=0       ! YZAVG_DOC: $\left< s \right>_{yz}$
  integer :: idiag_ppmx=0       ! YZAVG_DOC: $\left< p \right>_{yz}$
  integer :: idiag_TTmx=0       ! YZAVG_DOC: $\left< T \right>_{yz}$
  integer :: idiag_uxTTmx=0     ! YZAVG_DOC: $\left< u_x T \right>_{yz}$
!
! y averaged diagnostics given in yaver.in
!
  integer :: idiag_TTmxz=0     ! YAVG_DOC: $\left< T \right>_{y}$
  integer :: idiag_ssmxz=0     ! YAVG_DOC: $\left< s \right>_{y}$
!
! z averaged diagnostics given in zaver.in
!
  integer :: idiag_TTmxy=0     ! ZAVG_DOC: $\left< T \right>_{z}$
  integer :: idiag_ssmxy=0     ! ZAVG_DOC: $\left< s \right>_{z}$
  integer :: idiag_uxTTmxy=0   ! ZAVG_DOC: $\left< u_x T \right>_{z}$
  integer :: idiag_uyTTmxy=0   ! ZAVG_DOC: $\left< u_y T \right>_{z}$
  integer :: idiag_uzTTmxy=0   ! ZAVG_DOC: $\left< u_z T \right>_{z}$
  integer :: idiag_fconvxy=0   ! ZAVG_DOC: $\left<\varrho u_x T \right>_{z}$
  integer :: idiag_fconvyxy=0  ! ZAVG_DOC: $\left<\varrho u_y T \right>_{z}$
  integer :: idiag_fconvzxy=0  ! ZAVG_DOC: $\left<\varrho u_z T \right>_{z}$
  integer :: idiag_fradxy_Kprof=0 ! ZAVG_DOC: $F_{\rm rad}$ ($z$-averaged, from Kprof)
  integer :: idiag_fradxy_kramers=0 ! ZAVG_DOC: $F_{\rm rad}$ ($z$-averaged,
                                    ! ZAVG_DOC: from Kramers' opacity)
  integer :: idiag_fturbxy=0   ! ZAVG_DOC: $\left<\varrho T \chi_t \nabla_x
                               ! ZAVG_DOC: s\right>_{z}$
  integer :: idiag_fturbrxy=0  ! ZAVG_DOC: $\left<\varrho T \chi_{ri} \nabla_i
                               ! ZAVG_DOC: s\right>_{z}$ \quad(radial part
                               ! ZAVG_DOC: of anisotropic turbulent heat flux)
  integer :: idiag_fturbthxy=0 ! ZAVG_DOC: $\left<\varrho T \chi_{\theta i}
                               ! ZAVG_DOC: \nabla_i s\right>_{z}$ \quad
                               ! ZAVG_DOC: (latitudinal part of anisotropic
                               ! ZAVG_DOC: turbulent heat flux)
  integer :: idiag_dcoolxy=0   ! ZAVG_DOC: surface cooling flux
!
  contains
!***********************************************************************
    subroutine register_entropy()
!
!  Initialise variables which should know that we solve an entropy
!  equation: iss, etc; increase nvar accordingly.
!
!  6-nov-01/wolf: coded
!
      use FArrayManager, only: farray_register_pde
      use SharedVariables, only: get_shared_variable
!
      integer :: ierr
!
      call farray_register_pde('ss',iss)
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  logical variable lpressuregradient_gas shared with hydro modules
!
      call get_shared_variable('lpressuregradient_gas',lpressuregradient_gas,ierr)
      if (ierr/=0) call fatal_error('register_entropy','lpressuregradient_gas')
!
    endsubroutine register_entropy
!***********************************************************************
    subroutine initialize_entropy(f,lstarting)
!
!  Called by run.f90 after reading parameters, but before the time loop.
!
!  21-jul-02/wolf: coded
!
      use BorderProfiles, only: request_border_driving
      use EquationOfState, only: cs0, get_soundspeed, get_cp1, &
                                 beta_glnrho_global, beta_glnrho_scaled, &
                                 mpoly, mpoly0, mpoly1, mpoly2, &
                                 select_eos_variable,gamma,gamma_m1
      use FArrayManager
      use Gravity, only: gravz, g0, compute_gravity_star
      use Initcond
      use SharedVariables, only: put_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      real, dimension (nx,3) :: glhc
      real, dimension (nx) :: hcond
      real :: beta1, cp1, beta0, TT_crit, star_cte
      integer :: i, ierr, q
      logical :: lnothing
      type (pencil_case) :: p
!
!  Check any module dependencies.
!
      if (.not. leos) then
        call fatal_error('initialize_entropy', &
            'EOS=noeos but entropy requires an EQUATION OF STATE for the fluid')
      endif
!
!  Tell the equation of state that we are here and what f variable we use.
!
      if (pretend_lnTT) then
        call select_eos_variable('lnTT',iss)
      else
        call select_eos_variable('ss',iss)
      endif
!
!  Radiative diffusion: initialize flux etc.
!
      hcond=0.0
!
!  Kbot and hcond0 are used interchangibly, so if one is
!  =impossible, set it to the other's value.
!
      if (hcond0==impossible) then
        if (Kbot==impossible) then
          hcond0=0.0
          Kbot=0.0
        else                    ! Kbot = possible
          hcond0=Kbot
        endif
      else                      ! hcond0 = possible
        if (Kbot==impossible) then
          Kbot=hcond0
        else
          call warning('initialize_entropy', &
              'You should not set Kbot and hcond0 at the same time')
        endif
      endif
!
!  hcond0 is given by mixinglength_flux in the MLT case
!  hcond1 and hcond2 follow below in the block 'if (lmultilayer) etc...'
!
      if (mixinglength_flux/=0.0) then
        Fbot=mixinglength_flux
        hcond0=-mixinglength_flux*(mpoly0+1.)*gamma_m1*gamma_inv/gravz
        hcond1 = (mpoly1+1.)/(mpoly0+1.)
        hcond2 = (mpoly2+1.)/(mpoly0+1.)
        Kbot=hcond0
        lmultilayer=.true.  ! just to be sure...
        if (lroot) print*, &
            'initialize_entropy: hcond0 given by mixinglength_flux=', &
            hcond0, Fbot
      endif
!
!  Freeze entropy.
!
      if (lfreeze_sint) lfreeze_varint(iss) = .true.
      if (lfreeze_sext) lfreeze_varext(iss) = .true.
!
!  Make sure the top boundary condition for temperature (if cT is used)
!  knows about the cooling function or vice versa (cs2cool will take over
!  if /=0).
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
!  Settings for fluxes.
!
        if (lmultilayer) then
!
!  Calculate hcond1,hcond2 if they have not been set in run.in.
!
          if (hcond1==impossible) hcond1 = (mpoly1+1.)/(mpoly0+1.)
          if (hcond2==impossible) hcond2 = (mpoly2+1.)/(mpoly0+1.)
!
!  Calculate Fbot if it has not been set in run.in.
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
          if (hcond0*hcond1/=0.0) then
            FbotKbot=Fbot/(hcond0*hcond1)
            hcondzbot=hcond0*hcond1
          else
            FbotKbot=0.0
          endif
!
!  Calculate Ftop if it has not been set in run.in.
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
          if (hcond0*hcond2/=0.0) then
            FtopKtop=Ftop/(hcond0*hcond2)
            hcondztop=hcond0*hcond2
          else
            FtopKtop=0.0
          endif
!
        else  ! lmultilayer=.false.
!
!  NOTE: in future we should define chiz=chi(z) or Kz=K(z) here.
!  Calculate hcond and FbotKbot=Fbot/K
!  (K=hcond is radiative conductivity).
!
!  Calculate Fbot if it has not been set in run.in.
!
          if (Fbot==impossible) then
            if (bcz1(iss)=='c1') then
              Fbot=-gamma/(gamma-1)*hcond0*gravz/(mpoly+1)
              if (lroot) print*, &
                  'initialize_entropy: Calculated Fbot = ', Fbot
!
              Kbot=gamma_m1/gamma*(mpoly+1.)*Fbot
              FbotKbot=gamma/gamma_m1/(mpoly+1.)
              if (lroot) print*, 'initialize_entropy: Calculated Fbot,Kbot=', &
                  Fbot, Kbot
            ! else
            !! Don't need Fbot in this case (?)
            !  Fbot=-gamma/(gamma-1)*hcond0*gravz/(mpoly+1)
            !  if (lroot) print*, &
            !       'initialize_entropy: Calculated Fbot = ', Fbot
            endif
          else
!
!  Define hcond0 from the given value of Fbot.
!
            if (bcz1(iss)=='c1') then
              hcond0=-gamma_m1/gamma*(mpoly0+1.)*Fbot/gravz
              Kbot=hcond0
              FbotKbot=-gravz*gamma/gamma_m1/(mpoly0+1.)
              if (lroot) print*, &
                   'initialize_entropy: Calculated hcond0 = ', hcond0
            endif
          endif
!
!  Calculate Ftop if it has not been set in run.in.
!
          if (Ftop==impossible) then
            if (bcz2(iss)=='c1') then
              Ftop=-gamma/(gamma-1)*hcond0*gravz/(mpoly+1)
              if (lroot) print*, &
                      'initialize_entropy: Calculated Ftop = ', Ftop
              Ktop=gamma_m1/gamma*(mpoly+1.)*Ftop
              FtopKtop=gamma/gamma_m1/(mpoly+1.)
              if (lroot) print*,'initialize_entropy: Ftop,Ktop=',Ftop,Ktop
            ! else
            !! Don't need Ftop in this case (?)
            !  Ftop=0.
            endif
          endif
!
        endif
! FbotKbot is required here for boundcond
      elseif (lgravx) then
        if (coord_system=='spherical'.or.lconvection_gravx.and.(.not.lhcond_global)) then
             FbotKbot=Fbot/hcond0
             hcondxbot = hcond0
        endif
      endif
!
!  Make sure all relevant parameters are set for spherical shell problems.
!
      select case (initss(1))
        case ('geo-kws','geo-benchmark','shell_layers')
          if (lroot) then
            print*, 'initialize_entropy: '// &
                'set boundary temperatures for spherical shell problem'
          endif
!
!  Calculate temperature gradient from polytropic index.
!
          call get_cp1(cp1)
          beta1=cp1*g0/(mpoly+1)*gamma/gamma_m1
!
!  Temperatures at shell boundaries.
!
          if (initss(1) == 'shell_layers') then
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
            lmultilayer=.false.    ! to ensure that hcond=cte
            if (T0/=0.0) then
              TT_ext=T0            ! T0 defined in start.in for geodynamo
            else
              TT_ext=cs20/gamma_m1 ! T0 defined from cs20
            endif
            if (coord_system=='spherical')then
              r_ext=x(l2)
              r_int=x(l1)
            else
            endif
            TT_int=TT_ext*(1.+beta1*(r_ext/r_int-1.))
          endif
!
!  Set up cooling parameters for spherical shell in terms of sound speeds
!
          call get_soundspeed(TT_ext,cs2_ext)
          call get_soundspeed(TT_int,cs2_int)
          cs2cool=cs2_ext
          if (lroot) then
            print*,'initialize_entropy: g0,mpoly,beta1',g0,mpoly,beta1
            print*,'initialize_entropy: TT_int, TT_ext=',TT_int,TT_ext
            print*,'initialize_entropy: cs2_ext, cs2_int=',cs2_ext, cs2_int
          endif
!
        case ('star_heat')
          if (hcond1==impossible) hcond1=(mpoly1+1.)/(mpoly0+1.)
          if (hcond2==impossible) hcond2=(mpoly2+1.)/(mpoly0+1.)
          if (lroot) print*,'initialize_entropy: set cs2cool=cs20'
          cs2cool=cs20
          cs2_ext=cs20
          if (rcool==0.) rcool=r_ext
!
!  Compute the gravity profile inside the star.
!
          star_cte=(mpoly0+1.)/hcond0*(1.-gamma_inv)
          call compute_gravity_star(f, wheat, luminosity, star_cte)
!
        case ('cylind_layers')
          if (bcx1(iss)=='c1') then
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
!  gradient dlnrho/dr = beta/H.
!
      if (maxval(abs(beta_glnrho_global))/=0.0) then
        beta_glnrho_scaled=beta_glnrho_global*Omega/cs0
        if (lroot) print*, 'initialize_entropy: Global density gradient '// &
            'with beta_glnrho_global=', beta_glnrho_global
      endif
!
!  Turn off pressure gradient term and advection for 0-D runs.
!
      if (nwgrid==1) then
        lpressuregradient_gas=.false.
        ladvection_entropy=.false.
        print*, 'initialize_entropy: 0-D run, turned off pressure gradient term'
        print*, 'initialize_entropy: 0-D run, turned off advection of entropy'
      endif
!
!  Turn off advection of entropy when fargo is used
!
      if (lfargo_advection) then
        ladvection_entropy=.false.
        if (lroot) print*,&
             'initialize_entropy: fargo used, turned off advection of entropy'
      endif
!
!  Possible to calculate pressure gradient directly from the pressure, in which
!  case we must open an auxiliary slot in f to store the pressure. This method
!  is not compatible with non-blocking communication of ghost zones, so we turn
!  this behaviour off.
!
      if (lfpres_from_pressure) then
        call farray_register_auxiliary('ppaux',ippaux)
        test_nonblocking=.true.
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
      lheatc_kramers=.false.
      lheatc_corona=.false.
      lheatc_shock=.false.
      lheatc_hyper3ss=.false.
      lheatc_hyper3ss_polar=.false.
      lheatc_hyper3ss_mesh=.false.
      lheatc_hyper3ss_aniso=.false.
!
      lnothing=.false.
!
!  Select which radiative heating we are using.
!
      if (lroot) print*,'initialize_entropy: nheatc_max,iheatcond=',nheatc_max,iheatcond(1:nheatc_max)
      do i=1,nheatc_max
        select case (iheatcond(i))
        case ('K-profile')
          lheatc_Kprof=.true.
          if (lroot) print*, 'heat conduction: K-profile'
        case ('K-const')
          lheatc_Kconst=.true.
          if (lroot) print*, 'heat conduction: K=cte'
        case ('chi-const')
          lheatc_chiconst=.true.
          if (lroot) print*, 'heat conduction: constant chi'
        case ('chi-therm')
          lheatc_chitherm=.true.
          if (lroot) print*, 'heat conduction: chi scaled with c_s'
        case ('sqrtrhochi-const')
          lheatc_sqrtrhochiconst=.true.
          if (lroot) print*, 'heat conduction: constant sqrt(rho)*chi'
        case ('tensor-diffusion')
          lheatc_tensordiffusion=.true.
          if (lroot) print*, 'heat conduction: tensor diffusion'
        case ('spitzer')
          lheatc_spitzer=.true.
          if (lroot) print*, 'heat conduction: spitzer'
        case ('hubeny')
          lheatc_hubeny=.true.
          if (lroot) print*, 'heat conduction: hubeny'
        case ('kramers')
          lheatc_kramers=.true.
          if (lroot) print*, 'heat conduction: kramers'
        case ('corona')
          lheatc_corona=.true.
          if (lroot) print*, 'heat conduction: corona'
        case ('shock')
          lheatc_shock=.true.
          if (lroot) print*, 'heat conduction: shock'
        case ('hyper3_ss','hyper3')
          lheatc_hyper3ss=.true.
          if (lroot) print*, 'heat conduction: hyperdiffusivity of ss'
       case ('hyper3_aniso','hyper3-aniso')
          if (lroot) print*, 'heat conduction: anisotropic '//&
               'hyperdiffusivity of ss'
          lheatc_hyper3ss_aniso=.true.
        case ('hyper3_cyl','hyper3-cyl','hyper3-sph','hyper3_sph')
          lheatc_hyper3ss_polar=.true.
          if (lroot) print*, 'heat conduction: hyperdiffusivity of ss'
        case ('hyper3-mesh','hyper3_mesh')
          lheatc_hyper3ss_mesh=.true.
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
!
!  A word of warning...
!
      if (lheatc_Kprof .and. hcond0==0.0) then
        call warning('initialize_entropy', 'hcond0 is zero!')
      endif
      if (lheatc_chiconst .and. (chi==0.0 .and. chi_t==0.0)) then
        call warning('initialize_entropy','chi and chi_t are zero!')
      endif
      if (lheatc_chitherm .and. (chi_th==0.0 .and. chi_t==0.0)) then
        call warning('initialize_entropy','chi_th and chi_t are zero!')
      endif
      if (lheatc_sqrtrhochiconst .and. (chi_rho==0.0 .and. chi_t==0.0)) then
        call warning('initialize_entropy','chi_rho and chi_t are zero!')
      endif
      if (all(iheatcond=='nothing') .and. hcond0/=0.0) then
        call warning('initialize_entropy', &
            'No heat conduction, but hcond0 /= 0')
      endif
      if (lheatc_Kconst .and. Kbot==0.0) then
        call warning('initialize_entropy','Kbot is zero!')
      endif
      if ((lheatc_spitzer.or.lheatc_corona) .and. (Kgpara==0.0 .or. Kgperp==0.0) ) then
        call warning('initialize_entropy','Kgperp or Kgpara is zero!')
      endif
      if (lheatc_kramers .and. hcond0_kramers==0.0) then
        call warning('initialize_entropy','hcond0_kramers is zero!')
      endif
      if (lheatc_hyper3ss .and. chi_hyper3==0.0) &
           call warning('initialize_entropy','chi_hyper3 is zero!')
      if (lheatc_hyper3ss_polar .and. chi_hyper3==0.0) &
           call warning('initialize_entropy','chi_hyper3 is zero!')
      if (lheatc_hyper3ss_mesh .and. chi_hyper3_mesh==0.0) &
           call warning('initialize_entropy','chi_hyper3_mesh is zero!')
      if ( (lheatc_hyper3ss_aniso) .and.  &
           ((chi_hyper3_aniso(1)==0. .and. nxgrid/=1 ).or. &
            (chi_hyper3_aniso(2)==0. .and. nygrid/=1 ).or. &
            (chi_hyper3_aniso(3)==0. .and. nzgrid/=1 )) ) &
           call fatal_error('initialize_entropy', &
           'A diffusivity coefficient of chi_hyper3 is zero!')
      if (lheatc_shock .and. chi_shock==0.0) then
        call warning('initialize_entropy','chi_shock is zero!')
      endif
!
!  Heat conduction calculated globally: if lhcond_global=.true., we
!  compute hcond and glhc and put the results in global arrays.
!
      if (lhcond_global) then
        call farray_register_auxiliary('hcond',iglobal_hcond)
        call farray_register_auxiliary('glhc',iglobal_glhc,vector=3)
        if (coord_system=='spherical' .or. lconvection_gravx) then
          call read_hcond(hcond,glhc)
          do q=n1,n2; do m=m1,m2
            f(l1:l2,m,q,iglobal_hcond)=hcond
            f(l1:l2,m,q,iglobal_glhc:iglobal_glhc+2)=glhc
          enddo; enddo
          FbotKbot=Fbot/hcond(1)
          FtopKtop=Ftop/hcond(nx)
          hcondxbot=hcond(1)
          hcondxtop=hcond(nx)
        else
          do n=n1,n2; do m=m1,m2
            if (lgravz) then
              p%z_mn=spread(z(n),1,nx)
            else
              p%r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
            endif
            call heatcond(hcond,p)
            call gradloghcond(glhc,p)
            f(l1:l2,m,n,iglobal_hcond)=hcond
            f(l1:l2,m,n,iglobal_glhc:iglobal_glhc+2)=glhc
          enddo; enddo
          hcondxbot=hcond(1)
          hcondxtop=hcond(nx)
        endif
      endif
!
!  Tell the BorderProfiles module if we intend to use border driving, so
!  that the module can request the right pencils.
!
      if (borderss/='nothing') call request_border_driving(borderss)
!
!  Shared variables.
!
      call put_shared_variable('hcond0',hcond0,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting hcond0')
      call put_shared_variable('hcond1',hcond1,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting hcond1')
      call put_shared_variable('Kbot',Kbot,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
           'there was a problem when putting Kbot')
      call put_shared_variable('hcondxbot',hcondxbot,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting hcondxbot')
      call put_shared_variable('hcondxtop',hcondxtop,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting hcondxtop')
      call put_shared_variable('hcondzbot',hcondzbot,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting hcondzbot')
      call put_shared_variable('hcondztop',hcondztop,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting hcondztop')
      call put_shared_variable('Fbot',Fbot,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting Fbot')
      call put_shared_variable('Ftop',Ftop,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting Ftop')
      call put_shared_variable('FbotKbot',FbotKbot,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting FbotKbot')
      call put_shared_variable('FtopKtop',FtopKtop,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting FtopKtop')
      call put_shared_variable('chi',chi,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting chi')
      call put_shared_variable('chi_t',chi_t,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting chi_t')
!      call put_shared_variable('chi_th',chi_th,ierr)
!      if (ierr/=0) call fatal_error('initialize_entropy', &
!          'there was a problem when putting chi_th')
!      call put_shared_variable('chi_rho',chi_rho,ierr)
!      if (ierr/=0) call fatal_error('initialize_entropy', &
!          'there was a problem when putting chi_rho')
      call put_shared_variable('lmultilayer',lmultilayer,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting lmultilayer')
      call put_shared_variable('lheatc_chiconst',lheatc_chiconst,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting lcalc_heatcond_constchi')
!      call put_shared_variable('lheatc_chitherm',lheatc_chitherm,ierr)
!      if (ierr/=0) call fatal_error('initialize_entropy', &
!          'there was a problem when putting lcalc_heatcond_chitherm')
!      call put_shared_variable('lheatc_sqrtrhochiconst',lheatc_sqrtrhochiconst,ierr)
!      if (ierr/=0) call fatal_error('initialize_entropy', &
!          'there was a problem when putting lcalc_heatcond_sqrtrhochiconst')
      call put_shared_variable('lviscosity_heat',lviscosity_heat,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting lviscosity_heat')
      call put_shared_variable('hcond0_kramers',hcond0_kramers,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting hcond0_kramers')
      call put_shared_variable('nkramers',nkramers,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting nkramers')
      call put_shared_variable('lheatc_kramers',lheatc_kramers,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting lheatc_kramers')
!
      call keep_compiler_quiet(lstarting)
!
      endsubroutine initialize_entropy
!***********************************************************************
    subroutine read_entropy_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=entropy_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=entropy_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_entropy_init_pars
!***********************************************************************
    subroutine write_entropy_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=entropy_init_pars)
!
    endsubroutine write_entropy_init_pars
!***********************************************************************
    subroutine read_entropy_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=entropy_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=entropy_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_entropy_run_pars
!***********************************************************************
    subroutine write_entropy_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=entropy_run_pars)
!
    endsubroutine write_entropy_run_pars
!***********************************************************************
    subroutine init_ss(f)
!
!  Initial condition for entropy.
!
!  07-nov-2001/wolf: coded
!  24-nov-2002/tony: renamed for consistancy (i.e. init_[variable name])
!
      use EquationOfState,  only: isothtop, &
                                mpoly0, mpoly1, mpoly2, cs2cool, cs0, &
                                rho0, lnrho0, isothermal_entropy, &
                                isothermal_lnrho_ss, eoscalc, ilnrho_pp, &
                                eosperturb
      use General, only: itoa
      use Gravity
      use Initcond
      use InitialCondition, only: initial_condition_ss
      use Sub
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
!  Select different initial conditions.
!
        select case (initss(j))
!
          case ('zero', '0'); f(:,:,:,iss) = 0.
          case ('const_ss'); f(:,:,:,iss)=f(:,:,:,iss)+ss_const
          case ('gaussian-noise'); call gaunoise(ampl_ss,f,iss,iss)
          case ('blob')
            call blob(ampl_ss,f,iss,radius_ss,center1_x,center1_y,center1_z)
          case ('blob_radeq')
            call blob_radeq(ampl_ss,f,iss,radius_ss,center1_x,center1_y,center1_z)
          case ('isothermal'); call isothermal_entropy(f,T0)
          case ('isothermal_lnrho_ss')
            print*, 'init_ss: Isothermal density and entropy stratification'
            call isothermal_lnrho_ss(f,T0,rho0)
          case ('hydrostatic-isentropic')
            call hydrostatic_isentropic(f,lnrho_bot,ss_const)
          case ('wave')
            do n=n1,n2; do m=m1,m2
              f(l1:l2,m,n,iss)=f(l1:l2,m,n,iss)+ss_const+ampl_ss*sin(kx_ss*x(l1:l2)+pi)
            enddo; enddo
          case('Ferriere'); call ferriere(f)
          case('Ferriere-hs'); call ferriere_hs(f,rho0hs)
          case('thermal-hs')
            call thermal_hs_equilibrium_ism(f)
          case('Galactic-hs'); call galactic_hs(f,rho0hs,cs0hs,H0hs)
          case('xjump'); call jump(f,iss,ss_left,ss_right,widthss,'x')
          case('yjump'); call jump(f,iss,ss_left,ss_right,widthss,'y')
          case('zjump'); call jump(f,iss,ss_left,ss_right,widthss,'z')
          case('sinxsinz'); call sinxsinz(ampl_ss,f,iss,kx_ss,ky_ss,kz_ss)
          case('hor-fluxtube')
            call htube(ampl_ss,f,iss,iss,radius_ss,epsilon_ss, &
                center1_x,center1_z)
          case ('hor-tube')
            call htube2(ampl_ss,f,iss,iss,radius_ss,epsilon_ss)
          case ('mixinglength')
             call mixinglength(mixinglength_flux,f)
             hcond0=-mixinglength_flux*(mpoly0+1.)*gamma_m1*gamma_inv/gravz
             print*,'init_ss : Fbot, hcond0=', Fbot, hcond0
          case ('sedov')
            if (lroot) print*,'init_ss: sedov - thermal background with gaussian energy burst'
            call blob(thermal_peak,f,iss,radius_ss, &
                center1_x,center1_y,center1_z)
          case ('sedov-dual')
            if (lroot) print*,'init_ss: sedov - thermal background with gaussian energy burst'
            call blob(thermal_peak,f,iss,radius_ss, &
                center1_x,center1_y,center1_z)
            call blob(thermal_peak,f,iss,radius_ss, &
                center2_x,center2_y,center2_z)
          case ('shock2d')
            call shock2d(f)
          case ('isobaric')
!
!  ss = - ln(rho/rho0)
!
            if (lroot) print*,'init_ss: isobaric stratification'
            if (pp_const==0.) then
              f(:,:,:,iss) = -(f(:,:,:,ilnrho)-lnrho0)
            else
              pp=pp_const
              do n=n1,n2; do m=m1,m2
                lnrho=f(l1:l2,m,n,ilnrho)
                call eoscalc(ilnrho_pp,lnrho,pp,ss=ss)
                f(l1:l2,m,n,iss)=ss
              enddo; enddo
            endif
          case ('isentropic', '1')
!
!  ss = const.
!
            if (lroot) print*,'init_ss: isentropic stratification'
            f(:,:,:,iss)=0.
            if (ampl_ss/=0.) then
              print*, 'init_ss: put bubble: radius_ss,ampl_ss=', &
                  radius_ss, ampl_ss
              do n=n1,n2; do m=m1,m2
                tmp=x(l1:l2)**2+y(m)**2+z(n)**2
                f(l1:l2,m,n,iss)=f(l1:l2,m,n,iss)+ &
                    ampl_ss*exp(-tmp/max(radius_ss**2-tmp,1e-20))
              enddo; enddo
            endif
          case ('linprof', '2')
!
!  Linear profile of ss, centered around ss=0.
!
            if (lroot) print*,'init_ss: linear entropy profile'
            do n=n1,n2; do m=m1,m2
              f(l1:l2,m,n,iss) = grads0*z(n)
            enddo; enddo
          case ('isentropic-star')
!
!  Isentropic/isothermal hydrostatic sphere:
!    ss  = 0       for r<R,
!    cs2 = const   for r>R
!
!  Only makes sense if both initlnrho=initss='isentropic-star'.
!
            if (.not. ldensity) &
                call fatal_error('isentropic-star','requires density.f90')
            if (lgravr) then
              if (lroot) print*, &
                  'init_lnrho: isentropic star with isothermal atmosphere'
              do n=n1,n2; do m=m1,m2
                r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
                call potential(POT=pot,POT0=pot0,RMN=r_mn) ! gravity potential
!
!  rho0, cs0,pot0 are the values in the centre.
!
                if (gamma/=1.0) then
!
!  Note:
!  (a) `where' is expensive, but this is only done at
!      initialization.
!  (b) Comparing pot with pot_ext instead of r with r_ext will
!      only work if grav_r<=0 everywhere -- but that seems
!      reasonable.
!
                  call potential(R=r_ext,POT=pot_ext) ! get pot_ext=pot(r_ext)
                  cs2_ext   = cs20*(1 - gamma_m1*(pot_ext-pot0)/cs20)
!
!  Make sure init_lnrho (or start.in) has already set cs2cool.
!
                  if (cs2cool==0.0) call fatal_error('init_ss', &
                      'inconsistency - cs2cool ca not be 0')
                  ss_ext = 0.0 + log(cs2cool/cs2_ext)
                  where (pot <= pot_ext) ! isentropic for r<r_ext
                    f(l1:l2,m,n,iss) = 0.0
                  elsewhere              ! isothermal for r>r_ext
                    f(l1:l2,m,n,iss) = ss_ext + gamma_m1*(pot-pot_ext)/cs2cool
                  endwhere
                else         ! gamma=1 --> simply isothermal (I guess [wd])
                  ! [NB: Never tested this..]
                  f(l1:l2,m,n,iss) = -gamma_m1/gamma*(f(l1:l2,m,n,ilnrho)-lnrho0)
                endif
              enddo; enddo
            endif
          case ('piecew-poly', '4')
!
!  Piecewise polytropic convection setup.
!  cs0, rho0 and ss0=0 refer to height z=zref
!
            if (lroot) print*, &
                'init_ss: piecewise polytropic vertical stratification (ss)'
!
            cs2int = cs0**2
            ss0 = 0.              ! reference value ss0 is zero
            ssint = ss0
            f(:,:,:,iss) = 0.    ! just in case
!  Top layer.
            call polytropic_ss_z(f,mpoly2,zref,z2,z0+2*Lz, &
                                 isothtop,cs2int,ssint)
!  Unstable layer.
            call polytropic_ss_z(f,mpoly0,z2,z1,z2,0,cs2int,ssint)
!  Stable layer.
            call polytropic_ss_z(f,mpoly1,z1,z0,z1,0,cs2int,ssint)
          case ('piecew-disc', '41')
!
!  Piecewise polytropic convective disc.
!  cs0, rho0 and ss0=0 refer to height z=zref.
!
            if (lroot) print*, 'init_ss: piecewise polytropic disc'
!
            ztop = xyz0(3)+Lxyz(3)
            cs2int = cs0**2
            ss0 = 0.              ! reference value ss0 is zero
            ssint = ss0
            f(:,:,:,iss) = 0.    ! just in case
!  Bottom (middle) layer.
            call polytropic_ss_disc(f,mpoly1,zref,z1,z1, &
                                 0,cs2int,ssint)
!  Unstable layer.
            call polytropic_ss_disc(f,mpoly0,z1,z2,z2,0,cs2int,ssint)
!  Stable layer (top).
            call polytropic_ss_disc(f,mpoly2,z2,ztop,ztop,&
                                 isothtop,cs2int,ssint)
          case ('polytropic', '5')
!
!  Polytropic stratification.
!  cs0, rho0 and ss0=0 refer to height z=zref.
!
            if (lroot) print*,'init_ss: polytropic vertical stratification'
!
            cs20 = cs0**2
            ss0 = 0.             ! reference value ss0 is zero
            f(:,:,:,iss) = ss0   ! just in case
            cs2int = cs20
            ssint = ss0
!  Only one layer.
            call polytropic_ss_z(f,mpoly0,zref,z0,z0+2*Lz,0,cs2int,ssint)
!  Reset mpoly1, mpoly2 (unused) to make IDL routine `thermo.pro' work.
            mpoly1 = mpoly0
            mpoly2 = mpoly0
          case ('geo-kws')
!
!  Radial temperature profiles for spherical shell problem.
!
            if (lroot) print*,'init_ss: kws temperature in spherical shell'
            call shell_ss(f)
          case ('geo-benchmark')
!
!  Radial temperature profiles for spherical shell problem.
!
            if (lroot) print*, 'init_ss: '// &
                'benchmark temperature in spherical shell'
            call shell_ss(f)
          case ('shell_layers')
!
!  Radial temperature profiles for spherical shell problem.
!
            call information('init_ss', &
                'two polytropic layers in a spherical shell')
            call shell_ss_layers(f)
          case ('star_heat')
!
!  Radial temperature profiles for spherical shell problem.
!
            call information('init_ss', &
                'two polytropic layers with a central heating')
            call star_heat(f)
          case ('cylind_layers')
            call cylind_layers(f)
          case ('polytropic_simple')
!
!  Vertical temperature profiles for convective layer problem.
!
            call layer_ss(f)
          case ('blob_hs')
            print*, 'init_ss: put blob in hydrostatic equilibrium: '// &
                'radius_ss,ampl_ss=', radius_ss, ampl_ss
            call blob(ampl_ss,f,iss,radius_ss,center1_x,center1_y,center1_z)
            call blob(-ampl_ss,f,ilnrho,radius_ss,center1_x,center1_y,center1_z)
          case ('single_polytrope')
            call single_polytrope(f)
          case default
!
!  Catch unknown values.
!
            write(unit=errormsg,fmt=*) 'No such value for initss(' &
                //trim(iinit_str)//'): ',trim(initss(j))
            call fatal_error('init_ss',errormsg)
        endselect
!
        if (lroot) print*, 'init_ss: initss('//trim(iinit_str)//') = ', &
            trim(initss(j))
!
      enddo
!
      if (lnothing.and.lroot) print*, 'init_ss: nothing'
!
!  Add perturbation(s).
!
      select case (pertss)
!
      case ('zero', '0')
!
!  Do not perturb.
!
      case ('hexagonal', '1')
!
!  Hexagonal perturbation.
!
        if (lroot) print*, 'init_ss: adding hexagonal perturbation to ss'
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,iss) = f(l1:l2,m,n,iss) &
              + ampl_ss*(2*cos(sqrt(3.)*0.5*khor_ss*x(l1:l2)) &
                          *cos(0.5*khor_ss*y(m)) &
                         + cos(khor_ss*y(m))) * cos(pi*z(n))
        enddo; enddo
      case default
!
!  Catch unknown values.
!
        write (unit=errormsg,fmt=*) 'No such value for pertss:', pertss
        call fatal_error('init_ss',errormsg)
!
      endselect
!
!  Interface fow user's own initial condition.
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
    endsubroutine init_ss
!***********************************************************************
    subroutine blob_radeq(ampl,f,i,radius,xblob,yblob,zblob)
!
!  Add blob-like perturbation in radiative pressure equilibrium.
!
!   6-may-07/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, optional :: xblob,yblob,zblob
      real :: ampl,radius,x01,y01,z01
      integer :: i
!
!  Single blob.
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
      use Gravity, only: gravz
      use Sub, only: step
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
!  Smoothly blend the old value (above zblend) and the new one (below
!  zblend) for the two regions.
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
      use Gravity, only: gravz, nu_epicycle
      use Sub, only: step
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mz) :: stp
      real :: tmp,mpoly,zint,zbot,zblend,beta1,cs2int,ssint, nu_epicycle2
      integer :: isoth
!
!  Warning: beta1 is here not dT/dz, but dcs2/dz = (gamma-1)c_p dT/dz.
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
!  Smoothly blend the old value (above zblend) and the new one (below
!  zblend) for the two regions.
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
!  20-feb-04/tobi: coded
!
      use EquationOfState, only: eoscalc, ilnrho_ss
      use Gravity, only: gravz
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, intent(in) :: lnrho_bot,ss_const
      real :: cs2_,lnrho,lnrho_m
!
      if (.not. lgravz) then
        call fatal_error('hydrostatic_isentropic', &
            'Currently only works for vertical gravity field')
      endif
!
!  In case this processor is not located at the very bottom perform
!  integration through lower lying processors.
!
      lnrho=lnrho_bot
      do n=1,nz*ipz
        call eoscalc(ilnrho_ss,lnrho,ss_const,cs2=cs2_)
        lnrho_m=lnrho+dz*gravz/cs2_/2
        call eoscalc(ilnrho_ss,lnrho_m,ss_const,cs2=cs2_)
        lnrho=lnrho+dz*gravz/cs2_
      enddo
!
!  Do the integration on this processor
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
!  Entropy is simply constant.
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
      use EquationOfState, only: gamma_m1, rho0, lnrho0, &
                                 cs20, cs2top, eoscalc, ilnrho_lnTT
      use General, only: safe_character_assign
      use Gravity, only: z1
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nzgrid) :: tempm,lnrhom
      real :: zm,ztop,mixinglength_flux,lnrho,ss,lnTT
      real :: zbot,rbot,rt_old,rt_new,rb_old,rb_new,crit, &
              rhotop,rhobot
      integer :: iter
      character (len=120) :: wfile
!
      if (headtt) print*,'init_ss : mixinglength stratification'
      if (.not.lgravz) then
        call fatal_error('mixinglength','works only for vertical gravity')
      endif
!
!  Do the calculation on all processors, and then put the relevant piece
!  into the f array.
!  Choose value zbot where rhobot should be applied and give two first
!  estimates for rhotop.
!
      ztop=xyz0(3)+Lxyz(3)
      zbot=z1
      rbot=1.
      rt_old=.1*rbot
      rt_new=.12*rbot
!
!  Need to iterate for rhobot=1.
!  Produce first estimate.
!
      rhotop=rt_old
      cs2top=cs20  ! just to be sure...
      call strat_MLT (rhotop,mixinglength_flux,lnrhom,tempm,rhobot)
      rb_old=rhobot
!
!  Next estimate.
!
      rhotop=rt_new
      call strat_MLT (rhotop,mixinglength_flux,lnrhom,tempm,rhobot)
      rb_new=rhobot
!
      do 10 iter=1,10
!
!  New estimate.
!
        rhotop=rt_old+(rt_new-rt_old)/(rb_new-rb_old)*(rbot-rb_old)
!
!  Check convergence.
!
        crit=abs(rhotop/rt_new-1.)
        if (crit<=1e-4) goto 20
!
        call strat_MLT (rhotop,mixinglength_flux,lnrhom,tempm,rhobot)
!
!  Update new estimates.
!
        rt_old=rt_new
        rb_old=rb_new
        rt_new=rhotop
        rb_new=rhobot
   10 continue
   20 if (lfirst_proc_z) print*,'- iteration completed: rhotop,crit=',rhotop,crit
!
!  Redefine rho0 and lnrho0 as we don't have rho0=1 at the top.
!  (important for eoscalc!)
!  Put density and entropy into f-array.
!  Write the initial stratification in data/proc*/zprof_stratMLT.dat.
!
      rho0=rhotop
      lnrho0=log(rhotop)
      print*,'new rho0 and lnrho0=',rho0,lnrho0
!
      call safe_character_assign(wfile,trim(directory)//'/zprof_stratMLT.dat')
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
!  a spherical shell.
!
!  20-oct-03/dave -- coded
!  21-aug-08/dhruba: added spherical coordinates
!
      use EquationOfState, only: eoscalc, ilnrho_lnTT, mpoly, get_cp1
      use Gravity, only: g0
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
!  Set initial condition.
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
        call fatal_error('shell_ss','not valid in cylindrical coordinates')
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
!  Compute perturbation to initial temperature profile.
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
!  Initial state entropy profile used for initss='polytropic_simple',
!  for `conv_slab' style runs, with a layer of polytropic gas in [z0,z1].
!  generalised for cp/=1.
!
      use EquationOfState, only: eoscalc, ilnrho_lnTT, mpoly, get_cp1
      use Gravity, only: gravz, zinfty
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
!  Set initial condition (first in terms of TT, and then in terms of ss).
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
!  Density profile from K. Ferriere, ApJ 497, 759, 1998,
!   eqns (6), (7), (9), (13), (14) [with molecular H, n_m, neglected]
!   at solar radius.  (for interstellar runs)
!  Entropy is set via pressure, assuming a constant T for each gas component
!   (cf. eqn 15) and crudely compensating for non-thermal sources.
!  [An alternative treatment of entropy, based on hydrostatic equilibrium,
!   might be preferable. This was implemented in serial (in e.g. r1.59)
!   but abandoned as overcomplicated to adapt for nprocz /= 0.]
!
!  AJ: PLEASE IDENTIFY AUTHOR
!
      use EquationOfState, only: eoscalc, ilnrho_pp, eosperturb
      use Mpicomm, only: mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(nx) :: absz
      double precision, dimension(nx) :: n_c,n_w,n_i,n_h
!  T in K, k_B s.t. pp is in code units ( = 9.59e-15 erg/cm/s^2)
!  (i.e. k_B = 1.381e-16 (erg/K) / 9.59e-15 (erg/cm/s^2) )
      real, parameter :: T_c_cgs=500.0,T_w_cgs=8.0e3,T_i_cgs=8.0e3,T_h_cgs=1.0e6
      real :: T_c,T_w,T_i,T_h
      real, dimension(nx) :: rho,pp,lnrho
      real, dimension(1) :: fmpi1
      real :: kpc
      double precision ::  rhoscale
!
      if (lroot) print*,'ferriere: Ferriere density and entropy profile'
!
!  First define reference values of pp, cs2, at midplane.
!  Pressure is set to 6 times thermal pressure, this factor roughly
!  allowing for other sources, as modelled by Ferriere.
!
      kpc = 3.086D21 / unit_length
      rhoscale = 1.36 * m_p * unit_length**3
      print*, 'ferrier: kpc, rhoscale =', kpc, rhoscale
      T_c=T_c_cgs/unit_temperature
      T_w=T_w_cgs/unit_temperature
      T_i=T_i_cgs/unit_temperature
      T_h=T_h_cgs/unit_temperature
!
      do n=n1,n2            ! nb: don't need to set ghost-zones here
      absz=abs(z(n))
      do m=m1,m2
!
!  Cold gas profile n_c (eq 6).
!
        n_c=0.340*(0.859*exp(-(z(n)/(0.127*kpc))**2) +         &
                   0.047*exp(-(z(n)/(0.318*kpc))**2) +         &
                   0.094*exp(-absz/(0.403*kpc)))
!
!  Warm gas profile n_w (eq 7).
!
        n_w=0.226*(0.456*exp(-(z(n)/(0.127*kpc))**2) +  &
                   0.403*exp(-(z(n)/(0.318*kpc))**2) +  &
                   0.141*exp(-absz/(0.403*kpc)))
!
!  Ionized gas profile n_i (eq 9).
!
        n_i=0.0237*exp(-absz/kpc) + 0.0013* exp(-absz/(0.150*kpc))
!
!  Hot gas profile n_h (eq 13).
!
        n_h=0.00048*exp(-absz/(1.5*kpc))
!
!  Normalised s.t. rho0 gives mid-plane density directly (in 10^-24 g/cm^3).
!
        rho=real((n_c+n_w+n_i+n_h)*rhoscale)
        lnrho=log(rho)
        f(l1:l2,m,n,ilnrho)=lnrho
!
!  Define entropy via pressure, assuming fixed T for each component.
!
        if (lentropy) then
!
!  Thermal pressure (eq 15).
!
          pp=real(k_B*unit_length**3 * &
             (1.09*n_c*T_c + 1.09*n_w*T_w + 2.09*n_i*T_i + 2.27*n_h*T_h))
!
          call eosperturb(f,nx,pp=pp)
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
      if (lroot) print*, 'ferriere: cs2bot=', cs2bot, ' cs2top=', cs2top
!
    endsubroutine ferriere
!***********************************************************************
    subroutine ferriere_hs(f,rho0hs)
!
!  Density and isothermal entropy profile in hydrostatic equilibrium
!  with the Ferriere profile set in gravity_simple.f90.
!  Use gravz_profile='Ferriere'(gravity) and initlnrho='Galactic-hs'
!  both in grav_init_pars and in entropy_init_pars to obtain hydrostatic
!  equilibrium. Constants g_A..D from gravz_profile.
!
!  12-feb-11/fred: older subroutine now use thermal_hs_equilibrium_ism.
!
      use EquationOfState , only: eosperturb, getmu
      use Mpicomm, only: mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(nx) :: rho,pp,lnrho
      real, dimension(1) :: fmpi1
      real :: rho0hs,muhs
      real :: g_A, g_C
      real, parameter :: g_A_cgs=4.4e-9, g_C_cgs=1.7e-9
      double precision :: g_B, g_D
      double precision, parameter :: g_B_cgs=6.172D20, g_D_cgs=3.086D21
!
!  Set up physical units.
!
      if (unit_system=='cgs') then
          g_A = g_A_cgs/unit_velocity*unit_time
          g_B = g_B_cgs/unit_length
          g_C = g_C_cgs/unit_velocity*unit_time
          g_D = g_D_cgs/unit_length
      else if (unit_system=='SI') then
        call fatal_error('initialize_gravity', &
            'SI unit conversions not inplemented')
      endif
!
!  Uses gravity profile from K. Ferriere, ApJ 497, 759, 1998, eq (34)
!  at solar radius.  (for interstellar runs)
!
      call getmu(f,muhs)
!
      if (lroot) print*, &
          'Ferriere-hs: hydrostatic equilibrium density and entropy profiles'
      T0=0.8    ! cs20/gamma_m1
      do n=n1,n2
      do m=m1,m2
        rho=rho0hs*exp(-m_u*muhs/T0/k_B*(-g_A*g_B+g_A*sqrt(g_B**2 + z(n)**2)+g_C/g_D*z(n)**2/2.))
        lnrho=log(rho)
        if (ldensity_nolog) then
          f(l1:l2,m,n,irho)=rho
        else
          f(l1:l2,m,n,ilnrho)=lnrho
        endif
        if (lentropy) then
!  Isothermal
          pp=rho*gamma_m1/gamma*T0
          call eosperturb(f,nx,pp=pp)
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
      if (lroot) print*, 'Ferriere-hs: cs2bot=',cs2bot, ' cs2top=',cs2top
!
    endsubroutine ferriere_hs
!***********************************************************************
    subroutine thermal_hs_equilibrium_ism(f)
!
!  This routine calculates a vertical profile for density for an appropriate
!  isothermal entropy designed to balance the vertical 'Ferriere' gravity
!  profile used in interstellar.
!  T0 and rho0 are chosen to ensure uv-heating approx 0.0147 at z=0.
!  Initial thermal & hydrostatice equilibrium is achieved by ensuring
!  Lambda*rho(z)=Gamma(z) in interstellar.
!
!  Requires gravz_profile='Ferriere' in gravity_simple.f90,
!  initlnrho='Galactic-hs' in density.f90 and heating_select='thermal-hs'
!  in interstellar.f90. Constants g_A..D from gravz_profile.
!
!  22-mar-10/fred: coded
!  12-aug-10/fred: moved to interstellar and added shared variables
!
      use EquationOfState , only: eoscalc, ilnrho_lnTT, getmu
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real :: lnrho,ss,lnTT,TT
      real, dimension(:), pointer :: zrho
      real, pointer :: T0hs
      real :: g_A, g_C
      real, parameter ::  g_A_cgs=4.4e-9, g_C_cgs=1.7e-9
      double precision :: g_B ,g_D
      double precision, parameter :: g_B_cgs=6.172D20 , g_D_cgs=3.086D21
      integer :: ierr
!
!  identifier
!
      if (headtt) print*,'thermal_hs_equilibrium_ism: ENTER'
!
!  Set up physical units.
!
      if (unit_system=='cgs') then
        g_A = g_A_cgs/unit_velocity*unit_time
        g_C = g_C_cgs/unit_velocity*unit_time
        g_D = g_D_cgs/unit_length
        g_B = g_B_cgs/unit_length
      else if (unit_system=='SI') then
        call fatal_error('initialize_entopy', &
            'SI unit conversions not inplemented')
      endif
!
!  Obtain vertical density profile and isothermal temperature from
!  interstellar: thermal_hs
!
      call get_shared_variable('zrho', zrho, ierr)
      if (ierr/=0) call fatal_error('thermal_hs_equilibrium_ism', &
          'there was a problem when getting zrho')
      call get_shared_variable('T0hs', T0hs, ierr)
      if (ierr/=0) call fatal_error('thermal_hs_equilibrium_ism', &
          'there was a problem when getting T0hs')
      if (lroot) print*, &
          'thermal_hs_equilibrium_ism: zrho received', &
          ' from interstellar, T0hs =',T0hs
!
!  Allocate density profile to f and derive entropy profile from
!  temperature and density
!
      do n=n1,n2
        lnrho=log(zrho(n))
        if (ldensity_nolog) then
          f(:,:,n,irho)=exp(lnrho)
        else
          f(:,:,n,ilnrho)=lnrho
        endif
!
        TT=T0hs/(g_A*g_B)* &
            (g_A*sqrt(g_B**2+(z(n))**2)+0.5*g_C*(z(n))**2/g_D)
        lnTT=log(TT)
!
        call eoscalc(ilnrho_lnTT,lnrho,lnTT,ss=ss)
!
        f(:,:,n,iss)=ss
!
      enddo
!
    endsubroutine thermal_hs_equilibrium_ism
!***********************************************************************
    subroutine galactic_hs(f,rho0hs,cs0hs,H0hs)
!
!   Density and isothermal entropy profile in hydrostatic equilibrium
!   with the galactic-hs-gravity profile set in gravity_simple.
!   Parameters cs0hs and H0hs need to be set consistently
!   both in grav_init_pars and in entropy_init_pars to obtain hydrostatic
!   equilibrium.
!
!   22-jan-10/fred
!
      use EquationOfState, only: eosperturb
      use Mpicomm, only: mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(nx) :: rho,pp,lnrho,ss
      real, dimension(1) :: fmpi1
      real :: rho0hs,cs0hs,H0hs
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
          if (ldensity_nolog) then
            ss=log(f(l1:l2,m,n,irho))
          else
            ss=f(l1:l2,m,n,ilnrho)
          endif
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
      if (lroot) print*, 'Galactic-hs: cs2bot=', cs2bot, ' cs2top=', cs2top
!
    endsubroutine galactic_hs
!***********************************************************************
    subroutine shock2d(f)
!
!  AJ: DOCUMENT ME!
!  AJ: MOVE ME TO SUBMODULE OR SPECIAL OR INITIAL CONDITION!
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
!  First quadrant:
!
        rpp(1) = 1.5d0
        rpr(1) = 1.5d0
        rpu(1) = 0.d0
        rpv(1) = 0.d0
!
!  Second quadrant:
!
        rpp(2) = 0.3d0
        rpr(2) = 0.532258064516129d0
        rpu(2) = 1.206045378311055d0
        rpv(2) = 0.0d0
!
!  Third quadrant:
!
        rpp(3) = 0.029032258064516d0
        rpr(3) = 0.137992831541219d0
        rpu(3) = 1.206045378311055d0
        rpv(3) = 1.206045378311055d0
!
!  Fourth quadrant:
!
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
    subroutine pencil_criteria_entropy()
!
!  All pencils that the Entropy module depends on are specified here.
!
!  20-11-04/anders: coded
!
      use EquationOfState, only: beta_glnrho_scaled
!
      if (lheatc_Kconst .or. lheatc_chiconst .or. lheatc_Kprof .or. &
          tau_cor>0 .or. lheatc_chitherm.or. &
          lheatc_sqrtrhochiconst) lpenc_requested(i_cp1)=.true.
      if (ldt) lpenc_requested(i_cs2)=.true.
      if (lpressuregradient_gas) lpenc_requested(i_fpres)=.true.
      if (ladvection_entropy) then
        if (lweno_transport) then
          lpenc_requested(i_rho1)=.true.
          lpenc_requested(i_transprho)=.true.
          lpenc_requested(i_transprhos)=.true.
        else
          lpenc_requested(i_ugss)=.true.
        endif
      endif
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
          cool_uniform/=0 .or. cool_newton/=0 .or. &
          (cool_ext/=0 .and. cool_int /= 0) .or. lturbulent_heat) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_TT1)=.true.
      endif
      if (pretend_lnTT) then
         lpenc_requested(i_uglnTT)=.true.
         lpenc_requested(i_divu)=.true.
         lpenc_requested(i_cv1)=.true.
         lpenc_requested(i_cp1)=.true.
      endif
      if (lgravr) then
        if (lcylindrical_coords) then
          lpenc_requested(i_rcyl_mn)=.true.
        else
          lpenc_requested(i_r_mn)=.true.
        endif
      endif
      if (lheatc_Kconst) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
      endif
      if (lheatc_chitherm) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
      endif
      if (lheatc_sqrtrhochiconst) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
      endif
      if (lheatc_Kprof) then
        if (hcond0/=0) then
          lpenc_requested(i_rho1)=.true.
          lpenc_requested(i_glnTT)=.true.
          lpenc_requested(i_del2lnTT)=.true.
          lpenc_requested(i_cp1)=.true.
          if (lmultilayer) then
            if (lgravz) then
              lpenc_requested(i_z_mn)=.true.
            elseif (lgravx) then
            else
              lpenc_requested(i_r_mn)=.true.
            endif
          endif
          if (l1davgfirst) then
             lpenc_requested(i_TT)=.true.
             lpenc_requested(i_gss)=.true.
             lpenc_requested(i_rho)=.true.
          endif
        endif
        if (chi_t/=0) then
           lpenc_requested(i_del2ss)=.true.
           lpenc_requested(i_glnrho)=.true.
           lpenc_requested(i_gss)=.true.
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
        lpenc_requested(i_cp)=.true.
      endif
      if (lheatc_hubeny) then
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_TT)=.true.
      endif
      if (lheatc_kramers) then
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_cv1)=.true.
        lpenc_requested(i_del2lnTT)=.true.
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
        if (chi_t/=0.) then
           lpenc_requested(i_gss)=.true.
           lpenc_requested(i_del2ss)=.true.
        endif
      endif
      if (lheatc_chitherm) then
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_glnrho)=.true.
        if (chi_t/=0.) then
           lpenc_requested(i_gss)=.true.
           lpenc_requested(i_del2ss)=.true.
        endif
      endif
      if (lheatc_sqrtrhochiconst) then
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_glnrho)=.true.
        if (chi_t/=0.) then
           lpenc_requested(i_gss)=.true.
           lpenc_requested(i_del2ss)=.true.
        endif
      endif
      if (lheatc_tensordiffusion) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_bij)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_cp)=.true.
      endif
      if (lheatc_shock) then
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_gss)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_gshock)=.true.
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_glnTT)=.true.
      endif
      if (lheatc_hyper3ss) lpenc_requested(i_del6ss)=.true.
      if (cooltype=='shell' .and. deltaT_poleq/=0.) then
        lpenc_requested(i_z_mn)=.true.
        lpenc_requested(i_rcyl_mn)=.true.
      endif
      if (tau_cool/=0.0 .or. cool_uniform/=0.0) then
        lpenc_requested(i_cp)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_rho)=.true.
      endif
      if (tau_cool2/=0.0) lpenc_requested(i_rho)=.true.
      if (cool_newton/=0.0) lpenc_requested(i_TT)=.true.
!
      if (maxval(abs(beta_glnrho_scaled))/=0.0) lpenc_requested(i_cs2)=.true.
!
!  magnetic chi-quenching
!
      if (chiB/=0.0) lpenc_requested(i_b2)=.true.
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
          idiag_ssmy/=0 .or. idiag_ssmx/=0 .or. idiag_ssmr/=0) &
           lpenc_diagnos(i_ss)=.true.
      if (idiag_ppmx/=0 .or. idiag_ppmy/=0 .or. idiag_ppmz/=0) &
         lpenc_diagnos(i_pp)=.true.
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
      if (idiag_fconvxy/=0 .or. idiag_fconvyxy/=0 .or. idiag_fconvzxy/=0) then
          lpenc_diagnos2d(i_rho)=.true.
          lpenc_diagnos2d(i_TT)=.true.
      endif
      if (idiag_fturbxy/=0 .or. idiag_fturbrxy/=0 .or. &
          idiag_fturbthxy/=0) then
          lpenc_diagnos2d(i_rho)=.true.
          lpenc_diagnos2d(i_TT)=.true.
      endif
      if (idiag_fradxy_kramers/=0) then
          lpenc_diagnos2d(i_rho)=.true.
          lpenc_diagnos2d(i_TT)=.true.
      endif
      if (idiag_TTm/=0 .or. idiag_TTmx/=0 .or. idiag_TTmy/=0 .or. &
          idiag_TTmz/=0 .or. idiag_TTmr/=0 .or. idiag_TTmax/=0 .or. &
          idiag_TTmin/=0 .or. idiag_uxTTmz/=0 .or.idiag_uyTTmz/=0 .or. &
          idiag_uzTTmz/=0 .or. idiag_TT2mz/=0 .or. idiag_uxTTmx/=0) &
          lpenc_diagnos(i_TT)=.true.
      if (idiag_TTmxy/=0 .or. idiag_TTmxz/=0 .or. idiag_uxTTmxy/=0 .or. &
          idiag_uyTTmxy/=0 .or. idiag_uzTTmxy/=0) &
          lpenc_diagnos2d(i_TT)=.true.
      if (idiag_yHm/=0 .or. idiag_yHmax/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_dtc/=0) lpenc_diagnos(i_cs2)=.true.
      if (idiag_TTp/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_cs2)=.true.
        lpenc_diagnos(i_rcyl_mn)=.true.
      endif
      if (idiag_cs2mphi/=0) lpenc_diagnos(i_cs2)=.true.
!
!  diagnostics for baroclinic term
!
      if (idiag_gTrms/=0.or.idiag_gTxgsrms/=0) lpenc_diagnos(i_gTT)=.true.
      if (idiag_gsrms/=0.or.idiag_gTxgsrms/=0) lpenc_diagnos(i_gss)=.true.
!
!  Cooling for cold core collapse
!
      if (lprestellar_cool_iso) then
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_divu)=.true.
        lpenc_requested(i_pp)=.true.
      endif
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
      if (lpencil_in(i_uglnTT)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_glnTT)=.true.
      endif
      if (lpencil_in(i_Ma2)) then
        lpencil_in(i_u2)=.true.
        lpencil_in(i_cs2)=.true.
      endif
      if (lpencil_in(i_fpres)) then
        if (lfpres_from_pressure) then
          lpencil_in(i_rho1)=.true.
        else
          lpencil_in(i_cs2)=.true.
          lpencil_in(i_glnrho)=.true.
          if (leos_idealgas) then
            lpencil_in(i_glnTT)=.true.
          else
            lpencil_in(i_cp1tilde)=.true.
            lpencil_in(i_gss)=.true.
          endif
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
      use EquationOfState, only: gamma_inv
      use Sub, only: u_dot_grad, grad
      use WENO_transport, only: weno_transp
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      real, dimension(nx,3) :: uu
      integer :: j
!
      intent(in) :: f
      intent(inout) :: p
! Ma2
      if (lpencil(i_Ma2)) p%Ma2=p%u2/p%cs2
! ugss
      if (lpencil(i_ugss)) then
        uu = p%uu
        if (lconst_advection) uu = uu + spread(u0_advec,1,nx)
        call u_dot_grad(f,iss,p%gss,uu,p%ugss,UPWIND=lupw_ss)
      endif
! for pretend_lnTT
      if (lpencil(i_uglnTT)) &
          call u_dot_grad(f,iss,p%glnTT,p%uu,p%uglnTT,UPWIND=lupw_ss)
! fpres
      if (lpencil(i_fpres)) then
        if (lfpres_from_pressure) then
          call grad(f,ippaux,p%fpres)
          do j=1,3
            p%fpres(:,j)=-p%rho1*p%fpres(:,j)
          enddo
        else
          if (leos_idealgas) then
            do j=1,3
              p%fpres(:,j)=-p%cs2*(p%glnrho(:,j) + p%glnTT(:,j))*gamma_inv
            enddo
!  TH: The following would work if one uncomments the intrinsic operator
!  extensions in sub.f90. Please Test.
!          p%fpres      =-p%cs2*(p%glnrho + p%glnTT)*gamma_inv
          else
            do j=1,3
              p%fpres(:,j)=-p%cs2*(p%glnrho(:,j) + p%cp1tilde*p%gss(:,j))
            enddo
          endif
        endif
      endif
!  transprhos
      if (lpencil(i_transprhos)) &
          call weno_transp(f,m,n,iss,irho,iux,iuy,iuz,p%transprhos,dx_1,dy_1,dz_1)
!
    endsubroutine calc_pencils_entropy
!***********************************************************************
    subroutine dss_dt(f,df,p)
!
!  Calculate right hand side of entropy equation,
!  ds/dt = -u.grads + [H-C + div(K*gradT) + mu0*eta*J^2 + ...]
!
!  17-sep-01/axel: coded
!   9-jun-02/axel: pressure gradient added to du/dt already here
!   2-feb-03/axel: added possibility of ionization
!
      use Diagnostics
      use EquationOfState, only: beta_glnrho_global, beta_glnrho_scaled, &
                                 gamma_inv, cs0
      use Interstellar, only: calc_heat_cool_interstellar
      use Special, only: special_calc_entropy
      use Sub
      use Viscosity, only: calc_viscous_heat
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: Hmax,gT2,gs2,gTxgs2
      real, dimension (nx,3) :: gTxgs
      real :: ztop,xi,profile_cor,uT,fradz,TTtop
      real, dimension(nx) :: ufpres, uduu
      integer :: j,ju
      integer :: i
!
      intent(inout)  :: f,p
      intent(out) :: df
!
      Hmax = 0.0
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dss_dt: SOLVE dss_dt'
      if (headtt) call identify_bcs('ss',iss)
      if (headtt) print*,'dss_dt: lnTT,cs2,cp1=', p%lnTT(1), p%cs2(1), p%cp1(1)
!
!  ``cs2/dx^2'' for timestep
!
      if (lhydro.and.lfirst.and.ldt) advec_cs2=p%cs2*dxyz_2
      if (headtt.or.ldebug) &
          print*, 'dss_dt: max(advec_cs2) =', maxval(advec_cs2)
!
!  Pressure term in momentum equation (setting lpressuregradient_gas to
!  .false. allows suppressing pressure term for test purposes).
!
      if (lhydro) then
        if (lpressuregradient_gas) then
          do j=1,3
            ju=j+iuu-1
            df(l1:l2,m,n,ju) = df(l1:l2,m,n,ju) + p%fpres(:,j)
          enddo
        endif
!
!  Add pressure force from global density gradient.
!
        if (maxval(abs(beta_glnrho_global))/=0.0) then
          if (headtt) print*, 'dss_dt: adding global pressure gradient force'
          do j=1,3
            df(l1:l2,m,n,(iux-1)+j) = df(l1:l2,m,n,(iux-1)+j) &
                - p%cs2*beta_glnrho_scaled(j)
          enddo
        endif
!
!  Velocity damping in the coronal heating zone.
!
        if (tau_cor>0) then
          ztop=xyz0(3)+Lxyz(3)
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
!  Advection of entropy.
!  If pretend_lnTT=.true., we pretend that ss is actually lnTT
!  Otherwise, in the regular case with entropy, s is the dimensional
!  specific entropy, i.e. it is not divided by cp.
!  NOTE: in the entropy module is it lnTT that is advanced, so
!  there are additional cv1 terms on the right hand side.
!
      if (ladvection_entropy) then
        if (pretend_lnTT) then
          df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) - p%divu*gamma_m1-p%uglnTT
          if (lconst_advection) &
            call fatal_error('dss_dt', 'uniform advection is not implemented with lnTT')
        else
          if (lweno_transport) then
            df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
                - p%transprhos*p%rho1 + p%ss*p%rho1*p%transprho
            if (lconst_advection) &
              call fatal_error('dss_dt', 'uniform advection is not implemented with WENO transport')
          else
            df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) - p%ugss
          endif
        endif
      endif
!
!  Calculate viscous contribution to entropy.
!
      if (lviscosity .and. lviscosity_heat) call calc_viscous_heat(f,df,p,Hmax)
!
!  Entry possibility for "personal" entries.
!  In that case you'd need to provide your own "special" routine.
!
      if (lspecial) call special_calc_entropy(f,df,p)
!
!  Thermal conduction delegated to different subroutines.
!
      if (lheatc_Kprof)    call calc_heatcond(f,df,p)
      if (lheatc_Kconst)   call calc_heatcond_constK(df,p)
      if (lheatc_chiconst) call calc_heatcond_constchi(df,p)
      if (lheatc_chitherm) call calc_heatcond_thermchi(df,p)
      if (lheatc_sqrtrhochiconst) call calc_heatcond_sqrtrhochi(df,p)
      if (lheatc_shock)    call calc_heatcond_shock(df,p)
      if (lheatc_hyper3ss) call calc_heatcond_hyper3(df,p)
      if (lheatc_spitzer)  call calc_heatcond_spitzer(df,p)
      if (lheatc_hubeny)   call calc_heatcond_hubeny(df,p)
      if (lheatc_kramers)  call calc_heatcond_kramers(df,p)
      if (lheatc_corona) then
        call calc_heatcond_spitzer(df,p)
        call newton_cool(df,p)
        call calc_heat_cool_RTV(df,p)
      endif
      if (lheatc_tensordiffusion) call calc_heatcond_tensor(df,p)
      if (lheatc_hyper3ss_polar) call calc_heatcond_hyper3_polar(f,df)
      if (lheatc_hyper3ss_mesh)  call calc_heatcond_hyper3_mesh(f,df)
      if (lheatc_hyper3ss_aniso) call calc_heatcond_hyper3_aniso(f,df)
!
!  Explicit heating/cooling terms.
!
      if ((luminosity/=0.0) .or. (cool/=0.0) .or. &
          (tau_cor/=0.0) .or. (tauheat_buffer/=0.0) .or. &
          heat_uniform/=0.0 .or. tau_cool/=0.0 .or. &
          cool_uniform/=0.0 .or. cool_newton/=0.0 .or. &
          (cool_ext/=0.0 .and. cool_int/=0.0) .or. lturbulent_heat .or. &
          (tau_cool2 /=0)) &
          call calc_heat_cool(df,p,Hmax)
      if (tdown/=0.0) call newton_cool(df,p)
      if (cool_RTV/=0.0) call calc_heat_cool_RTV(df,p)
      if (lprestellar_cool_iso) call calc_heat_cool_prestellar(f,df,p)
!
!  Interstellar radiative cooling and UV heating.
!
      if (linterstellar) call calc_heat_cool_interstellar(f,df,p,Hmax)
!
!  Possibility of entropy relaxation in exterior region.
!
      if (tau_ss_exterior/=0.0) call calc_tau_ss_exterior(df,p)
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
!      if (lfirst.and.ldt) dt1_max=max(dt1_max,Hmax/p%ee/cdts)
!
!  Calculate entropy related diagnostics.
!
      if (ldiagnos) then
        !uT=unit_temperature !(define shorthand to avoid long lines below)
        uT=1. !(AB: for the time being; to keep compatible with auto-test
        if (idiag_ssmax/=0)  call max_mn_name(p%ss*uT,idiag_ssmax)
        if (idiag_ssmin/=0)  call max_mn_name(-p%ss*uT,idiag_ssmin,lneg=.true.)
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
        if (idiag_ufpresm/=0) then
            ufpres=0.
            do i = 1, 3
              ufpres=ufpres+p%uu(:,i)*p%fpres(:,i)
            enddo
            call sum_mn_name(p%rho*ufpres,idiag_ufpresm)
        endif
        if (idiag_uduum/=0) then
            uduu=0.
            do i = 1, 3
              uduu=uduu+p%uu(:,i)*df(l1:l2,m,n,iux-1+i)
            enddo
            call sum_mn_name(p%rho*uduu,idiag_uduum)
        endif
!
!  Analysis for the baroclinic term.
!
        if (idiag_gTrms/=0) then
          call dot2(p%gTT,gT2)
          call sum_mn_name(gT2,idiag_gTrms,lsqrt=.true.)
        endif
!
        if (idiag_gsrms/=0) then
          call dot2(p%gss,gs2)
          call sum_mn_name(gs2,idiag_gsrms,lsqrt=.true.)
        endif
!
        if (idiag_gTxgsrms/=0) then
          call cross(p%gTT,p%gss,gTxgs)
          call dot2(gTxgs,gTxgs2)
          call sum_mn_name(gTxgs2,idiag_gTxgsrms,lsqrt=.true.)
        endif
!
!  Radiative heat flux at the bottom (assume here that hcond=hcond0=const).
!
        if (idiag_fradbot/=0) then
          if (lfirst_proc_z.and.n==n1) then
            fradz=sum(-hcond0*p%TT*p%glnTT(:,3)*dsurfxy)
          else
            fradz=0.
          endif
          call surf_mn_name(fradz,idiag_fradbot)
        endif
!
!  Radiative heat flux at the top (assume here that hcond=hcond0=const).
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
!  Mean temperature at the top.
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
!  Calculate integrated temperature in in limited radial range.
!
        if (idiag_TTp/=0) call sum_lim_mn_name(p%rho*p%cs2*gamma_inv,idiag_TTp,p)
      endif
!
!  1-D averages.
!
      if (l1davgfirst) then
        call xysum_mn_name_z(-hcond0*p%TT*p%glnTT(:,3),idiag_fradz)
        call xysum_mn_name_z(p%rho*p%uu(:,3)*p%TT,idiag_fconvz)
        call yzsum_mn_name_x(p%ss,idiag_ssmx)
        call xzsum_mn_name_y(p%ss,idiag_ssmy)
        call xysum_mn_name_z(p%ss,idiag_ssmz)
        call yzsum_mn_name_x(p%pp,idiag_ppmx)
        call xzsum_mn_name_y(p%pp,idiag_ppmy)
        call xysum_mn_name_z(p%pp,idiag_ppmz)
        call yzsum_mn_name_x(p%TT,idiag_TTmx)
        call xzsum_mn_name_y(p%TT,idiag_TTmy)
        call xysum_mn_name_z(p%TT,idiag_TTmz)
        call xysum_mn_name_z(p%TT**2,idiag_TT2mz)
        call xysum_mn_name_z(p%uu(:,1)*p%TT,idiag_uxTTmz)
        call yzsum_mn_name_x(p%uu(:,1)*p%TT,idiag_uxTTmx)
        call xysum_mn_name_z(p%uu(:,2)*p%TT,idiag_uyTTmz)
        call xysum_mn_name_z(p%uu(:,3)*p%TT,idiag_uzTTmz)
        if (idiag_ssmr/=0)  call phizsum_mn_name_r(p%ss,idiag_ssmr)
        if (idiag_TTmr/=0)  call phizsum_mn_name_r(p%TT,idiag_TTmr)
      endif
!
!  2-D averages.
!
      if (l2davgfirst) then
        if (idiag_TTmxy/=0) call zsum_mn_name_xy(p%TT,idiag_TTmxy)
        if (idiag_TTmxz/=0) call ysum_mn_name_xz(p%TT,idiag_TTmxz)
        if (idiag_ssmxy/=0) call zsum_mn_name_xy(p%ss,idiag_ssmxy)
        if (idiag_ssmxz/=0) call ysum_mn_name_xz(p%ss,idiag_ssmxz)
        if (idiag_uxTTmxy/=0) &
            call zsum_mn_name_xy(p%uu(:,1)*p%TT,idiag_uxTTmxy)
        if (idiag_uyTTmxy/=0) &
            call zsum_mn_name_xy(p%uu(:,2)*p%TT,idiag_uyTTmxy)
        if (idiag_uzTTmxy/=0) &
            call zsum_mn_name_xy(p%uu(:,3)*p%TT,idiag_uzTTmxy)
        if (idiag_fconvxy/=0) &
            call zsum_mn_name_xy(p%rho*p%uu(:,1)*p%TT,idiag_fconvxy)
        if (idiag_fconvyxy/=0) &
            call zsum_mn_name_xy(p%rho*p%uu(:,2)*p%TT,idiag_fconvyxy)
        if (idiag_fconvzxy/=0) &
            call zsum_mn_name_xy(p%rho*p%uu(:,3)*p%TT,idiag_fconvzxy)
      endif
!
    endsubroutine dss_dt
!***********************************************************************
    subroutine calc_lentropy_pars(f)
!
!  Calculate <s>, which is needed for diffusion with respect to xy-flucts.
!
!  17-apr-10/axel: adapted from calc_lmagnetic_pars
!
      use Deriv, only: der_z, der2_z
      use Mpicomm, only: mpiallreduce_sum
      use EquationOfState, only : lnrho0, cs20, get_cv1
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: nxy=nxgrid*nygrid
      integer :: nyz=nygrid*nzgrid
      integer :: l,m,n
      real :: fact, cv1
      real, dimension (mz) :: cs2mz_tmp, ssmz1_tmp
      real, dimension (mx) :: cs2mx_tmp, ssmx1_tmp
      real, dimension (mx,my) :: cs2mxy_tmp, ssmxy1_tmp
!
      intent(in) :: f
!
!  Compute horizontal average of entropy. Include the ghost zones,
!  because they have just been set.
!
      if (lcalc_ssmean) then
        fact=1./nxy
        do n=1,mz
          ssmz(n)=fact*sum(f(l1:l2,m1:m2,n,iss))
        enddo
!
!  Communication over all processors in the xy plane.
!
        if (nprocx>1.or.nprocy>1) then
          call mpiallreduce_sum(ssmz,ssmz1_tmp,mz,idir=12)
          ssmz=ssmz1_tmp
        endif
!
!  Compute first and second derivatives.
!
        gssmz(:,1:2)=0.
        call der_z(ssmz,gssmz(:,3))
        call der2_z(ssmz,del2ssmz)
      endif
!
!  Compute yz- and z-averages of entropy.
!
      if (lcalc_ssmeanxy) then
!
!  Radial (yz) average
!
        fact=1./nyz
        do l=1,mx
          ssmx(l)=fact*sum(f(l,m1:m2,n1:n2,iss))
        enddo
!
!  Communication over all processors in the yz plane.
!
        if (nprocy>1.or.nprocz>1) then
          call mpiallreduce_sum(ssmx,ssmx1_tmp,mx,idir=23)
          ssmx=ssmx1_tmp
        endif
!
!  Azimuthal (z) average
!
        fact=1./nzgrid
        do l=1,mx
          do m=1,my
            ssmxy(l,m)=fact*sum(f(l,m,n1:n2,iss))
          enddo
        enddo
!
        if (nprocz>1) then
!
          call mpiallreduce_sum(ssmxy,ssmxy1_tmp,(/mx,my/),idir=3)
          ssmxy = ssmxy1_tmp
!
        endif
      endif
!
!  Compute average sound speed cs2(x) and cs2(x,y)
!
      if (lcalc_cs2mean) then
        call get_cv1(cv1)
        fact=1./nyz
        do l=1,mx
          cs2mx(l)=fact*sum(cs20*exp(gamma_m1*(f(l,m1:m2,n1:n2,ilnrho) &
              -lnrho0)+cv1*f(l,m1:m2,n1:n2,iss)))
        enddo
!
!  Communication over all processors in the yz plane.
!
        if (nprocy>1.or.nprocz>1) then
          call mpiallreduce_sum(cs2mx,cs2mx_tmp,mx,idir=23)
          cs2mx=cs2mx_tmp
        endif
!
!  Do 2D averages at the same time
!
        fact=1./nzgrid
        cs2mxy = 0.
!
        do l=1,mx
          do m=1,my
            cs2mxy(l,m)=fact*sum(cs20*exp(gamma_m1*(f(l,m,n1:n2,ilnrho) &
                -lnrho0)+cv1*f(l,m,n1:n2,iss)))
          enddo
        enddo
!
        if (nprocz>1) then
!
          call mpiallreduce_sum(cs2mxy,cs2mxy_tmp,(/mx,my/),idir=3)
          cs2mxy = cs2mxy_tmp
!
        endif
!
      endif
!
!  Compute average sound speed cs2(z)
!
      if (lcalc_cs2mz_mean) then
        fact=1./nxy
        cs2mz=0.
        if (ldensity_nolog) then
          do n=1,mz
            cs2mz(n)=fact*sum(cs20*exp(gamma_m1*(alog(f(l1:l2,m1:m2,n,irho)) &
                -lnrho0)+cv1*f(l1:l2,m1:m2,n,iss)))
          enddo
        else
          do n=1,mz
            cs2mz(n)=fact*sum(cs20*exp(gamma_m1*(f(l1:l2,m1:m2,n,ilnrho) &
                -lnrho0)+cv1*f(l1:l2,m1:m2,n,iss)))
          enddo
        endif
!
!  communicate over x and y directions
!
        if (nprocx>1.or.nprocy>1) then
          call mpiallreduce_sum(cs2mz,cs2mz_tmp,mz,idir=12)
          cs2mz=cs2mz_tmp
        endif
      endif
!
!  Compute volume average of entropy.
!
      if (lcalc_ss_volaverage) then
        ss_volaverage_tmp=sum(f(l1:l2,m1:m2,n1:n2,iss))/nwgrid
        call mpiallreduce_sum(ss_volaverage_tmp,ss_volaverage)
      endif
!
    endsubroutine calc_lentropy_pars
!***********************************************************************
    subroutine set_border_entropy(f,df,p)
!
!  Calculates the driving term for the border profile
!  of the ss variable.
!
!  28-jul-06/wlad: coded
!
      use BorderProfiles, only: border_driving, set_border_initcond
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(nx) :: f_target
!
      select case (borderss)
!
      case ('zero','0')
        f_target=0.0
      case ('constant')
        f_target=ss_const
      case ('initial-condition')
        call set_border_initcond(f,iss,f_target)
      case ('nothing')
        if (lroot.and.ip<=5) &
            print*, "set_border_entropy: borderss='nothing'"
      case default
        write(unit=errormsg,fmt=*) &
            'set_border_entropy: No such value for borderss: ', trim(borderss)
        call fatal_error('set_border_entropy',errormsg)
      endselect
!
      if (borderss/='nothing') then
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
      use Diagnostics, only: max_mn_name
      use Sub, only: dot
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: thdiff,g2
!
      intent(out) :: df
!
!  Check that chi is ok.
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
!  Add heat conduction to entropy equation.
!
      df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)+thdiff
      if (headtt) print*,'calc_heatcond_constchi: added thdiff'
!
!  Check maximum diffusion from thermal diffusion.
!  With heat conduction, the second-order term for entropy is
!  gamma*chi*del2ss.
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
    subroutine calc_heatcond_thermchi(df,p)
!
!  Adapted from Heat conduction for constant value to
!  include temperature dependence to handle high temperatures
!  in hot diffuse cores of SN remnants in interstellar chi propto sqrt(T)
!  This routine also adds in turbulent diffusion, if chi_t /= 0.
!  Ds/Dt = ... + 1/(rho*T) grad(flux), where
!  flux = chi_th*rho*gradT + chi_t*rho*T*grads
!  This routine is currently not correct when ionization is used.
!
!  19-mar-10/fred: adapted from calc_heatcond_constchi - still need to test
!                  physics
!  12-mar-06/axel: used p%glnTT and p%del2lnTT, so that general cp work ok
!
      use Diagnostics, only: max_mn_name
      use Sub, only: dot
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: thdiff,g2,thchi
!
      intent(out) :: df
!
!  Check that chi is ok.
!
      if (headtt) print*,'calc_heatcond_thermchi: chi_th=',chi_th
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
!  Note: need thermally sensitive diffusion without magnetic field
!  for interstellar hydro runs to contrain SNr core temp
!
!
      thchi=chi_th*sqrt(exp(p%lnTT))
!      thchi=max(chi_th*(exp(p%lnTT))**0.5,dxmax*0.5)
      if (pretend_lnTT) then
        call dot(p%glnrho+p%glnTT,p%glnTT,g2)
        thdiff=gamma*thchi*(p%del2lnTT+g2)
        if (chi_t/=0.) then
          call dot(p%glnrho+p%glnTT,p%gss,g2)
          thdiff=thdiff+chi_t*(p%del2ss+g2)
        endif
      else
        call dot(p%glnrho+p%glnTT,p%glnTT,g2)
!AB:  divide by p%cp1, since we don't have cp here.
        thdiff=thchi*(p%del2lnTT+g2)/p%cp1
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
!  Add heat conduction to entropy equation.
!
      df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)+thdiff
      if (headtt) print*,'calc_heatcond_thermchi: added thdiff'
!
!  Check maximum diffusion from thermal diffusion.
!  With heat conduction, the second-order term for entropy is
!  gamma*chi*del2ss.
!
      if (lfirst.and.ldt) then
        if (leos_idealgas) then
          diffus_chi=diffus_chi+(gamma*thchi+chi_t)*dxyz_2
        else
          diffus_chi=diffus_chi+(thchi+chi_t)*dxyz_2
        endif
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_thermchi
!***********************************************************************
    subroutine calc_heatcond_sqrtrhochi(df,p)
!
!  Adapted from Heat conduction for constant value to
!  include temperature dependence to handle high temperatures
!  in hot diffuse cores of SN remnants in interstellar chi propto sqrt(rho)
!  This routine also adds in turbulent diffusion, if chi_t /= 0.
!  Ds/Dt = ... + 1/(rho*T) grad(flux), where
!  flux = chi_th*rho*gradT + chi_t*rho*T*grads
!  This routine is currently not correct when ionization is used.
!
!  19-mar-10/fred: adapted from calc_heatcond_constchi - still need to test
!                  physics
!  12-mar-06/axel: used p%glnTT and p%del2lnTT, so that general cp work ok
!
      use Diagnostics, only: max_mn_name
      use Sub, only: dot
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: thdiff,g2,rhochi
!
      intent(out) :: df
!
!  Check that chi is ok.
!
      if (headtt) print*,'calc_heatcond_sqrtrhochi: chi_rho=',chi_rho
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
!  Note: need thermally sensitive diffusion without magnetic field
!  for interstellar hydro runs to contrain SNr core temp
!
      rhochi=chi_rho*sqrt(p%rho1)
      if (pretend_lnTT) then
        call dot(p%glnrho+p%glnTT,p%glnTT,g2)
        thdiff=gamma*rhochi*(p%del2lnTT+g2)
        if (chi_t/=0.) then
          call dot(p%glnrho+p%glnTT,p%gss,g2)
          thdiff=thdiff+chi_t*(p%del2ss+g2)
        endif
      else
        call dot(p%glnrho+p%glnTT,p%glnTT,g2)
!AB:  divide by p%cp1, since we don't have cp here.
        thdiff=rhochi*(p%del2lnTT+g2)/p%cp1
        if (chi_t/=0.) then
          call dot(p%glnrho+p%glnTT,p%gss,g2)
!
!  Provisional expression for magnetic chi_t quenching;
!  derivatives of B are still missing.
!
          if (chiB==0.) then
            thdiff=thdiff+chi_t*(p%del2ss+g2)
          else
            thdiff=thdiff+chi_t*(p%del2ss+g2)/(1.+chiB*p%b2)
          endif
        endif
      endif
!
!  Add heat conduction to entropy equation.
!
      df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)+thdiff
      if (headtt) print*,'calc_heatcond_sqrtrhochi: added thdiff'
!
!  Check maximum diffusion from thermal diffusion.
!  With heat conduction, the second-order term for entropy is
!  gamma*chi*del2ss.
!
      if (lfirst.and.ldt) then
        if (leos_idealgas) then
          diffus_chi=diffus_chi+(gamma*rhochi+chi_t)*dxyz_2
        else
          diffus_chi=diffus_chi+(rhochi+chi_t)*dxyz_2
        endif
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_sqrtrhochi
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
!  Check that chi_hyper3 is ok.
!
      if (headtt) print*, 'calc_heatcond_hyper3: chi_hyper3=', chi_hyper3
!
!  Heat conduction.
!
      thdiff = chi_hyper3 * p%del6ss
!
!  Add heat conduction to entropy equation.
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
      if (headtt) print*,'calc_heatcond_hyper3: added thdiff'
!
!  Check maximum diffusion from thermal diffusion.
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
!  Check that chi_hyper3_aniso is ok.
!
      if (headtt) print*, &
           'calc_heatcond_hyper3_aniso: chi_hyper3_aniso=', &
           chi_hyper3_aniso
!
!  Heat conduction.
!
      call del6fj(f,chi_hyper3_aniso,iss,tmp)
      thdiff = tmp
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
!
      if (lfirst.and.ldt) diffus_chi3=diffus_chi3+ &
           (chi_hyper3_aniso(1)*dx_1(l1:l2)**6 + &
            chi_hyper3_aniso(2)*dy_1(  m  )**6 + &
            chi_hyper3_aniso(3)*dz_1(  n  )**6)
!
    endsubroutine calc_heatcond_hyper3_aniso
!***********************************************************************
    subroutine calc_heatcond_hyper3_polar(f,df)
!
!  Naive hyperdiffusivity of entropy in polar coordinates.
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
      thdiff=0.
      do j=1,3
        call der6(f,iss,tmp,j,IGNOREDX=.true.)
        thdiff = thdiff + chi_hyper3*pi4_1/60.*tmp*dline_1(:,j)**2
      enddo
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
!
      if (headtt) print*,'calc_heatcond_hyper3: added thdiff'
!
      if (lfirst.and.ldt) &
          diffus_chi3=diffus_chi3+chi_hyper3*pi4_1*dxyz_2
!
    endsubroutine calc_heatcond_hyper3_polar
!***********************************************************************
    subroutine calc_heatcond_hyper3_mesh(f,df)
!
!  Naive resolution-independent hyperdiffusivity of entropy
!
!  25-dec-10/wlad: coded
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
      thdiff=0.0
      do j=1,3
        call der6(f,iss,tmp,j,IGNOREDX=.true.)
        thdiff = thdiff + chi_hyper3_mesh*pi5_1/60.*tmp*dline_1(:,j)
      enddo
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
!
      if (headtt) print*,'calc_heatcond_hyper3: added thdiff'
!
      if (lfirst.and.ldt) &
          advec_hypermesh_ss=chi_hyper3_mesh*pi5_1*sqrt(dxyz_2)
!
    endsubroutine calc_heatcond_hyper3_mesh
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
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: thdiff,g2,gshockglnTT
!
      intent(in) :: p
      intent(out) :: df
!
!  Check that chi is ok.
!
      if (headtt) print*,'calc_heatcond_shock: chi_t,chi_shock=',chi_t,chi_shock
!
!  Calculate terms for shock diffusion:
!  Ds/Dt = ... + chi_shock*[del2ss + (glnchi_shock+glnpp).gss]
!
      call dot(p%gshock,p%glnTT,gshockglnTT)
      call dot(p%glnrho+p%glnTT,p%glnTT,g2)
!
!  Shock entropy diffusivity.
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
!  Add heat conduction to entropy equation.
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
      if (headtt) print*,'calc_heatcond_shock: added thdiff'
!
!  Check maximum diffusion from thermal diffusion.
!  With heat conduction, the second-order term for entropy is
!  gamma*chi*del2ss.
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
!  Heat conduction.
!
!   8-jul-02/axel: adapted from Wolfgang's more complex version
!  30-mar-06/ngrs: simplified calculations using p%glnTT and p%del2lnTT
!
      use Diagnostics, only: max_mn_name
      use Sub, only: dot
!
      type (pencil_case) :: p
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: chix
      real, dimension (nx) :: thdiff,g2
      real, dimension (nx) :: hcond
!
      intent(in) :: p
      intent(out) :: df
!
!  This particular version assumes a simple polytrope, so mpoly is known.
!
      if (tau_diff==0) then
        hcond=Kbot
      else
        hcond=Kbot/tau_diff
      endif
!
      if (headtt) then
        print*,'calc_heatcond_constK: hcond=', maxval(hcond)
      endif
!
!  Heat conduction
!  Note: these routines require revision when ionization turned on
!
! NB: the following left in for the record, but the version below,
!     using del2lnTT & glnTT, is simpler
!
!chix = p%rho1*hcond*p%cp1                    ! chix = K/(cp rho)
!glnT = gamma*p%gss*spread(p%cp1,2,3) + gamma_m1*p%glnrho ! grad ln(T)
!glnThcond = glnT !... + glhc/spread(hcond,2,3)    ! grad ln(T*hcond)
!call dot(glnT,glnThcond,g2)
!thdiff =  p%rho1*hcond * (gamma*p%del2ss*p%cp1 + gamma_m1*p%del2lnrho + g2)
!
!  diffusion of the form:
!  rho*T*Ds/Dt = ... + nab.(K*gradT)
!        Ds/Dt = ... + K/rho*[del2lnTT+(glnTT)^2]
!
! NB: chix = K/(cp rho) is needed for diffus_chi calculation
!
!  Put empirical heat transport suppression by the B-field.
!
      if (chiB==0.) then
        chix = p%rho1*hcond*p%cp1
      else
        chix = p%rho1*hcond*p%cp1/(1.+chiB*p%b2)
      endif
      call dot(p%glnTT,p%glnTT,g2)
!
      if (pretend_lnTT) then
         thdiff = gamma*chix * (p%del2lnTT + g2)
      else
         thdiff = p%rho1*hcond * (p%del2lnTT + g2)
      endif
!
!  Add heat conduction to entropy equation.
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
      if (headtt) print*,'calc_heatcond_constK: added thdiff'
!
!  Check maximum diffusion from thermal diffusion.
!  With heat conduction, the second-order term for entropy is
!  gamma*chix*del2ss.
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
!  to magnetic field lines.
!
!  See: Solar MHD; Priest 1982
!
!  10-feb-04/bing: coded
!
      use EquationOfState, only: gamma, gamma_m1
      use IO, only: output_pencil
      use Sub, only: dot2_mn, multsv_mn, tensor_diffusion_coef
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: gvKpara,gvKperp,tmpv,tmpv2
      real, dimension (nx) :: bb2,thdiff,b1
      real, dimension (nx) :: tmps,vKpara,vKperp
!
      integer ::i,j
      type (pencil_case) :: p
!
      intent(in) :: p
      intent(out) :: df
!
!  Calculate variable diffusion coefficients along pencils.
!
      call dot2_mn(p%bb,bb2)
      b1=1./max(tiny(bb2),bb2)
!
      vKpara = Kgpara * exp(p%lnTT)**3.5
      vKperp = Kgperp * b1*exp(2*p%lnrho+0.5*p%lnTT)
!
!  Calculate gradient of variable diffusion coefficients.
!
      tmps = 3.5 * vKpara
      call multsv_mn(tmps,p%glnTT,gvKpara)
!
      do i=1,3
         tmpv(:,i)=0.
         do j=1,3
            tmpv(:,i)=tmpv(:,i) + p%bb(:,j)*p%bij(:,j,i)
         enddo
      enddo
      call multsv_mn(2*b1,tmpv,tmpv2)
      tmpv=2.*p%glnrho+0.5*p%glnTT-tmpv2
      call multsv_mn(vKperp,tmpv,gvKperp)
!
!  Calculate diffusion term.
!
      call tensor_diffusion_coef(p%glnTT,p%hlnTT,p%bij,p%bb,vKperp,vKpara,thdiff,GVKPERP=gvKperp,GVKPARA=gvKpara)
!
      if (lfirst .and. ip == 13) then
        call output_pencil(trim(directory)//'/spitzer.dat',thdiff,1)
        call output_pencil(trim(directory)//'/viscous.dat',p%visc_heat*exp(p%lnrho),1)
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
!  to magnetic field lines.
!
!  24-aug-09/bing: moved from dss_dt to here
!
      use Diagnostics, only: max_mn_name
      use Sub, only: tensor_diffusion_coef, dot, dot2
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
      where (p%rho <= tiny(0.0)) p%rho1=0.0
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
      use EquationOfState, only: gamma, gamma_m1
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
      if (pretend_lnTT) call fatal_error('calc_heatcond_hubeny', &
          'not implemented when pretend_lnTT = T')
!
      kappa0_cgs=2.0e-6  !cm2/g
      kappa0=kappa0_cgs*unit_density/unit_length
      kappa=kappa0*p%TT**2
!
!  Optical Depth tau=kappa*rho*H
!  If we are using 2D, the pencil value p%rho is actually
!  sigma, the column density, sigma=rho*2*H.
!
      if (nzgrid==1) then
        tau = 0.5*kappa*p%rho
!
!  Analytical gray description of Hubeny (1990).
!  a1 is the optically thick contribution,
!  a3 the optically thin one.
!
        a1=0.375*tau ; a2=0.433013 ; a3=0.25/tau
!
        cooling = 2*sigmaSB*p%TT**4/(a1+a2+a3)
!
        df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) - cool_fac*cooling
!
      else
        call fatal_error('calc_heat_hubeny', &
            'opacity not yet implemented for 3D')
      endif
!
    endsubroutine calc_heatcond_hubeny
!***********************************************************************
    subroutine calc_heatcond_kramers(df,p)
!
!  Heat conduction using Kramers' opacity law
!
!  23-feb-11/pete: coded
!
      use Diagnostics
      use IO, only: output_pencil
      use Sub, only: dot, notanumber
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: thdiff, chix, g2
!
      intent(in) :: p
      intent(out) :: df
!
!  Diffusion of the form
!      rho*T*Ds/Dt = ... + nab.(K*gradT),
!  where
!      K = K_0*(T**6.5/rho**2)**n.
!  In reality n=1, but we may need to use n\=1 for numerical reasons.
!
!  Here chix = K/(cp rho) is needed for diffus_chi calculation.
!
      chix = p%cp1*hcond0_kramers*p%rho1**(2.*nkramers+1.)*p%TT**(6.5*nkramers)
      call dot(-2.*nkramers*p%glnrho+6.5*nkramers*p%glnTT,p%glnTT,g2)
      if (pretend_lnTT) then
        thdiff = p%cv1*hcond0_kramers*p%rho1**(2.*nkramers+1.)*p%TT**(6.5*nkramers)*(p%del2lnTT+g2)
      else
        thdiff = hcond0_kramers*p%rho1**(2.*nkramers+1.)*p%TT**(6.5*nkramers)*(p%del2lnTT+g2)
      endif
!
!  Write radiative flux array.
!
      if (l1davgfirst) then
        call xysum_mn_name_z(-chix*p%rho*p%TT*p%glnTT(:,3)/p%cp1,idiag_fradz_kramers)
      endif
!
!  2d-averages
!
      if (l2davgfirst) then
        if (idiag_fradxy_kramers/=0) call zsum_mn_name_xy(-chix*p%rho*p%TT*p%glnTT(:,1)/p%cp1,idiag_fradxy_kramers)
      endif
!
!  Check for NaNs initially.
!
      if (headt .and. (hcond0_kramers/=0.0)) then
        if (notanumber(p%rho1))    print*,'calc_heatcond_kramers: NaNs in rho1'
        if (notanumber(chix))      print*,'calc_heatcond_kramers: NaNs in chix'
        if (notanumber(p%del2ss))  print*,'calc_heatcond_kramers: NaNs in del2ss'
        if (notanumber(p%TT))      print*,'calc_heatcond_kramers: NaNs in TT'
        if (notanumber(p%glnTT))   print*,'calc_heatcond_kramers: NaNs in glnT'
        if (notanumber(g2))        print*,'calc_heatcond_kramers: NaNs in g2'
        if (notanumber(thdiff))    print*,'calc_heatcond_kramers: NaNs in thdiff'
!
!  Most of these should trigger the following trap.
!
        if (notanumber(thdiff)) then
          print*, 'calc_heatcond_kramers: m,n,y(m),z(n)=', m, n, y(m), z(n)
          call fatal_error('calc_heatcond_kramers','NaNs in thdiff')
        endif
      endif
!
!  At the end of this routine, add all contribution to
!  thermal diffusion on the rhs of the entropy equation,
!  so Ds/Dt = ... + thdiff = ... + (...)/(rho*T)
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
!
      if (headtt) print*,'calc_heatcond_kramers: added thdiff'
!
!  Check maximum diffusion from thermal diffusion.
!  NB: With heat conduction, the second-order term for entropy is
!    gamma*chix*del2ss.
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+(gamma*chix+chi_t)*dxyz_2
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_kramers
!***********************************************************************
    subroutine calc_heatcond(f,df,p)
!
!  In this routine general heat conduction profiles are being provided.
!
!  17-sep-01/axel: coded
!  14-jul-05/axel: corrected expression for chi_t diffusion.
!  30-mar-06/ngrs: simplified calculations using p%glnTT and p%del2lnTT
!
      use Diagnostics
      use IO, only: output_pencil
      use Sub, only: dot, notanumber, g2ij, write_zprof, write_zprof_once
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: glnThcond,glhc,glchit_prof,gss1,glchit_aniso_prof
      real, dimension (nx) :: chix
      real, dimension (nx) :: thdiff,g2,del2ss1
      real, dimension (nx) :: hcond,chit_prof,chit_aniso_prof
      real, dimension (nx,3,3) :: tmp
      !real, save :: z_prev=-1.23e20
      real :: s2,c2,sc
      integer :: j,ix
!
      save :: hcond, glhc, chit_prof, glchit_prof, chit_aniso_prof, glchit_aniso_prof
!
      intent(in) :: p
      intent(out) :: df
!
!  Heat conduction / entropy diffusion
!
      if (hcond0 == 0) then
        chix = 0.0
        thdiff = 0.0
        hcond = 0.0
        glhc = 0.0
      else
        if (headtt) then
          print*,'calc_heatcond: hcond0=',hcond0
          print*,'calc_heatcond: lgravz=',lgravz
          if (lgravz) print*,'calc_heatcond: Fbot,Ftop=',Fbot,Ftop
        endif
!
!  Assume that a vertical K profile is given if lgravz is true.
!
        if (lgravz) then
!
! DM+GG Added routines to compute hcond and gradloghcond_zprof.
! When called for the first time calculate z dependent profile of
! heat conductivity. For all other times use this stored arrays.
! We also write the z dependent profile of heatconduction and
! gradient-logarithm of heat conduction.
!
          if (lfirstcall_hcond.and.lfirstpoint) then
            if (.not.lhcond_global) then
              call get_gravz_heatcond()
              call write_zprof_once('hcond',hcond_zprof)
              call write_zprof_once('gloghcond',gradloghcond_zprof(:,3))
            endif
            if (chi_t/=0.0) call get_gravz_chit()
            lfirstcall_hcond=.false.
          endif
!
          if (lhcond_global) then
            hcond=f(l1:l2,m,n,iglobal_hcond)
            glhc=f(l1:l2,m,n,iglobal_glhc:iglobal_glhc+2)
          else
            hcond=hcond_zprof(n)
              do ix=1,nx
                glhc(ix,:)=gradloghcond_zprof(n,:)
              enddo
          endif
          if (chi_t/= 0.0) then
            chit_prof=chit_zprof(n)
              do ix=1,nx
                glchit_prof(ix,:)=gradlogchit_zprof(n,:)
              enddo
          endif
! If not gravz, using or not hcond_global
        elseif (lgravx) then
          if (lfirstcall_hcond.and.lfirstpoint) then
            if (.not.lhcond_global) then
               call get_gravx_heatcond()
            endif
            lfirstcall_hcond=.false.
          endif
          if (lhcond_global) then
            hcond=f(l1:l2,m,n,iglobal_hcond)
            glhc=f(l1:l2,m,n,iglobal_glhc:iglobal_glhc+2)
          else
            hcond = hcond_xprof(l1:l2)
            glhc = gradloghcond_xprof(l1:l2,:)
          endif
          if (chi_t/=0.0) then
            call chit_profile(chit_prof)
            call gradlogchit_profile(glchit_prof)
          endif
          if (chit_aniso/=0.0) then
            call chit_aniso_profile(chit_aniso_prof)
            call gradlogchit_aniso_profile(glchit_aniso_prof)
          endif
        else
          if (lhcond_global) then
            hcond=f(l1:l2,m,n,iglobal_hcond)
            glhc=f(l1:l2,m,n,iglobal_glhc:iglobal_glhc+2)
          else
            call heatcond(hcond,p)
            call gradloghcond(glhc,p)
          endif
          if (chi_t/=0.0) then
            call chit_profile(chit_prof)
            call gradlogchit_profile(glchit_prof)
          endif
          if (chit_aniso/=0.0) then
            call chit_aniso_profile(chit_aniso_prof)
            call gradlogchit_aniso_profile(glchit_aniso_prof)
          endif
        endif
!
!  Diffusion of the form
!
!  rho*T*Ds/Dt = ... + nab.(K*gradT)
!        Ds/Dt = ... + K/rho*[del2lnTT+(glnTT+glnhcond).glnTT]
!
!  where chix = K/(cp rho) is needed for diffus_chi calculation.
!
        chix = p%rho1*hcond*p%cp1
        glnThcond = p%glnTT + glhc/spread(hcond,2,3)    ! grad ln(T*hcond)
        call dot(p%glnTT,glnThcond,g2)
          if (pretend_lnTT) then
            thdiff = p%cv1*p%rho1*hcond * (p%del2lnTT + g2)
          else
            thdiff = p%rho1*hcond * (p%del2lnTT + g2)
          endif
      endif  ! hcond0/=0
!
!  Write out hcond z-profile (during first time step only).
!
!DM+GG This is done earlier now. The profile writing done earlier in this
! code also includes the ghost zones. This commented line may stay longer
! than the ones above.
!      if (lgravz) call write_zprof('hcond',hcond)
!
!  Write radiative flux array.
!
      if (l1davgfirst) then
        call xysum_mn_name_z(-hcond*p%TT*p%glnTT(:,3),idiag_fradz_Kprof)
        call xysum_mn_name_z(-chi_t*chit_prof*p%rho*p%TT*p%gss(:,3),idiag_fturbz)
      endif
!
!  2d-averages
!
      if (l2davgfirst) then
        if (idiag_fradxy_Kprof/=0) call zsum_mn_name_xy(-hcond*p%TT*p%glnTT(:,1),idiag_fradxy_Kprof)
        if (idiag_fturbxy/=0) call zsum_mn_name_xy(-chi_t*chit_prof*p%rho*p%TT*p%gss(:,1),idiag_fturbxy)
        if (idiag_fturbrxy/=0) &
            call zsum_mn_name_xy(-chi_t*chit_aniso_prof*chit_aniso*p%rho*p%TT* &
            (costh(m)**2*p%gss(:,1)-sinth(m)*costh(m)*p%gss(:,2)),idiag_fturbrxy)
        if (idiag_fturbthxy/=0) &
            call zsum_mn_name_xy(-chi_t*chit_aniso_prof*chit_aniso*p%rho*p%TT* &
            (-sinth(m)*costh(m)*p%gss(:,1)+sinth(m)**2*p%gss(:,2)),idiag_fturbthxy)
      endif
!
!  "Turbulent" entropy diffusion.
!
!  Should only be present if g.gradss > 0 (unstable stratification).
!  But this is not curently being checked.
!
      if (chi_t/=0.) then
        if (headtt) then
          print*,'calc_headcond: "turbulent" entropy diffusion: chi_t=',chi_t
          if (hcond0 /= 0) then
            call warning('calc_heatcond', &
                'hcond0 and chi_t combined do not seem to make sense')
          endif
        endif
!
!  ... + div(rho*T*chi*grads) = ... + chi*[del2s + (glnrho+glnTT+glnchi).grads]
!
        if (lcalc_ssmean) then
          do j=1,3; gss1(:,j)=p%gss(:,j)-gssmz(n-n1+1,j); enddo
          del2ss1=p%del2ss-del2ssmz(n-n1+1)
          call dot(p%glnrho+p%glnTT,gss1,g2)
          thdiff=thdiff+chi_t*chit_prof*(del2ss1+g2)
          call dot(glchit_prof,gss1,g2)
          thdiff=thdiff+chi_t*g2
        else
          call dot(p%glnrho+p%glnTT,p%gss,g2)
          thdiff=thdiff+chi_t*chit_prof*(p%del2ss+g2)
          call dot(glchit_prof,p%gss,g2)
          thdiff=thdiff+chi_t*g2
        endif
      endif
!
!  Turbulent entropy diffusion with rotational anisotropy:
!    chi_ij = chi_t*(delta_{ij} + chit_aniso*Om_i*Om_j),
!  where chit_aniso=chi_Om/chi_t. The first term is the isotropic part,
!  which is already dealt with above. The current formulation works only
!  in spherical coordinates. Here we assume axisymmetry so all off-diagonals
!  involving phi-indices are neglected.
!
      if (chit_aniso/=0.0 .and. lspherical_coords) then
        if (headtt) then
          print*, 'calc_headcond: '// &
              'anisotropic "turbulent" entropy diffusion: chit_aniso=', &
               chit_aniso
          if (hcond0/=0.0) then
            call warning('calc_heatcond', &
                'hcond0 and chi_t combined do not seem to make sense')
          endif
        endif
!
        sc=sinth(m)*costh(m)
        c2=costh(m)**2
        s2=sinth(m)**2
!
!  Possibility to use a simplified version where only chi_{theta r} is
!  affected by rotation (cf. Brandenburg, Moss & Tuominen 1992).
!
        if (lchit_aniso_simplified) then
!
          call g2ij(f,iss,tmp)
          thdiff=thdiff+chi_t*chit_aniso*chit_aniso_prof* &
              ((1.-3.*c2)*r1_mn*p%gss(:,1) + &
              (-sc*tmp(:,1,2))+(p%glnrho(:,2)+p%glnTT(:,2))*(-sc*p%gss(:,1)))
        else
!
!  Otherwise use the full formulation.
!
          thdiff=thdiff+chi_t*chit_aniso* &
              ((glchit_aniso_prof(:,1)*c2+chit_aniso_prof/x(l1:l2))*p%gss(:,1)+ &
              sc*(chit_aniso_prof/x(l1:l2)-glchit_prof(:,1))*p%gss(:,2))
!
          call g2ij(f,iss,tmp)
          thdiff=thdiff+chi_t*chit_aniso_prof*chit_aniso* &
              ((-sc*(tmp(:,1,2)+tmp(:,2,1))+c2*tmp(:,1,1)+s2*tmp(:,2,2))+ &
              ((p%glnrho(:,1)+p%glnTT(:,1))*(c2*p%gss(:,1)-sc*p%gss(:,2))+ &
              ( p%glnrho(:,2)+p%glnTT(:,2))*(s2*p%gss(:,2)-sc*p%gss(:,1))))
!
        endif
      endif
!
!  Check for NaNs initially.
!
      if (headt .and. (hcond0/=0.0)) then
        if (notanumber(glhc))      print*,'calc_heatcond: NaNs in glhc'
        if (notanumber(p%rho1))    print*,'calc_heatcond: NaNs in rho1'
        if (notanumber(hcond))     print*,'calc_heatcond: NaNs in hcond'
        if (notanumber(chix))      print*,'calc_heatcond: NaNs in chix'
        if (notanumber(p%del2ss))    print*,'calc_heatcond: NaNs in del2ss'
!        if (notanumber(p%del2lnrho)) print*,'calc_heatcond: NaNs in del2lnrho'
        if (notanumber(glhc))      print*,'calc_heatcond: NaNs in glhc'
        if (notanumber(1/hcond))   print*,'calc_heatcond: NaNs in 1/hcond'
        if (notanumber(p%glnTT))      print*,'calc_heatcond: NaNs in glnT'
        if (notanumber(glnThcond)) print*,'calc_heatcond: NaNs in glnThcond'
        if (notanumber(g2))        print*,'calc_heatcond: NaNs in g2'
        if (notanumber(thdiff))    print*,'calc_heatcond: NaNs in thdiff'
!
!  Most of these should trigger the following trap.
!
        if (notanumber(thdiff)) then
          print*, 'calc_heatcond: m,n,y(m),z(n)=', m, n, y(m), z(n)
          call fatal_error('calc_heatcond','NaNs in thdiff')
        endif
      endif
      if (headt .and. lfirst .and. ip == 13) then
         call output_pencil(trim(directory)//'/heatcond.dat',thdiff,1)
      endif
      if (lwrite_prof .and. ip<=9) then
        call output_pencil(trim(directory)//'/chi.dat',chix,1)
        call output_pencil(trim(directory)//'/hcond.dat',hcond,1)
        call output_pencil(trim(directory)//'/glhc.dat',glhc,3)
      endif
!
!  At the end of this routine, add all contribution to
!  thermal diffusion on the rhs of the entropy equation,
!  so Ds/Dt = ... + thdiff = ... + (...)/(rho*T)
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
!
      if (headtt) print*,'calc_heatcond: added thdiff'
!
!  Check maximum diffusion from thermal diffusion.
!  NB: With heat conduction, the second-order term for entropy is
!    gamma*chix*del2ss.
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+(gamma*chix+chi_t)*dxyz_2
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond
!***********************************************************************
    subroutine calc_heat_cool(df,p,Hmax)
!
!  Add combined heating and cooling to entropy equation.
!
!  02-jul-02/wolf: coded
!
      use Diagnostics, only: sum_mn_name
      use Gravity, only: z2
      use IO, only: output_pencil
      use Sub, only: step, cubic_step, write_zprof
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: Hmax
!
      real, dimension (nx) :: heat
      real :: profile_buffer
!
      intent(in) :: p
      intent(out) :: df
      intent(inout) :: Hmax
!
      if (pretend_lnTT) call fatal_error('calc_heat_cool', &
          'not implemented when pretend_lnTT = T')
!
!  Vertical gravity determines some heat/cool models.
!
      if (headtt) print*, 'calc_heat_cool: '// &
          'lgravz, lgravr, lgravx, lspherical_coords=', &
          lgravz, lgravr, lgravx, lspherical_coords
!
!  Initialize heating/cooling term.
!
      heat=0.0
!
!  General spatially distributed cooling profiles (independent of gravity).
!
      if (lcooling_general) call get_heat_cool_general(heat,p)
!
!  Vertical gravity case: Heat at bottom, cool top layers
!
      if (lgravz .and. ( (luminosity/=0.) .or. (cool/=0.) ) ) &
          call get_heat_cool_gravz(heat,p)
!
!  Spherical gravity case: heat at centre, cool outer layers.
!
      if (lgravr .and. (.not.lspherical_coords) .and. &
           !several possible heating/cooling sources used
           (luminosity/=0. .or. cool/=0. .or. cool_int/=0. .or. cool_ext/=0) ) &
           call get_heat_cool_gravr(heat,p)
!
!  (also see the comments inside the above subroutine to apply it to
!  spherical coordinates.)
!
!  Spherical gravity in spherical coordinate case:
!  heat at centre, cool outer layers.
!
      if (lgravx .and. lspherical_coords .and. (luminosity/=0 .or. cool/=0)) &
          call get_heat_cool_gravx_spherical(heat,p)
!
!  In Cartesian coordinates, but with the gravity in the
!  x-direction the same module may be used.
!
      if (lconvection_gravx) call get_heat_cool_gravx_cartesian(heat,p)
!
!  Add spatially uniform heating.
!
      if (heat_uniform/=0.0) heat=heat+heat_uniform
!
!  Add spatially uniform cooling in Ds/Dt equation.
!
      if (cool_uniform/=0.0) heat=heat-cool_uniform*p%rho*p%cp*p%TT
      if (cool_newton/=0.0) heat=heat-cool_newton*p%TT**4
!
!  Add cooling with constant time-scale to TTref_cool.
!
      if (tau_cool/=0.0) then
        if (.not.ltau_cool_variable) then
          heat=heat-p%rho*p%cp*gamma_inv*(p%TT-TTref_cool)/tau_cool
        else
          call calc_heat_cool_variable(heat,p)
        endif
      endif
!
!  Cooling/heating with respect to cs2mz(n)-cs2cool.
!  There is also the possibility to do cooling with respect to
!  the horizontal mean of cs2, cs2mz(n), but that has led to
!  secular instabilities for reasons that are not yet well understood.
!  A test directory is in axel/forced/BruntVaisala/tst.
!
      if (tau_cool2/=0.0) then
        if (lcalc_cs2mz_mean) then
          heat=heat-p%rho*(cs2mz(n)-cs2cool)/tau_cool2
          if (ip<12) then
            if (m==m1.and.n==50) print*,n,cs2mz(n),p%cs2(1),cs2cool
          endif
        else
          heat=heat-p%rho*(p%cs2-cs2cool)/tau_cool2
        endif
      endif
!
!  Add "coronal" heating (to simulate a hot corona).
!
      if (tau_cor>0) call get_heat_cool_corona(heat,p)
!
!  Add heating and cooling to a reference temperature in a buffer
!  zone at the z boundaries. Only regions in |z| > zheat_buffer are affected.
!  Inverse width of the transition is given by dheat_buffer1.
!
      if (tauheat_buffer/=0.) then
        profile_buffer=0.5*(1.+tanh(dheat_buffer1*(z(n)-zheat_buffer)))
        heat=heat+profile_buffer*p%ss* &
            (TTheat_buffer-1/p%TT1)/(p%rho1*tauheat_buffer)
      endif
!
!  Add heating/cooling to entropy equation.
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + p%TT1*p%rho1*heat
      if (lfirst.and.ldt) Hmax=Hmax+heat*p%rho1
!
!  Volume heating/cooling term on the mean entropy with respect to ss_const.
!
      if (lcalc_ss_volaverage.and.tau_cool_ss/=0.) then
        df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)-(ss_volaverage-ss_const)/tau_cool_ss
      endif
!
!  Heating/cooling related diagnostics.
!
      if (ldiagnos) then
        if (idiag_heatm/=0) call sum_mn_name(heat,idiag_heatm)
      endif
!
    endsubroutine calc_heat_cool
!***********************************************************************
    subroutine calc_heat_cool_variable(heat,p)
!
! Thermal relaxation for radially stratified global Keplerian disks
!
      real, dimension(nx), intent(inout) :: heat
      real, dimension (nx) :: period,rr1
      type (pencil_case), intent(in) :: p
!
      if (lcartesian_coords.or.lcylindrical_coords) then
        rr1=p%rcyl_mn1
      elseif (lspherical_coords) then
        rr1=p%r_mn1
      endif
!
      period=2*pi*rr1**(-1.5)
      heat=heat-p%rho*p%cp*gamma_inv*&
           (p%TT-TTref_cool*rr1**TT_powerlaw)/(tau_cool*period)
!
    endsubroutine calc_heat_cool_variable
!***********************************************************************
    subroutine get_heat_cool_general(heat,p)
!
      type (pencil_case) :: p
!
      real, dimension (nx) :: heat,prof
      intent(in) :: p
!
!  Subroutine to do volume heating and cooling in a layer independent of
!  gravity.
!
      select case (cooling_profile)
      case ('gaussian-z')
        prof=spread(exp(-0.5*((zcool-z(n))/wcool)**2), 1, l2-l1+1)
!
!  Cooling with a profile linear in z (unstable).
!
      case ('lin-z')
        prof=spread(z(n)/wcool,1,l2-l1+1)
      endselect
!
!  Note: the cooltype 'Temp' used below was introduced by Axel for the
!  aerosol runs. Although this 'Temp' does not match with the cooltype
!  'Temp' used in other parts of this subroutine. I have introduced
!  'Temp2' which is the same as the 'Temp' elsewhere. Later runs
!  will clarify this. - Dhruba
!
     select case (cooltype)
     case ('Temp')
       if (headtt) print*,'calc_heat_cool: cs20,cs2cool=',cs20,cs2cool
       heat=heat-cool*(p%cs2-(cs20-prof*cs2cool))/cs2cool
     case('Temp2')
       heat=heat-cool*prof*(p%cs2-cs2cool)/cs2cool
     case default
       call fatal_error('get_heat_cool_general','please select a cooltype')
     endselect
!
    endsubroutine get_heat_cool_general
!***********************************************************************
    subroutine get_heat_cool_gravz(heat,p)
!
!  Subroutine to calculate the heat/cool term for gravity along the z direction.
!
!  17-jul-07/axel: coded
!
      use Diagnostics, only: sum_mn_name, xysum_mn_name_z
      use Gravity, only: z2
      use Sub, only: step, cubic_step, write_zprof
!
      type (pencil_case) :: p
      real, dimension (nx) :: heat,prof
      real :: zbot,ztop
      intent(in) :: p
!
!  Define top and bottom positions of the box.
!
      zbot=xyz0(3)
      ztop=xyz0(3)+Lxyz(3)
!
!  Add heat near bottom (we start by default from heat=0.0)
!
!  Heating profile, normalised, so volume integral = 1.
!  Note: only the 3-D version is coded (Lx *and* Ly /= 0)
!
      if (luminosity/=0.) then
        prof = spread(exp(-0.5*((z(n)-zbot)/wheat)**2), 1, l2-l1+1) &
             /(sqrt(pi/2.)*wheat*Lx*Ly)
        heat = luminosity*prof
!
!  Smoothly switch on heating if required.
!
        if ((ttransient>0) .and. (t<ttransient)) then
           heat = heat * t*(2*ttransient-t)/ttransient**2
        endif
      endif
!
!  Allow for different cooling profile functions.
!  The gaussian default is rather broad and disturbs the entire interior.
!
      if (headtt) print*, 'cooling_profile: cooling_profile,z2,wcool,cs2cool=',cooling_profile, z2, wcool, cs2cool
      select case (cooling_profile)
      case ('gaussian')
        prof = spread(exp(-0.5*((ztop-z(n))/wcool)**2), 1, l2-l1+1)
      case ('step')
        prof = step(spread(z(n),1,nx),z2,wcool)
      case ('cubic_step')
        prof = cubic_step(spread(z(n),1,nx),z2,wcool)
!
!  Cooling with a profile linear in z (unstable).
!
      case ('lin-z')
        prof=spread(z(n)/wcool, 1, l2-l1+1)
      case default
        call fatal_error('get_heat_cool_gravz','please select a cooltype')
      endselect
      heat = heat - cool*prof*(p%cs2-cs2cool)/cs2cool
!
!  Write out cooling profile (during first time step only) and apply.
!
      call write_zprof('cooling_profile',prof)
!
!  Write divergence of cooling flux.
!
      if (l1davgfirst) then
        call xysum_mn_name_z(heat,idiag_dcoolz)
      endif
!
    endsubroutine get_heat_cool_gravz
!***********************************************************************
    subroutine get_heat_cool_gravr(heat,p)
!
!  Subroutine to calculate the heat/cool term for radial gravity
!  in cartesian and cylindrical coordinate, used in 'star-in-box' type of
!  simulations (including the sample run of geodynamo).
!  Note that this may actually work for the spherical coordinate too
!  because p%r_mn is set to be the correct quantity for each coordinate
!  system in subroutine grid. But the spherical part has not been tested
!  in spherical coordinates. At present (May 2010)
!  get_heat_cool_gravr_spherical is recommended for the spherical coordinates.
!  Normalised central heating profile so volume integral = 1
!
!  13-sep-07/boris: coded
!
      use IO, only: output_pencil
      use Sub, only: step
!
      type (pencil_case) :: p
      real, dimension (nx) :: heat,prof,theta_profile
!      real :: zbot,ztop
      intent(in) :: p
!
      if (nzgrid == 1) then
        prof = exp(-0.5*(p%r_mn/wheat)**2) * (2*pi*wheat**2)**(-1.)  ! 2-D heating profile
      else
        prof = exp(-0.5*(p%r_mn/wheat)**2) * (2*pi*wheat**2)**(-1.5) ! 3-D one
      endif
      heat = luminosity*prof
      if (headt .and. lfirst .and. ip<=9) &
           call output_pencil(trim(directory)//'/heat.dat',heat,1)
!
!  Surface cooling: entropy or temperature
!  Cooling profile; maximum = 1
!
!       prof = 0.5*(1+tanh((r_mn-1.)/wcool))
!
      if (rcool==0.) rcool=r_ext
      if (lcylindrical_coords) then
        prof = step(p%rcyl_mn,rcool,wcool)
      else
        prof = step(p%r_mn,rcool,wcool)
      endif
!
!  Pick type of cooling.
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
!  Possibility of a latitudinal heating profile.
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
          prof = 1. - step(p%r_mn,r_int,wcool)  ! inner heating/cooling step
          heat = heat - cool_int*prof*(p%cs2-cs2_int)/cs2_int*theta_profile
        else
          prof = step(p%r_mn,r_ext,wcool)     ! outer heating/cooling step
          heat = heat - cool_ext*prof*(p%cs2-cs2_ext)/cs2_ext
          prof = 1. - step(p%r_mn,r_int,wcool) ! inner heating/cooling step
          heat = heat - cool_int*prof*(p%cs2-cs2_int)/cs2_int
        endif
!
      case default
        write(unit=errormsg,fmt=*) &
             'calc_heat_cool: No such value for cooltype: ', trim(cooltype)
        call fatal_error('calc_heat_cool',errormsg)
      endselect
!
    endsubroutine get_heat_cool_gravr
!***********************************************************************
    subroutine get_heat_cool_gravx_spherical(heat,p)
!
!  Subroutine to calculate the heat/cool term in spherical coordinates
!  with gravity along x direction. At present (May 2010) this is the
!  recommended routine to use heating and cooling in spherical coordinates.
!  This is the one being used in the solar convection runs with a cooling
!  layer.
!
!   1-sep-08/dhruba: coded
!
      use Diagnostics, only: sum_mn_name, zsum_mn_name_xy
      use IO, only: output_pencil
      use Sub, only: step
!
      type (pencil_case) :: p
!
      real, dimension (nx) :: heat,prof,prof2
      real :: prof1
!
      intent(in) :: p
!
      r_ext=x(l2)
      r_int=x(l1)
!
!  Normalised central heating profile so volume integral = 1.
!
      if (nzgrid==1) then
!
!  2-D heating profile
!
        prof = exp(-0.5*(x(l1:l2)/wheat)**2) * (2*pi*wheat**2)**(-1.)
      else
        prof = exp(-0.5*(x(l1:l2)/wheat)**2) * (2*pi*wheat**2)**(-1.5)
!
!  3-D heating profile
!
      endif
      heat = luminosity*prof
      if (headt .and. lfirst .and. ip<=9) &
          call output_pencil(trim(directory)//'/heat.dat',heat,1)
!
!  Surface cooling: entropy or temperature
!  Cooling profile; maximum = 1
!
!  Pick type of cooling.
!
      select case (cooltype)
!
!  Heating/cooling at shell boundaries.
!
      case ('shell')
        if (rcool==0.0) rcool=r_ext
        prof = step(x(l1:l2),rcool,wcool)
        heat = heat - cool*prof*(p%cs2-cs2cool)/cs2cool
!
!  Latitude dependent heating/cooling: imposes a latitudinal variation
!  of temperature proportional to cos(theta) at each depth. deltaT gives
!  the amplitude of the variation between theta_0 and the equator.
!
      case ('latheat')
        if (rcool==0.0) rcool=r_ext
        prof = step(x(l1:l2),rcool1,wcool)-step(x(l1:l2),rcool2,wcool)
        prof1= 1.+deltaT*cos(2.*pi*(y(m)-y0)/Ly)
        heat = heat - cool*prof*(cs2mxy(l1:l2,m)-prof1*cs2mx(l1:l2))
!
!  Latitude dependent heating/cooling (see above) plus additional cooling
!  layer on top.
!
      case ('shell+latheat')
        if (rcool==0.0) rcool=r_ext
        prof = step(x(l1:l2),rcool1,wcool)-step(x(l1:l2),rcool2,wcool)
        prof1= 1.+deltaT*cos(2.*pi*(y(m)-y0)/Ly)
        heat = heat - cool*prof*(cs2mxy(l1:l2,m)-prof1*cs2mx(l1:l2))
!
        prof2= step(x(l1:l2),rcool,wcool)
        heat = heat - cool*prof2*(p%cs2-cs2cool)/cs2cool
!
!  Enforce latitudinal gradient of entropy plus additional cooling
!  layer on top.
!
      case ('shell+latss')
        if (rcool==0.0) rcool=r_ext
        prof = step(x(l1:l2),rcool1,wcool)-step(x(l1:l2),rcool2,wcool)
        prof1= 1.+deltaT*cos(2.*pi*(y(m)-y0)/Ly)
        heat = heat - cool*prof*(ssmxy(l1:l2,m)-prof1*ssmx(l1:l2))
!
        prof2= step(x(l1:l2),rcool,wcool)
        heat = heat - cool*prof2*(p%cs2-cs2cool)/cs2cool
!
      case default
        write(unit=errormsg,fmt=*) &
            'calc_heat_cool: No such value for cooltype: ', trim(cooltype)
        call fatal_error('calc_heat_cool',errormsg)
      endselect
!
!  Write divergence of cooling flux.
!
      if (l2davgfirst) then
        if (idiag_dcoolxy/=0) call zsum_mn_name_xy(heat,idiag_dcoolxy)
      endif
!
    endsubroutine get_heat_cool_gravx_spherical
!***********************************************************************
    subroutine get_heat_cool_gravx_cartesian(heat,p)
!
!  Subroutine to calculate the heat/cool term in cartesian coordinates
!  with gravity along x direction.
!
!  This is equivalent to the revious rutine.
!
!   4-aug-10/gustavo: coded
!
      use IO, only: output_pencil
      use Sub, only: step
!
      type (pencil_case) :: p
!
      real, dimension (nx) :: heat,prof
!
      intent(in) :: p
!
      r_ext=x(l2)
      r_int=x(l1)
!
!  Normalised central heating profile so volume integral = 1.
!
      if (nzgrid == 1) then
        prof = exp(-0.5*(x(l1:l2)/wheat)**2) * (2*pi*wheat**2)**(-1.)
!
!  2-D heating profile.
!
      else
!
!  3-D heating profile.
!
        prof = exp(-0.5*(x(l1:l2)/wheat)**2) * (2*pi*wheat**2)**(-1.5)
      endif
      heat = luminosity*prof
      if (headt .and. lfirst .and. ip<=9) &
          call output_pencil(trim(directory)//'/heat.dat',heat,1)
!
!  Surface cooling: entropy or temperature
!  Cooling profile; maximum = 1
!
!  Pick type of cooling.
!
      select case (cooltype)
!
!  Heating/cooling at shell boundaries.
!
      case ('top_layer')
        if (rcool==0.) rcool=r_ext
        prof = step(x(l1:l2),rcool,wcool)
        heat = heat - cool*prof*(p%cs2-cs2cool)/cs2cool
      case default
        write(unit=errormsg,fmt=*) &
            'calc_heat_cool: No such value for cooltype: ', trim(cooltype)
        call fatal_error('calc_heat_cool',errormsg)
      endselect
!
    endsubroutine get_heat_cool_gravx_cartesian
!***********************************************************************
    subroutine get_heat_cool_corona(heat,p)
!
!  Subroutine to calculate the heat/cool term for hot corona
!  Assume a linearly increasing reference profile.
!  This 1/rho1 business is clumpsy, but so would be obvious alternatives...
!
!  17-may-10/dhruba: coded
!
      type (pencil_case) :: p
      real, dimension (nx) :: heat
      real :: ztop,xi,profile_cor
      intent(in) :: p
!
      ztop=xyz0(3)+Lxyz(3)
      if (z(n)>=z_cor) then
        xi=(z(n)-z_cor)/(ztop-z_cor)
        profile_cor=xi**2*(3-2*xi)
        heat=heat+profile_cor*(TT_cor-1/p%TT1)/(p%rho1*tau_cor*p%cp1)
      endif
!
    endsubroutine get_heat_cool_corona
!***********************************************************************
    subroutine calc_heat_cool_RTV(df,p)
!
!  calculate cool term:  C = ne*ni*Q(T)
!  with ne*ni = 1.2*np^2 = 1.2*rho^2/(1.4*mp)^2
!  Q(T) = H*T^B is piecewice poly
!  [Q] = [v]^3 [rho] [l]^5
!
!  15-dec-04/bing: coded
!
      use IO, only: output_pencil
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: lnQ,rtv_cool,lnTT_SI,lnneni
      integer :: i,imax
      real :: unit_lnQ
!
!  Parameters for subroutine cool_RTV in SI units (from Cook et al. 1989).
!
      double precision, parameter, dimension (10) :: &
          intlnT_1 =(/4.605, 8.959, 9.906, 10.534, 11.283, 12.434, 13.286, 14.541, 17.51, 20.723 /)
      double precision, parameter, dimension (9) :: &
          lnH_1 = (/ -542.398,  -228.833, -80.245, -101.314, -78.748, -53.88, -80.452, -70.758, -91.182/), &
       B_1   = (/     50.,      15.,      0.,      2.0,      0.,    -2.,      0., -0.6667,    0.5 /)
!
!  A second set of parameters for cool_RTV (from interstellar.f90).
!
      double precision, parameter, dimension(7) :: &
          intlnT_2 = (/ 5.704,7.601 , 8.987 , 11.513 , 17.504 , 20.723, 24.0 /)
      double precision, parameter, dimension(6) :: &
          lnH_2 = (/-102.811, -99.01, -111.296, -70.804, -90.934, -80.572 /)
      double precision, parameter, dimension(6) :: &
          B_2   = (/    2.0,     1.5,   2.867,  -0.65,   0.5, 0.0 /)
!
      intent(in) :: p
      intent(out) :: df
!
      if (pretend_lnTT) call fatal_error('calc_heat_cool_RTV', &
          'not implemented when pretend_lnTT = T')
!
!     All is in SI units and has to be rescaled to PENCIL units.
!
      unit_lnQ=3*alog(real(unit_velocity))+5*alog(real(unit_length))+alog(real(unit_density))
      if (unit_system == 'cgs') unit_lnQ = unit_lnQ+alog(1.e13)
!
      lnTT_SI = p%lnTT + alog(real(unit_temperature))
!
!  Calculate ln(ne*ni) : ln(ne*ni) = ln( 1.2*rho^2/(1.4*mp)^2)
!
      lnneni = 2*p%lnrho + alog(1.2) - 2*alog(1.4*real(m_p))
!
!  First set of parameters.
!
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
!
!  Second set of parameters
!
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
!  Add to entropy equation.
!
      if (lfirst .and. ip == 13) &
           call output_pencil(trim(directory)//'/rtv.dat',rtv_cool*exp(p%lnrho+p%lnTT),1)
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss)-rtv_cool
!
      if (lfirst.and.ldt) then
         dt1_max=max(dt1_max,maxval(rtv_cool*gamma)/(cdts))
      endif
!
    endsubroutine calc_heat_cool_RTV
!***********************************************************************
    subroutine calc_tau_ss_exterior(df,p)
!
!  Entropy relaxation to zero on time scale tau_ss_exterior within
!  exterior region. For the time being this means z > zgrav.
!
!  29-jul-02/axel: coded
!
      use Gravity, only: zgrav
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real :: scl
!
      intent(in) :: p
      intent(out) :: df
!
      if (pretend_lnTT) call fatal_error('calc_tau_ss_exterior', &
          'not implemented when pretend_lnTT = T')
!
      if (headtt) print*,'calc_tau_ss_exterior: tau=',tau_ss_exterior
      if (z(n)>zgrav) then
        scl=1./tau_ss_exterior
        df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)-scl*p%ss
      endif
!
    endsubroutine calc_tau_ss_exterior
!***********************************************************************
    subroutine calc_heat_cool_prestellar(f,df,p)
!
!  Removes the heating caused by the work done to the system. Use for the
!  pre-stellar cloud simulations.
!
!  12-nov-10/mvaisala: adapted from cacl_heat_cool
!
      use IO, only: output_pencil
      use Sub, only: step, cubic_step, write_zprof
      use EquationOfState, only: eoscalc,ilnrho_ss
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: pp
      type (pencil_case) :: p
!
      intent(in) :: p
      intent(out) :: df
!
!  Add heating/cooling to entropy equation.
!
      call eoscalc(ilnrho_ss,f(l1:l2,m,n,ilnrho),f(l1:l2,m,n,iss),pp=pp)
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + 0.5*pp*p%divu
!
!  0.5 Because according to the virial theorem a collapsin cloud loses approximately
!  half of its internal energy to radiation.
!
      !df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + p%pp*p%divu
!
    endsubroutine calc_heat_cool_prestellar
!***********************************************************************
    subroutine rprint_entropy(lreset,lwrite)
!
!  Reads and registers print parameters relevant to entropy.
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Diagnostics, only: parse_name
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      integer :: iname,inamez,inamey,inamex,inamexy,inamexz,irz,inamer
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtc=0; idiag_ethm=0; idiag_ethdivum=0; idiag_ssm=0; idiag_ss2m=0
        idiag_eem=0; idiag_ppm=0; idiag_csm=0; idiag_pdivum=0; idiag_heatm=0
        idiag_ugradpm=0; idiag_ethtot=0; idiag_dtchi=0; idiag_ssmphi=0
        idiag_fradbot=0; idiag_fradtop=0; idiag_TTtop=0
        idiag_yHmax=0; idiag_yHm=0; idiag_TTmax=0; idiag_TTmin=0; idiag_TTm=0
        idiag_ssmax=0; idiag_ssmin=0
        idiag_gTrms=0; idiag_gsrms=0; idiag_gTxgsrms=0
        idiag_fconvm=0; idiag_fconvz=0; idiag_dcoolz=0; idiag_fradz=0
        idiag_fturbz=0; idiag_ppmx=0; idiag_ppmy=0; idiag_ppmz=0
        idiag_ssmx=0; idiag_ssmy=0; idiag_ssmz=0; idiag_ssmr=0; idiag_TTmr=0
        idiag_TTmx=0; idiag_TTmy=0; idiag_TTmz=0; idiag_TTmxy=0; idiag_TTmxz=0
        idiag_uxTTmz=0; idiag_uyTTmz=0; idiag_uzTTmz=0; idiag_cs2mphi=0
        idiag_ssmxy=0; idiag_ssmxz=0; idiag_fradz_Kprof=0; idiag_uxTTmxy=0
        idiag_uyTTmxy=0; idiag_uzTTmxy=0; idiag_TT2mz=0; idiag_uxTTmx=0;
        idiag_fturbxy=0; idiag_fturbrxy=0; idiag_fturbthxy=0;
        idiag_fradxy_Kprof=0; idiag_fconvxy=0;
        idiag_fradz_kramers=0; idiag_fradxy_kramers=0;
        idiag_fconvyxy=0; idiag_fconvzxy=0; idiag_dcoolxy=0
        idiag_ufpresm=0; idiag_uduum=0
      endif
!
!  iname runs through all possible names that may be listed in print.in.
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
        call parse_name(iname,cname(iname),cform(iname),'ssmin',idiag_ssmin)
        call parse_name(iname,cname(iname),cform(iname),'ssmax',idiag_ssmax)
        call parse_name(iname,cname(iname),cform(iname),'gTrms',idiag_gTrms)
        call parse_name(iname,cname(iname),cform(iname),'gsrms',idiag_gsrms)
        call parse_name(iname,cname(iname),cform(iname),'gTxgsrms',idiag_gTxgsrms)
        call parse_name(iname,cname(iname),cform(iname),'TTp',idiag_TTp)
        call parse_name(iname,cname(iname),cform(iname),'fconvm',idiag_fconvm)
        call parse_name(iname,cname(iname),cform(iname),'ufpresm',idiag_ufpresm)
        call parse_name(iname,cname(iname),cform(iname),'uduum',idiag_uduum)
      enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ssmx',idiag_ssmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ppmx',idiag_ppmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'TTmx',idiag_TTmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'uxTTmx',idiag_uxTTmx)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'ssmy',idiag_ssmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'ppmy',idiag_TTmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'TTmy',idiag_TTmy)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fturbz',idiag_fturbz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fconvz',idiag_fconvz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'dcoolz',idiag_dcoolz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fradz',idiag_fradz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fradz_Kprof',idiag_fradz_Kprof)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fradz_kramers',idiag_fradz_kramers)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ssmz',idiag_ssmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TTmz',idiag_TTmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TT2mz',idiag_TT2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ppmz',idiag_ppmz)
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
!  Check for those quantities for which we want z-averages.
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'TTmxy',idiag_TTmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'ssmxy',idiag_ssmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'uxTTmxy',idiag_uxTTmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'uyTTmxy',idiag_uyTTmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'uzTTmxy',idiag_uzTTmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'fturbxy',idiag_fturbxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'fturbrxy',idiag_fturbrxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'fturbthxy',idiag_fturbthxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'fradxy_Kprof',idiag_fradxy_Kprof)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'fradxy_kramers',idiag_fradxy_kramers)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'fconvxy',idiag_fconvxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'fconvyxy',idiag_fconvyxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'fconvzxy',idiag_fconvzxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'dcoolxy',idiag_dcoolxy)
      enddo
!
!  Check for those quantities for which we want y-averages.
!
      do inamexz=1,nnamexz
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'TTmxz',idiag_TTmxz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'ssmxz',idiag_ssmxz)
      enddo
!
!  Check for those quantities for which we want phi-averages.
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'ssmphi',idiag_ssmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'cs2mphi',idiag_cs2mphi)
      enddo
!
!  Write column where which entropy variable is stored.
!
      if (lwr) then
        write(3,*) 'nname=',nname
        write(3,*) 'iss=',iss
        write(3,*) 'iyH=',iyH
        write(3,*) 'ilnTT=',ilnTT
      endif
!
    endsubroutine rprint_entropy
!***********************************************************************
    subroutine get_slices_entropy(f,slices)
!
!  Write slices for animation of Entropy variables.
!
!  26-jul-06/tony: coded
!
      use EquationOfState, only: eoscalc, ilnrho_ss, irho_ss
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      real :: tmpval
      integer :: l,idensity,ieosvars
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Entropy.
!
        case ('ss')
          slices%yz =f(ix_loc,m1:m2,n1:n2,iss)
          slices%xz =f(l1:l2,iy_loc,n1:n2,iss)
          slices%xy =f(l1:l2,m1:m2,iz_loc,iss)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,iss)
          if (lwrite_slice_xy3) slices%xy3=f(l1:l2,m1:m2,iz3_loc,iss)
          if (lwrite_slice_xy4) slices%xy4=f(l1:l2,m1:m2,iz4_loc,iss)
          slices%ready=.true.
!
!  Pressure.
!
        case ('pp')
          if (ldensity_nolog) then
            ieosvars=irho_ss
            idensity=irho
          else
            ieosvars=ilnrho_ss
            idensity=ilnrho
          endif
!
          do m=m1,m2; do n=n1,n2
            call eoscalc(ieosvars,f(ix_loc,m,n,idensity),f(ix_loc,m,n,iss),pp=tmpval)
            slices%yz(m-m1+1,n-n1+1)=tmpval
          enddo; enddo
          do l=l1,l2; do n=n1,n2
            call eoscalc(ieosvars,f(l,iy_loc,n,idensity),f(l,iy_loc,n,iss),pp=tmpval)
            slices%xz(l-l1+1,n-n1+1)=tmpval
          enddo; enddo
          do l=l1,l2; do m=m1,m2
            call eoscalc(ieosvars,f(l,m,iz_loc,idensity),f(l,m,iz_loc,iss),pp=tmpval)
            slices%xy(l-l1+1,m-m1+1)=tmpval
            call eoscalc(ieosvars,f(l,m,iz2_loc,idensity),f(l,m,iz2_loc,iss),pp=tmpval)
            slices%xy2(l-l1+1,m-m1+1)=tmpval
            call eoscalc(ieosvars,f(l,m,iz3_loc,idensity),f(l,m,iz3_loc,iss),pp=tmpval)
            slices%xy3(l-l1+1,m-m1+1)=tmpval
            call eoscalc(ieosvars,f(l,m,iz4_loc,idensity),f(l,m,iz4_loc,iss),pp=tmpval)
            slices%xy4(l-l1+1,m-m1+1)=tmpval
          enddo; enddo
          slices%ready=.true.
!
! Temperature
!
        case ('TT')
!
          if (ldensity_nolog) then
            ieosvars=irho_ss
            idensity=irho
          else
            ieosvars=ilnrho_ss
            idensity=ilnrho
          endif
!
          do m=m1,m2; do n=n1,n2
            call eoscalc(ieosvars,f(ix_loc,m,n,idensity),f(ix_loc,m,n,iss),lnTT=tmpval)
            slices%yz(m-m1+1,n-n1+1)=exp(tmpval)
          enddo; enddo
          do l=l1,l2; do n=n1,n2
            call eoscalc(ieosvars,f(l,iy_loc,n,idensity),f(l,iy_loc,n,iss),lnTT=tmpval)
            slices%xz(l-l1+1,n-n1+1)=exp(tmpval)
          enddo; enddo
          do l=l1,l2; do m=m1,m2
            call eoscalc(ieosvars,f(l,m,iz_loc,idensity),f(l,m,iz_loc,iss),lnTT=tmpval)
            slices%xy(l-l1+1,m-m1+1)=exp(tmpval)
            call eoscalc(ieosvars,f(l,m,iz2_loc,idensity),f(l,m,iz2_loc,iss),lnTT=tmpval)
            slices%xy2(l-l1+1,m-m1+1)=exp(tmpval)
            call eoscalc(ieosvars,f(l,m,iz3_loc,idensity),f(l,m,iz3_loc,iss),lnTT=tmpval)
            slices%xy3(l-l1+1,m-m1+1)=exp(tmpval)
            call eoscalc(ieosvars,f(l,m,iz4_loc,idensity),f(l,m,iz4_loc,iss),lnTT=tmpval)
            slices%xy4(l-l1+1,m-m1+1)=exp(tmpval)
          enddo; enddo
          slices%ready=.true.
!
      endselect
!
    endsubroutine get_slices_entropy
!***********************************************************************
    subroutine calc_heatcond_zprof(zprof_hcond,zprof_glhc)
!
!  Calculate z-profile of heat conduction for multilayer setup.
!
!  12-jul-05/axel: coded
!
      use Gravity, only: z1, z2
      use Sub, only: cubic_step, cubic_der_step
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
    subroutine get_gravz_heatcond()
!
! Calculate z dependent heat conductivity and its gradient
! and stores them in the arrays which are saved at the first time
! step and used in all the next.
!
!  25-feb-2011/gustavo+dhruba: stolen from heatcond below
!
      use Gravity, only: z1, z2
      use Sub, only: step,der_step
!
      if (.not.lmultilayer) call fatal_error('get_gravz_heatcond:',&
           'dont call if you have only one layer')
!
      hcond_zprof = 1. + (hcond1-1.)*step(z,z1,-widthss) &
                     + (hcond2-1.)*step(z,z2,widthss)
      hcond_zprof = hcond0*hcond_zprof
      gradloghcond_zprof(:,1:2) = 0.0
      gradloghcond_zprof(:,3) = (hcond1-1.)*der_step(z,z1,-widthss) &
                      + (hcond2-1.)*der_step(z,z2,widthss)
      gradloghcond_zprof(:,3) = hcond0*gradloghcond_zprof(:,3)
!
    endsubroutine get_gravz_heatcond
!***********************************************************************
    subroutine get_gravx_heatcond()
!
! Calculate x dependent heat conductivity and its gradient
! and stores them in the arrays which are saved at the first time
! step and used in all the next.
!
!  10-march-2011/gustavo: Adapted from idl script strat_poly.pro
!
      use SharedVariables, only: get_shared_variable
      use Sub, only: step, der_step
      use Mpicomm, only: stop_it
      real, pointer :: cv, gravx,xc,xb
      real, dimension (:), pointer :: gravx_xpencil
      real, dimension (mx) :: dTTdxc,glhc,mpoly_xprof,dmpoly_dx
      integer :: ierr
      real :: Lum
!
      call get_shared_variable('xb',xb,ierr)
      if (ierr/=0) call stop_it(" get_gravx_heatcond: "//&
           "there was a problem when getting xb")
      call get_shared_variable('xc',xc,ierr)
      if (ierr/=0) call stop_it(" get_gravx_heatcond: "//&
           "there was a problem when getting xc")
      if (.not.lmultilayer) call fatal_error('get_gravx_heatcond:',&
           'dont call if you have only one layer')
      call get_shared_variable('cv', cv, ierr)
      if (ierr/=0) call stop_it(" get_gravx_heatcond: "//&
           "there was a problem when getting cv")
      call get_shared_variable('gravx_xpencil', gravx_xpencil, ierr)
      if (ierr/=0) call stop_it(" get_gravx_heatcond: "//&
           "there was a problem when getting gravx_xpencil")
      call get_shared_variable('gravx', gravx, ierr)
      if (ierr/=0) call stop_it(" get_gravx_heatcond: "//&
           "there was a problem when getting gravx")
!
      if (.not.lmultilayer) call fatal_error('get_gravx_heatcond:',&
           'dont call if you have only one layer')
!
! Radial profile of the polytropic index
      mpoly_xprof = mpoly0 - (mpoly0-mpoly1)*step(x,xb,widthss) - &
           (mpoly1-mpoly2)*step(x,xc,widthss)
      dmpoly_dx = -(mpoly0-mpoly1)*der_step(x,xb,widthss) - &
           (mpoly1-mpoly2)*der_step(x,xc,widthss)
!
! Hydrostatic equilibrium relations
      dTTdxc = gravx_xpencil / (cv*gamma_m1*(mpoly_xprof+1.))
      Lum = Fbot * (4.*pi*xyz0(1)**2)
!
! Kappa and its gradient are computed here
      hcond_xprof = -Lum/(4.*pi*x**2*dTTdxc)
      glhc = Lum*cv*gamma_m1/(4.*pi*gravx) * dmpoly_dx
      gradloghcond_xprof(:,1) = glhc
      gradloghcond_xprof(:,2:3) = 0.
!
      if (lroot) print*, &
           ' get_gravx_heatcond: hcond computed from gravity dependent ' //&
           ' hydrostatic equilibrium relations'
!
    endsubroutine get_gravx_heatcond
!***********************************************************************
    subroutine heatcond(hcond,p)
!
!  Calculate the heat conductivity hcond along a pencil.
!  This is an attempt to remove explicit reference to hcond[0-2] from
!  code, e.g. the boundary condition routine.
!
!  NB: if you modify this profile, you *must* adapt gradloghcond below.
!
!  23-jan-2002/wolf: coded
!  18-sep-2002/axel: added lmultilayer switch
!  09-aug-2006/dintrans: added a radial profile hcond(r)
!
      use Gravity, only: z1, z2
      use Sub, only: step
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
!  Calculate grad(log hcond), where hcond is the heat conductivity
!  NB: *Must* be in sync with heatcond() above.
!
!  23-jan-2002/wolf: coded
!
      use Gravity, only: z1, z2
      use Sub, only: der_step
!
      real, dimension (nx,3) :: glhc
      real, dimension (nx)   :: dhcond
      type (pencil_case)     :: p
!
      if (lgravz) then
        if (lmultilayer) then
          glhc(:,1:2) = 0.0
          glhc(:,3) = (hcond1-1.)*der_step(p%z_mn,z1,-widthss) &
                      + (hcond2-1.)*der_step(p%z_mn,z2,widthss)
          glhc(:,3) = hcond0*glhc(:,3)
        else
          glhc = 0.0
        endif
      else
        if (lmultilayer) then
          dhcond=(hcond1-1.)*der_step(p%r_mn,r_bcz,-widthss) &
                 + (hcond2-1.)*der_step(p%r_mn,r_ext,widthss)
          dhcond=hcond0*dhcond
          glhc(:,1) = x(l1:l2)/p%r_mn*dhcond
          glhc(:,2) = y(m)/p%r_mn*dhcond
          if (lcylinder_in_a_box) then
            glhc(:,3) = 0.0
          else
            glhc(:,3) = z(n)/p%r_mn*dhcond
          endif
        else
          glhc = 0.0
        endif
      endif
!
    endsubroutine gradloghcond
!***********************************************************************
    subroutine get_gravz_chit()
!
!  Calculate z-dependent chi_t and its gradient and store them in
!  arrays which are saved at the first time step and used in all the next.
!
!  25-feb-2011/gustavo+dhruba: stolen from chit_profile below
!
      use Gravity, only: z1, z2
      use Sub, only: step,der_step
!
      real :: zbot, ztop
!
      if (.not.lmultilayer) call fatal_error('get_gravz_chit:',&
           'dont call if you have only one layer')
!
!  If zz1 and/or zz2 are not set, use z1 and z2 instead.
!
      if (zz1 == impossible) then
          zbot=z1
      else
          zbot=zz1
      endif
!
      if (zz2 == impossible) then
          ztop=z2
      else
          ztop=zz2
      endif
!
      chit_zprof = 1 + (chit_prof1-1)*step(z,zbot,-widthss) &
                     + (chit_prof2-1)*step(z,ztop,widthss)
      gradlogchit_zprof(:,1:2) = 0.
      gradlogchit_zprof(:,3) = (chit_prof1-1)*der_step(z,zbot,-widthss) &
           + (chit_prof2-1)*der_step(z,ztop,widthss)
!
    endsubroutine get_gravz_chit
!***********************************************************************
    subroutine chit_profile(chit_prof)
!
!  Calculate the chit_profile conductivity chit_prof along a pencil.
!  This is an attempt to remove explicit reference to chit_prof[0-2] from
!  code, e.g. the boundary condition routine.
!
!  NB: if you modify this profile, you *must* adapt gradlogchit_prof below.
!
!  23-jan-2002/wolf: coded
!  18-sep-2002/axel: added lmultilayer switch
!
      use Gravity, only: z1, z2
      use Sub, only: step
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
        chit_prof = 1 + (chit_prof1-1)*step(x(l1:l2),xbot,-widthss) &
                      + (chit_prof2-1)*step(x(l1:l2),xtop,widthss)
      endif
!
    endsubroutine chit_profile
!***********************************************************************
    subroutine gradlogchit_profile(glchit_prof)
!
!  Calculate grad(log chit_prof), where chit_prof is the heat conductivity
!  NB: *Must* be in sync with heatcond() above.
!
!  23-jan-2002/wolf: coded
!
      use Gravity, only: z1, z2
      use Sub, only: der_step
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
        glchit_prof(:,1) = (chit_prof1-1)*der_step(x(l1:l2),xbot,-widthss) &
                         + (chit_prof2-1)*der_step(x(l1:l2),xtop,widthss)
        glchit_prof(:,2:3) = 0.
      endif
!
    endsubroutine gradlogchit_profile
!***********************************************************************
    subroutine chit_aniso_profile(chit_aniso_prof)
!
!  Calculate the chit_profile conductivity along a pencil.
!  This is an attempt to remove explicit reference to chit_prof[0-2] from
!  code, e.g. the boundary condition routine.
!
!  NB: if you modify this profile, you *must* adapt gradlogchit_aniso_prof
!  below.
!
!  27-apr-2011/pete: adapted from chit_profile
!
      use Sub, only: step
!
      real, dimension (nx) :: chit_aniso_prof
!
      if (lspherical_coords) then
        chit_aniso_prof = &
            1 + (chit_aniso_prof1-1)*step(x(l1:l2),xbot_aniso,-widthss) &
            + (chit_aniso_prof2-1)*step(x(l1:l2),xtop_aniso,widthss)
      endif
!
    endsubroutine chit_aniso_profile
!***********************************************************************
    subroutine gradlogchit_aniso_profile(glchit_aniso_prof)
!
!  Calculate grad(log chit_prof), where chit_prof is the heat conductivity
!  NB: *Must* be in sync with heatcond() above.
!
!  27-apr-2011/pete: adapted from gradlogchit_profile
!
      use Sub, only: der_step
!
      real, dimension (nx,3) :: glchit_aniso_prof
!
      if (lspherical_coords) then
        glchit_aniso_prof(:,1) = &
            (chit_aniso_prof1-1)*der_step(x(l1:l2),xbot_aniso,-widthss) &
            + (chit_aniso_prof2-1)*der_step(x(l1:l2),xtop_aniso,widthss)
        glchit_aniso_prof(:,2:3) = 0.
      endif
!
    endsubroutine gradlogchit_aniso_profile
!***********************************************************************
    subroutine newton_cool(df,p)
!
!  Keeps the temperature in the lower chromosphere
!  at a constant level using newton cooling.
!
!  15-dec-2004/bing: coded
!  25-sep-2006/bing: updated, using external data
!
      use EquationOfState, only: lnrho0, gamma
      use Io, only:  output_pencil
      use Mpicomm, only: mpibcast_real
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: newton
      integer, parameter :: prof_nz=150
      real, dimension (prof_nz), save :: prof_lnT,prof_z
      real :: lnTTor
      integer :: i,lend
      type (pencil_case) :: p
!
      intent(inout) :: df
      intent(in) :: p
!
      character (len=*), parameter :: lnT_dat = 'driver/b_lnT.dat'
!
      if (pretend_lnTT) call fatal_error('newton_cool', &
          'not implemented when pretend_lnTT = T')
!
!  Initial temperature profile is given in ln(T) in [K] over z in [Mm].
!
      if (it == 1) then
        if (lroot) then
          inquire(IOLENGTH=lend) lnTTor
          open (10,file=lnT_dat,form='unformatted',status='unknown',recl=lend*prof_nz)
          read (10) prof_lnT
          read (10) prof_z
          close (10)
        endif
        call mpibcast_real(prof_lnT, prof_nz)
        call mpibcast_real(prof_z, prof_nz)
!
        prof_lnT = prof_lnT - alog(real(unit_temperature))
        if (unit_system == 'SI') then
          prof_z = prof_z * 1.e6 / unit_length
        elseif (unit_system == 'cgs') then
          prof_z = prof_z * 1.e8 / unit_length
        endif
      endif
!
!  Get reference temperature.
!
      if (z(n) < prof_z(1) ) then
        lnTTor = prof_lnT(1)
      elseif (z(n) >= prof_z(prof_nz)) then
        lnTTor = prof_lnT(prof_nz)
      else
        do i=1,prof_nz-1
          if (z(n) >= prof_z(i) .and. z(n) < prof_z(i+1)) then
!
!  Linear interpolation.
!
            lnTTor = (prof_lnT(i)*(prof_z(i+1) - z(n)) +   &
                prof_lnT(i+1)*(z(n) - prof_z(i)) ) / (prof_z(i+1)-prof_z(i))
            exit
          endif
        enddo
      endif
!
      if (dt/=0.0)  then
        newton = tdown/gamma/dt*((lnTTor - p%lnTT) )
        newton = newton * exp(allp*(-abs(p%lnrho-lnrho0)))
!
!  Add newton cooling term to entropy.
!
        if (lfirst .and. ip == 13) &
            call output_pencil(trim(directory)//'/newton.dat',newton*exp(p%lnrho+p%lnTT),1)
!
        df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + newton
!
        if (lfirst.and.ldt) then
          dt1_max=max(dt1_max,maxval(abs(newton)*gamma)/(cdts))
        endif
      endif
!
    endsubroutine newton_cool
!***********************************************************************
    subroutine strat_MLT(rhotop,mixinglength_flux,lnrhom,tempm,rhobot)
!
!  04-mai-06/dintrans: called by 'mixinglength' to iterate the MLT
!  equations until rho=1 at the bottom of convection zone (z=z1)
!  see Eqs. (20-21-22) in Brandenburg et al., AN, 326 (2005)
!
!      use EquationOfState, only: gamma, gamma_m1
      use Gravity, only: z1, z2, gravz
      use Sub, only: interp1
!
      real, dimension (nzgrid) :: zz, lnrhom, tempm
      real :: rhotop, zm, ztop, dlnrho, dtemp, &
              mixinglength_flux, lnrhobot, rhobot
      real :: del, delad, fr_frac, fc_frac, fc, polyad
!
!  Initial values at the top.
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
!
      do iz=2,nzgrid
        zm=ztop-(iz-1)*dz
        zz(iz)=zm
        if (zm<=z1) then
!
!  radiative zone=polytropic stratification
!
          del=1./(mpoly1+1.)
        else
          if (zm<=z2) then
!
!  convective zone=mixing-length stratification
!
            del=delad+1.5*(fc/ &
                (exp(lnrhom(iz-1))*(gamma_m1*tempm(iz-1))**1.5))**.6666667
          else
!
!  upper zone=isothermal stratification
!
            del=0.
          endif
        endif
        dtemp=gamma*polyad*gravz*del
        dlnrho=gamma*polyad*gravz*(1.-del)/tempm(iz-1)
        tempm(iz)=tempm(iz-1)-dtemp*dz
        lnrhom(iz)=lnrhom(iz-1)-dlnrho*dz
      enddo
!
!  Interpolate.
!
      lnrhobot=interp1(zz,lnrhom,nzgrid,z1,ldescending=.true.)
      rhobot=exp(lnrhobot)
!
    endsubroutine strat_MLT
!***********************************************************************
    subroutine shell_ss_layers(f)
!
!  Initialize entropy in a spherical shell using two polytropic layers.
!
!  09-aug-06/dintrans: coded
!
      use EquationOfState, only: eoscalc, ilnrho_lnTT, mpoly0, &
                                 mpoly1, lnrho0, get_cp1
      use Gravity, only: g0
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
!  Set initial condition.
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
!
        r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
!
!  Convective layer.
!
        where (r_mn < r_ext .AND. r_mn > r_bcz)
          TT = TT_ext+beta0*(r_mn-r_ext)
          f(l1:l2,m,n,ilnrho)=lnrho0+mpoly0*log(TT)-mpoly0*log(TT_ext)
        endwhere
!
!  Radiative layer.
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
    subroutine star_heat(f)
!
!  Initialize entropy for two superposed polytropes with a central heating
!
!  20-dec-06/dintrans: coded
!  28-nov-07/dintrans: merged with strat_heat_grav
!  05-jan-10/dintrans: now the gravity field is computed in gravity_r
!
    use EquationOfState, only: rho0, lnrho0, get_soundspeed,eoscalc, ilnrho_TT
    use Sub, only: step, interp1, erfunc
!
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f
    integer, parameter   :: nr=100
    integer              :: i,l,iter
    real, dimension (nr) :: r,lnrho,temp,hcond,g,flux
    real                 :: u,r_mn,lnrho_r,temp_r,cs2,ss,lumi,Rgas
    real                 :: rhotop, rbot,rt_old,rt_new,rhobot,rb_old,rb_new,crit,r_max
!
!  Define the radial grid r=[0,r_max] and needed radial profiles.
!
      if (nzgrid == 1) then
        r_max=sqrt(xyz1(1)**2+xyz1(2)**2)
      else
        r_max=sqrt(xyz1(1)**2+xyz1(2)**2+xyz1(3)**2)
      endif
      do i=1,nr
        r(i)=r_max*float(i-1)/(nr-1)
        u=r(i)/sqrt(2.)/wheat
        if (i==1) then
          flux(1)=0.
        else
          if (nzgrid==1) then
            lumi=luminosity*(1.-exp(-u**2))
            flux(i)=lumi/(2.*pi*r(i))
          else
            lumi=luminosity*(erfunc(u)-2.*u/sqrt(pi)*exp(-u**2))
            flux(i)=lumi/(4.*pi*r(i)**2)
          endif
        endif
      enddo
      Rgas=1.-gamma_inv
      g=-flux*Rgas*(mpoly0+1.)/hcond0
!
!  The bottom density value we want at r=r_bcz, actually given by rho0.
!
    rbot=rho0
    rt_old=0.1*rbot
    rt_new=0.12*rbot
!
!  Need to iterate for rhobot=1.
!  Produce first estimate.
!
    rhotop=rt_old
    call strat_heat(nr,r,flux,g,lnrho,temp,rhotop,rhobot)
    rb_old=rhobot
!
!  Next estimate.
!
    rhotop=rt_new
    call strat_heat(nr,r,flux,g,lnrho,temp,rhotop,rhobot)
    rb_new=rhobot
!
    do 10 iter=1,10
!
!  New estimate.
!
      rhotop=rt_old+(rt_new-rt_old)/(rb_new-rb_old)*(rbot-rb_old)
!
      crit=abs(rhotop/rt_new-1.)
      if (crit<=1e-4) goto 20
      call strat_heat(nr,r,flux,g,lnrho,temp,rhotop,rhobot)
!
!  Update new estimates.
!
      rt_old=rt_new ; rb_old=rb_new
      rt_new=rhotop ; rb_new=rhobot
 10 continue
 20 print*,'- iteration completed: rhotop,crit=',rhotop,crit
!
!  Redefine rho0 and lnrho0 (important for eoscalc!).
!
      rho0=rhotop
      lnrho0=log(rhotop)
      T0=cs20/gamma_m1
      print*,'final rho0, lnrho0, T0=',rho0, lnrho0, T0
!
!  Compute rho and ss arrays from interpolations.
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
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
          call get_soundspeed(temp(i),cs2)
          call eoscalc(ilnrho_TT,lnrho(i),temp(i),ss=ss)
          write(11,'(f6.3,4e14.5,1pe14.5)') r(i),exp(lnrho(i)),ss, &
          cs2,g(i),hcond(i)
        enddo
        close(11)
      endif
!
    endsubroutine star_heat
!***********************************************************************
    subroutine strat_heat(nr,r,flux,g,lnrho,temp,rhotop,rhobot)
!
!  Compute the radial stratification for two superposed polytropic
!  layers and a central heating.
!
!  17-jan-07/dintrans: coded
!
    use EquationOfState, only: lnrho0
    use Sub, only: step, erfunc, interp1
!
    integer :: i,nr
    integer, dimension(1) :: i_ext
    real, dimension (nr) :: r,flux,g,lnrho,temp
    real :: dtemp,dlnrho,dr,rhotop,rhobot, &
            polyad,delad,fr_frac,fc_frac,fc,del,Rgas
!
!  Needed for computing a MLT stratification.
!
    polyad=1./gamma_m1
    delad=1.-gamma_inv
    fr_frac=delad*(mpoly0+1.)
    fc_frac=1.-fr_frac
    Rgas=1.-gamma_inv
!
!  Start from surface values for rho and temp.
!
    temp(nr)=cs20/gamma_m1 ; lnrho(nr)=alog(rhotop)
    dr=r(2)
    do i=nr-1,1,-1
      if (r(i+1) > r_ext) then
!
!  Isothermal exterior for r > r_ext.
!
        del=0.
      elseif (r(i+1) > r_bcz) then
!
!  Convection zone for r_bcz < r < r_ext: MLT solution if alpha_MLT/=0.
!
        fc=fc_frac*flux(i+1)
        del=delad+alpha_MLT*(fc/ &
                (exp(lnrho(i+1))*(gamma_m1*temp(i+1))**1.5))**.6666667
      else
!
!  Radiative zone for r < r_bcz.
!
        del=1./(mpoly1+1.)
      endif
      dtemp=-g(i+1)/Rgas*del
      dlnrho=-g(i+1)/(Rgas*temp(i+1))*(1.-del)
      temp(i)=temp(i+1)+dtemp*dr
      lnrho(i)=lnrho(i+1)+dlnrho*dr
    enddo
!
!  Find the value of rhobot at the bottom of convection zone.
!
!    lnrhobot=interp1(r,lnrho,nr,r_bcz)
!    rhobot=exp(lnrhobot)
!  new constraint for the iterative computation of stratification:
!  --> we impose Mtot=1 in r=(0,r_ext)
!
    i_ext=minloc(abs(r-r_ext))
    if (nzgrid == 1) then
      rhobot=sum(exp(lnrho(2:i_ext(1)-1))*r(2:i_ext(1)-1))+ &
        0.5*exp(lnrho(i_ext(1)))*r(i_ext(1))
      rhobot=rhobot*2.*pi*dr
    else
      rhobot=sum(exp(lnrho(2:i_ext(1)-1))*r(2:i_ext(1)-1)**2)+ &
        0.5*exp(lnrho(i_ext(1)))*r(i_ext(1))**2
      rhobot=rhobot*4.*pi*dr
    endif
    print*, 'total mass=', rhobot
!
    endsubroutine strat_heat
!***********************************************************************
    subroutine cylind_layers(f)
!
!  Initialise ss in a cylindrical ring using 2 superposed polytropic layers.
!
!  17-mar-07/dintrans: coded
!
      use EquationOfState, only: lnrho0, cs20, gamma, gamma_m1, cs2top, &
                                 cs2bot, get_cp1, eoscalc, ilnrho_TT
      use Gravity, only: gravz, g0
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
!  Convective layer.
!
        where (rcyl_mn <= r_ext .AND. rcyl_mn > r_bcz)
          TT=TT_ext+beta0*(rcyl_mn-r_ext)
          lnrho=lnrho0+mpoly0*log(TT/TT_ext)
        endwhere
!
!  Radiative layer.
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
!  Note: both entropy and density are initialized there (compared to layer_ss)
!
!  06-sep-07/dintrans: coded a single polytrope of index mpoly0
!
      use EquationOfState, only: eoscalc, ilnrho_TT, get_cp1, &
                                 gamma_m1, lnrho0
      use Gravity, only: gravz, g0
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: lnrho, TT, ss, z_mn
      real :: beta, cp1, zbot, ztop, TT0
      real, pointer :: gravx
      integer :: ierr
!
!  beta is the (negative) temperature gradient
!  beta = (g/cp) 1./[(1-1/gamma)*(m+1)]
!
      call get_cp1(cp1)
      if (lcylindrical_coords) then
        call get_shared_variable('gravx', gravx, ierr)
        if (ierr/=0) call fatal_error('single_polytrope', &
            'there was a problem when getting gravx')
        beta=cp1*gamma/gamma_m1*gravx/(mpoly0+1)
      else
        beta=cp1*gamma/gamma_m1*gravz/(mpoly0+1)
        ztop=xyz0(3)+Lxyz(3)
        zbot=xyz0(3)
      endif
      TT0=cs20/gamma_m1
!
!  Set initial condition (first in terms of TT, and then in terms of ss).
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
    subroutine read_hcond(hcond,glhc)
!
!  Read radial profiles of hcond and glhc from an ascii-file.
!
!  11-jun-09/pjk: coded
!
      real, dimension(nx) :: hcond
      real, dimension(nx,3) :: glhc
      integer, parameter :: ntotal=nx*nprocx
      real, dimension(nx*nprocx) :: tmp1,tmp2
      real :: var1,var2
      logical :: exist
      integer :: stat
!
!  Read hcond and glhc and write into an array.
!  If file is not found in run directory, search under trim(directory).
!
      inquire(file='hcond_glhc.dat',exist=exist)
      if (exist) then
        open(31,file='hcond_glhc.dat')
      else
        inquire(file=trim(directory)//'/hcond_glhc.ascii',exist=exist)
        if (exist) then
          open(31,file=trim(directory)//'/hcond_glhc.ascii')
        else
          call fatal_error('read_hcond','*** error *** - no input file')
        endif
      endif
!
!  Read profiles.
!
      do n=1,ntotal
        read(31,*,iostat=stat) var1,var2
        if (stat>=0) then
          if (ip<5) print*,'hcond, glhc: ',var1,var2
          tmp1(n)=var1
          tmp2(n)=var2
        else
          exit
        endif
      enddo
!
!  Assuming no ghost zones in hcond_glhc.dat.
!
      do n=l1,l2
        hcond(n-nghost)=tmp1(ipx*nx+n-nghost)
        glhc(n-nghost,1)=tmp2(ipx*nx+n-nghost)
        glhc(n-nghost,2)=0.0
        glhc(n-nghost,3)=0.0
      enddo
!
      close(31)
!
    endsubroutine read_hcond
!***********************************************************************
    subroutine fill_farray_pressure(f)
!
!  Fill f array with the pressure, to be able to calculate pressure gradient
!  directly from the pressure.
!
!  You need to open an auxiliary slot in f for this to work. Add the line
!
!  !  MAUX CONTRIBUTION 1
!
!  to the header of cparam.local.
!
!  18-feb-10/anders: coded
!
      use EquationOfState, only: eoscalc
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (mx) :: pp
!
      if (lfpres_from_pressure) then
        do n=1,mz; do m=1,my
          call eoscalc(f,mx,pp=pp)
          f(:,m,n,ippaux)=pp
        enddo; enddo
      endif
!
    endsubroutine fill_farray_pressure
!***********************************************************************
    subroutine impose_energy_floor(f)
!
!  Impose a floor in minimum entropy.  Note that entropy_floor is
!  interpreted as minimum lnTT when pretend_lnTT is set true.
!
!  07-aug-11/ccyang: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      if (entropy_floor > 0.) where(f(:,:,:,iss) < entropy_floor) f(:,:,:,ieth) = entropy_floor
!
    endsubroutine impose_energy_floor
!***********************************************************************
    subroutine dynamical_thermal_diffusion(umax)
!
!  Dummy subroutine
!
      real, intent(in) :: umax
!
      call keep_compiler_quiet(umax)
      call fatal_error('dynamical_thermal_diffusion', 'not implemented yet')
!
    endsubroutine dynamical_thermal_diffusion
!***********************************************************************
endmodule Entropy
