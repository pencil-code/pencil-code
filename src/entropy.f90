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
! PENCILS PROVIDED ugss; Ma2; fpres(3); uglnTT; sglnTT(3); transprhos !,dsdr
! PENCILS PROVIDED initss; initlnrho, uuadvec_gss 
!
!***************************************************************
module Energy
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use EquationOfState, only: gamma, gamma_m1, gamma1, cs20, cs2top, cs2bot
  use Density, only: beta_glnrho_global, beta_glnrho_scaled
  use DensityMethods, only: putrho, putlnrho, getlnrho, getrho_s
  use Messages
!
  implicit none
!
  include 'energy.h'
!
  real :: entropy_floor = impossible, TT_floor = impossible
  real, dimension(ninit) :: radius_ss=0.1, radius_ss_x=1., ampl_ss=0.0
  real :: widthss=2*epsi, epsilon_ss=0.0
  real :: luminosity=0.0, wheat=0.1, cool=0.0, cool2=0.0, zcool=0.0
  real :: rcool=0.0, wcool=0.1, ppcool=1., zcool2=0.0
  real :: rcool1=0.0, rcool2=0.0, wcool2=0.1, deltaT=0.0, cs2cool2=0.0
  real :: TT_int, TT_ext, cs2_int, cs2_ext
  real :: cool_int=0.0, cool_ext=0.0, ampl_TT=0.0
  real :: chi_jump_shock=1.0, xchi_shock=0.0,widthchi_shock=0.02
  real, target :: chi=0.0, cs2cool=0., mpoly0=1.5, mpoly1=1.5, mpoly2=1.5
  real, pointer :: mpoly
  real :: chi_t=0.0, chi_shock=0.0, chi_hyper3=0.0, chi_cspeed=0.5
  real :: chi_shock2=0.0
  real :: chi_t0=0.0, chi_t1=0.0
  real :: chi_hyper3_mesh=5.0, chi_rho=0.0
  real :: Kgperp=0.0, Kgpara=0.0, tdown=0.0, allp=2.0, TT_powerlaw=1.0
  real :: ss_left=1.0, ss_right=1.0
  real :: khor_ss=1.0, ss_const=0.0
  real :: pp_const=0.0
  real :: tau_ss_exterior=0.0, T0=0.0
  real :: mixinglength_flux=0.0, entropy_flux=0.0
  real, dimension(ninit) :: center1_x=0.0, center1_y=0.0, center1_z=0.0
  real :: center2_x=0.0, center2_y=0.0, center2_z=0.0
  real :: kx_ss=1.0, ky_ss=1.0, kz_ss=1.0
  real :: thermal_background=0.0, thermal_peak=0.0, thermal_scaling=1.0
  real :: cool_fac=1.0, chiB=0.0
  real :: downflow_cs2cool_fac=1.0
  real, dimension(3) :: chi_hyper3_aniso=0.0
  real, dimension(3) :: gradS0_imposed=(/0.0,0.0,0.0/)
  real, target :: hcond0=impossible, hcond1=impossible
  real, target :: hcondxbot=impossible, hcondxtop=0.0
  real, target :: hcondzbot=impossible, hcondztop=impossible
  real, target :: Fbot=impossible, FbotKbot=0. !set default to 0 vs impossible
  ! FbotKbot normally overwritten, but to remain finite if not
  real, target :: Ftop=impossible, FtopKtop=0. !also 0 not impossible
  real, target :: chit_prof1=1.0, chit_prof2=1.0
  real, pointer :: reduce_cs2
  real :: Kbot=impossible, Ktop=impossible
  real :: hcond2=impossible
  real :: chit_aniso=0.0, chit_aniso_prof1=1.0, chit_aniso_prof2=1.0
  real :: chit_fluct_prof1=1.0, chit_fluct_prof2=1.0
  real :: tau_cor=0.0, TT_cor=0.0, z_cor=0.0
  real :: tauheat_buffer=0.0, TTheat_buffer=0.0
  real :: heat_gaussianz=0.0, heat_gaussianz_sigma=0.0
  real :: zheat_buffer=0.0, dheat_buffer1=0.0
  real :: heat_uniform=0.0, cool_uniform=0.0, cool_newton=0.0, cool_RTV=0.0
  real :: deltaT_poleq=0.0, r_bcz=0.0
  real :: tau_cool=0.0, tau_diff=0.0, TTref_cool=0.0, tau_cool2=0.0
  real :: tau_cool_ss=0.0
  real :: cs0hs=0.0, H0hs=0.0, rho0hs=0.0
  real :: ss_volaverage=0.
  real :: xbot=0.0, xtop=0.0, alpha_MLT=1.5, xbot_aniso=0.0, xtop_aniso=0.0
  real :: xbot_chit1=0.0, xtop_chit1=0.0
  real :: zz1=impossible, zz2=impossible
  real :: zz1_fluct=impossible, zz2_fluct=impossible
  real :: rescale_TTmeanxy=1.
  real :: Pres_cutoff=impossible
  real :: pclaw=0.0, xchit=0.
  real, target :: hcond0_kramers=0.0, nkramers=0.0
  real :: chimax_kramers=0., chimin_kramers=0.
  integer :: nsmooth_kramers=0
  real :: zheat_uniform_range=0.
  real :: peh_factor=1., heat_ceiling=-1.0
  real :: Pr_smag1=1.
  real :: cs2top_ini=impossible, dcs2top_ini=impossible
  integer, parameter :: nheatc_max=4
  integer :: iglobal_hcond=0
  integer :: iglobal_glhc=0
  integer :: ippaux=0
  integer, target :: isothtop=0
  integer :: cool_type=1
  logical :: lturbulent_heat=.false.
  logical :: lheatc_Kprof=.false., lheatc_Kconst=.false., lheatc_sfluct=.false.
  logical, target :: lheatc_chiconst=.false.
  logical :: lheatc_tensordiffusion=.false., lheatc_spitzer=.false.
  logical :: lheatc_hubeny=.false.,lheatc_sqrtrhochiconst=.false.
  logical, target :: lheatc_kramers=.false.
  logical :: lheatc_smagorinsky=.false., lheatc_chit=.false.
  logical :: lheatc_corona=.false.,lheatc_chi_cspeed=.false.
  logical :: lheatc_shock=.false., lheatc_shock2=.false., lheatc_hyper3ss=.false.
  logical :: lheatc_hyper3ss_polar=.false., lheatc_hyper3ss_aniso=.false.
  logical :: lheatc_hyper3ss_mesh=.false., lheatc_shock_profr=.false.
  logical :: lcooling_general=.false.
  logical :: lupw_ss=.false.
  logical :: lcalc_ssmean=.false., lcalc_ss_volaverage=.false.
  logical :: lcalc_cs2mean=.false., lcalc_cs2mz_mean=.false.
  logical :: lcalc_cs2mz_mean_diag=.false.
  logical :: lcalc_ssmeanxy=.false.
  logical, target :: lmultilayer=.true.
  logical :: ladvection_entropy=.true.
  logical, pointer :: lpressuregradient_gas
  logical :: reinitialize_ss=.false.
  logical :: lviscosity_heat=.true.
  logical :: lfreeze_sint=.false.,lfreeze_sext=.false.
  logical :: lhcond_global=.false.,lchit_aniso_simplified=.false.
  logical :: lchit_total=.false., lchit_mean=.false., lchit_fluct=.false.
  logical :: lchiB_simplified=.false.
  logical :: lfpres_from_pressure=.false.
  logical :: lconvection_gravx=.false.
  logical :: lread_hcond=.false.
  logical :: ltau_cool_variable=.false.
  logical :: lprestellar_cool_iso=.false.
  logical :: lphotoelectric_heating=.false.
  logical :: lphotoelectric_heating_radius=.false.
  logical, pointer :: lreduced_sound_speed
  logical, pointer :: lscale_to_cs2top
  logical :: lborder_heat_variable=.false.
  logical :: lchromospheric_cooling=.false., &
             lchi_shock_density_dep=.false., &
             lhcond0_density_dep=.false.
  logical :: lenergy_slope_limited=.false.
  logical :: limpose_heat_ceiling=.false.
  logical :: lthdiff_Hmax=.false.
  logical :: lchit_noT=.false.
  logical :: lss_running_aver_as_aux=.false.
  logical :: lFenth_as_aux=.false.
  logical :: lss_flucz_as_aux=.false.
  logical :: lTT_flucz_as_aux=.false.
  logical :: lchi_t1_noprof=.false.
  real :: h_slope_limited=0.
  character (len=labellen) :: islope_limiter=''
  character (len=labellen), dimension(ninit) :: initss='nothing'
  character (len=labellen) :: borderss='nothing'
  character (len=labellen) :: pertss='zero'
  character (len=labellen) :: cooltype='Temp',cooling_profile='gaussian'
  character (len=labellen), dimension(nheatc_max) :: iheatcond='nothing'
  character (len=labellen) :: ichit='nothing'
  character (len=intlen) :: iinit_str
  real, dimension (mz)  :: ss_mz
  real, dimension (nx) :: chit_aniso_prof, dchit_aniso_prof 
  real, dimension (:), allocatable :: hcond_prof,dlnhcond_prof
  real, dimension (:), allocatable :: chit_prof_stored,chit_prof_fluct_stored
  real, dimension (:), allocatable :: dchit_prof_stored,dchit_prof_fluct_stored
!
!  xy-averaged field
!
  real, dimension (mz) :: ssmz,cs2mz
  real, dimension (nz,3) :: gssmz
  real, dimension (nz) :: del2ssmz
  real, dimension (mx) :: ssmx
  real, dimension (nx,3) :: gssmx
  real, dimension (nx) :: cs2mx, del2ssmx
  real, dimension (nx,my) :: cs2mxy, ssmxy
!
!  Input parameters.
!
  namelist /entropy_init_pars/ &
      initss, pertss, grads0, radius_ss, radius_ss_x, ampl_ss, &
      widthss, epsilon_ss, &
      mixinglength_flux, entropy_flux, &
      chi_t, chi_rho, pp_const, ss_left, ss_right, &
      ss_const, mpoly0, mpoly1, mpoly2, isothtop, khor_ss, &
!      ss_const, mpoly0, mpoly1, mpoly2, khor_ss, &
      thermal_background, thermal_peak, thermal_scaling, cs2cool, cs2cool2, &
      center1_x, center1_y, center1_z, center2_x, center2_y, center2_z, T0, &
      ampl_TT, kx_ss, ky_ss, kz_ss, beta_glnrho_global, ladvection_entropy, &
      lviscosity_heat, r_bcz, luminosity, wheat, hcond0, tau_cool, &
      tau_cool_ss, cool2, TTref_cool, lhcond_global, cool_fac, cs0hs, H0hs, &
      rho0hs, tau_cool2, lconvection_gravx, Fbot, cs2top_ini, dcs2top_ini, &
      hcond0_kramers, nkramers, alpha_MLT, lprestellar_cool_iso, lread_hcond, &
      limpose_heat_ceiling, heat_ceiling, lcooling_ss_mz, lss_running_aver_as_aux, &
      lFenth_as_aux, lss_flucz_as_aux, lTT_flucz_as_aux
!
!  Run parameters.
!
  namelist /entropy_run_pars/ &
      hcond0, hcond1, hcond2, widthss, borderss, mpoly0, mpoly1, mpoly2, &
      luminosity, wheat, cooling_profile, cooltype, cool, cs2cool, rcool, &
      rcool1, rcool2, deltaT, cs2cool2, cool2, zcool, ppcool, wcool, wcool2, Fbot, &
      lcooling_general, gradS0_imposed, &
      ss_const, chi_t, chi_rho, chit_prof1, zcool2, &
      chit_prof2, chi_shock, chi_shock2, chi, iheatcond, Kgperp, Kgpara, cool_RTV, &
      tau_ss_exterior, lmultilayer, Kbot, tau_cor, TT_cor, z_cor, &
      tauheat_buffer, TTheat_buffer, zheat_buffer, dheat_buffer1, &
      heat_gaussianz, heat_gaussianz_sigma, cs2top_ini, dcs2top_ini, &
      chi_jump_shock, xchi_shock, widthchi_shock, &
      heat_uniform, cool_uniform, cool_newton, lupw_ss, cool_int, cool_ext, &
      chi_hyper3, chi_hyper3_mesh, lturbulent_heat, deltaT_poleq, tdown, allp, &
      beta_glnrho_global, ladvection_entropy, lviscosity_heat, r_bcz, &
      lcalc_ss_volaverage, lcalc_ssmean, lcalc_cs2mean, lcalc_cs2mz_mean, &
      lfreeze_sint, lfreeze_sext, lhcond_global, tau_cool, TTref_cool, &
      mixinglength_flux, chiB, lchiB_simplified, chi_hyper3_aniso, Ftop, xbot, xtop, tau_cool2, &
      tau_cool_ss, tau_diff, lfpres_from_pressure, chit_aniso, &
      chit_aniso_prof1, chit_aniso_prof2, lchit_aniso_simplified, &
      chit_fluct_prof1, chit_fluct_prof2, &
      lconvection_gravx, ltau_cool_variable, TT_powerlaw, lcalc_ssmeanxy, &
      hcond0_kramers, nkramers, chimax_kramers, chimin_kramers, nsmooth_kramers, &
      xbot_aniso, xtop_aniso, entropy_floor, w_sldchar_ene, &
      lprestellar_cool_iso, zz1, zz2, lphotoelectric_heating, TT_floor, &
      reinitialize_ss, initss, ampl_ss, radius_ss, radius_ss_x, &
      center1_x, center1_y, center1_z, &
      lborder_heat_variable, rescale_TTmeanxy, lread_hcond,&
      Pres_cutoff,lchromospheric_cooling,lchi_shock_density_dep,lhcond0_density_dep,&
      cool_type,ichit,xchit,pclaw,lenergy_slope_limited,h_slope_limited,islope_limiter, &
      zheat_uniform_range, peh_factor, lphotoelectric_heating_radius, &
      limpose_heat_ceiling, heat_ceiling, lthdiff_Hmax, zz1_fluct, zz2_fluct, &
      Pr_smag1, chi_t0, chi_t1, lchit_total, lchit_mean, lchit_fluct, &
      chi_cspeed,xbot_chit1, xtop_chit1, lchit_noT, downflow_cs2cool_fac, &
      lss_running_aver_as_aux, lFenth_as_aux, lss_flucz_as_aux, lTT_flucz_as_aux, &
      lcalc_cs2mz_mean_diag, lchi_t1_noprof
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
  integer :: idiag_ssruzm=0     ! DIAG_DOC: $\left<s \varrho u_z/c_p\right>$
  integer :: idiag_ssuzm=0      ! DIAG_DOC: $\left<s u_z/c_p\right>$
  integer :: idiag_ssm=0        ! DIAG_DOC: $\left<s/c_p\right>$
                                ! DIAG_DOC:   \quad(mean entropy)
  integer :: idiag_ss2m=0       ! DIAG_DOC: $\left<(s/c_p)^2\right>$
                                ! DIAG_DOC:   \quad(mean squared entropy)
  integer :: idiag_eem=0        ! DIAG_DOC: $\left<e\right>$
  integer :: idiag_ppm=0        ! DIAG_DOC: $\left<p\right>$
  integer :: idiag_csm=0        ! DIAG_DOC: $\left<c_{\rm s}\right>$
  integer :: idiag_csmax=0      ! DIAG_DOC: $\max (c_{\rm s})$
  integer :: idiag_cgam=0       ! DIAG_DOC: $\left<c_{\gamma}\right>$
  integer :: idiag_pdivum=0     ! DIAG_DOC: $\left<p\nabla\cdot\uv\right>$
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
  integer :: idiag_Hmax=0       ! DIAG_DOC: $H_{\rm max}$\quad(net heat sources
                                ! DIAG_DOC: summed see \S~\ref{time-step})
  integer :: idiag_tauhmin=0    ! DIAG_DOC: $H_{\rm max}$\quad(net heat sources
                                ! DIAG_DOC: summed see \S~\ref{time-step})
  integer :: idiag_dtH=0       ! DIAG_DOC: $\delta t / [c_{\delta t,{\rm s}}\, 
                                ! DIAG_DOC:  c_{\rm v}T /H_{\rm max}]$
                                ! DIAG_DOC:   \quad(time step relative to time
                                ! DIAG_DOC:   step based on heat sources;
                                ! DIAG_DOC:   see \S~\ref{time-step})
  integer :: idiag_ssmphi=0     ! PHIAVG_DOC: $\left<s\right>_\varphi$
  integer :: idiag_cs2mphi=0    ! PHIAVG_DOC: $\left<c^2_s\right>_\varphi$
  integer :: idiag_yHm=0        ! DIAG_DOC: mean hydrogen ionization
  integer :: idiag_yHmax=0      ! DIAG_DOC: max of hydrogen ionization
  integer :: idiag_TTm=0        ! DIAG_DOC: $\left<T\right>$
  integer :: idiag_TTmax=0      ! DIAG_DOC: $T_{\max}$
  integer :: idiag_TTmin=0      ! DIAG_DOC: $T_{\min}$
  integer :: idiag_gTmax=0      ! DIAG_DOC: $\max (|\nabla T|)$
  integer :: idiag_ssmax=0      ! DIAG_DOC: $s_{\max}$
  integer :: idiag_ssmin=0      ! DIAG_DOC: $s_{\min}$
  integer :: idiag_gTrms=0      ! DIAG_DOC: $(\nabla T)_{\rm rms}$
  integer :: idiag_gsrms=0      ! DIAG_DOC: $(\nabla s)_{\rm rms}$
  integer :: idiag_gTxgsrms=0   ! DIAG_DOC: $(\nabla T\times\nabla s)_{\rm rms}$
  integer :: idiag_fconvm=0     ! DIAG_DOC: $\left< c_p \varrho u_z T \right>$
  integer :: idiag_TTp=0        ! DIAG_DOC:
  integer :: idiag_ssmr=0       ! DIAG_DOC:
  integer :: idiag_TTmr=0       ! DIAG_DOC:
  integer :: idiag_ufpresm=0    ! DIAG_DOC: $\left< -u/\rho\nabla p \right>$
  integer :: idiag_uduum=0
  integer :: idiag_Kkramersm=0  ! DIAG_DOC: $\left< K_{\rm kramers} \right>$
!
! xy averaged diagnostics given in xyaver.in
!
  integer :: idiag_fradz=0      ! XYAVG_DOC: $\left<F_{\rm rad}\right>_{xy}$
  integer :: idiag_fconvz=0     ! XYAVG_DOC: $\left<c_p \varrho u_z T \right>_{xy}$
  integer :: idiag_Fenthz=0     ! XYAVG_DOC: $\left<c_p (\varrho u_z)' T' \right>_{xy}$
  integer :: idiag_Fenthupz=0   ! XYAVG_DOC: $\left<(c_p (\varrho u_z)' T')_\uparrow \right>_{xy}$
  integer :: idiag_Fenthdownz=0 ! XYAVG_DOC: $\left<(c_p (\varrho u_z)' T')_\downarrow \right>_{xy}$
  integer :: idiag_ssmz=0       ! XYAVG_DOC: $\left< s \right>_{xy}$
  integer :: idiag_ssupmz=0     ! XYAVG_DOC: $\left< s_\uparrow \right>_{xy}$
  integer :: idiag_ssdownmz=0   ! XYAVG_DOC: $\left< s_\downarrow \right>_{xy}$
  integer :: idiag_ss2mz=0      ! XYAVG_DOC: $\left< s^2 \right>_{xy}$
  integer :: idiag_ss2upmz=0    ! XYAVG_DOC: $\left< s^2_\uparrow \right>_{xy}$
  integer :: idiag_ss2downmz=0  ! XYAVG_DOC: $\left< s^2_\downarrow \right>_{xy}$
  integer :: idiag_ssf2mz=0     ! XYAVG_DOC: $\left< s'^2\right>_{xy}$
  integer :: idiag_ssf2upmz=0   ! XYAVG_DOC: $\left< s'^2_\uparrow \right>_{xy}$
  integer :: idiag_ssf2downmz=0 ! XYAVG_DOC: $\left< s'^2_\downarrow \right>_{xy}$
  integer :: idiag_ppmz=0       ! XYAVG_DOC: $\left< p \right>_{xy}$
  integer :: idiag_TTmz=0       ! XYAVG_DOC: $\left< T \right>_{xy}$
  integer :: idiag_TTdownmz=0   ! XYAVG_DOC: $\left< T_\downarrow \right>_{xy}$
  integer :: idiag_TTupmz=0     ! XYAVG_DOC: $\left< T_\uparrow \right>_{xy}$
  integer :: idiag_TT2mz=0      ! XYAVG_DOC: $\left< T^2 \right>_{xy}$
  integer :: idiag_TT2upmz=0    ! XYAVG_DOC: $\left< T^2_\uparrow \right>_{xy}$
  integer :: idiag_TT2downmz=0  ! XYAVG_DOC: $\left< T^2_\downarrow \right>_{xy}$
  integer :: idiag_TTf2mz=0     ! XYAVG_DOC: $\left< T'^2 \right>_{xy}$
  integer :: idiag_TTf2upmz=0   ! XYAVG_DOC: $\left< T'^2_\uparrow \right>_{xy}$
  integer :: idiag_TTf2downmz=0 ! XYAVG_DOC: $\left< T'^2_\downarrow \right>_{xy}$
  integer :: idiag_uxTTmz=0     ! XYAVG_DOC: $\left< u_x T \right>_{xy}$
  integer :: idiag_uyTTmz=0     ! XYAVG_DOC: $\left< u_y T \right>_{xy}$
  integer :: idiag_uzTTmz=0     ! XYAVG_DOC: $\left< u_z T \right>_{xy}$
  integer :: idiag_uzTTupmz=0   ! XYAVG_DOC: $\left< (u_z T)_\uparrow \right>_{xy}$
  integer :: idiag_uzTTdownmz=0 ! XYAVG_DOC: $\left< (u_z T)_\downarrow \right>_{xy}$
  integer :: idiag_gTxgsxmz=0   ! XYAVG_DOC: $\left<(\nabla T\times\nabla s)_x\right>_{xy}$
  integer :: idiag_gTxgsymz=0   ! XYAVG_DOC: $\left<(\nabla T\times\nabla s)_y\right>_{xy}$
  integer :: idiag_gTxgszmz=0   ! XYAVG_DOC: $\left<(\nabla T\times\nabla s)_z\right>_{xy}$
  integer :: idiag_gTxgsx2mz=0  ! XYAVG_DOC: $\left<(\nabla T\times\nabla s)^2_x\right>_{xy}$
  integer :: idiag_gTxgsy2mz=0  ! XYAVG_DOC: $\left<(\nabla T\times\nabla s)^2_y\right>_{xy}$
  integer :: idiag_gTxgsz2mz=0  ! XYAVG_DOC: $\left<(\nabla T\times\nabla s)^2_z\right>_{xy}$
  integer :: idiag_fradz_kramers=0 ! XYAVG_DOC: $F_{\rm rad}$ (from Kramers'
                                   ! XYAVG_DOC: opacity)
  integer :: idiag_fradz_Kprof=0 ! XYAVG_DOC: $F_{\rm rad}$ (from Kprof)
  integer :: idiag_fradz_constchi=0 ! XYAVG_DOC: $F_{\rm rad}$ (from chi_const)
  integer :: idiag_fturbz=0     ! XYAVG_DOC: $\left<\varrho T \chi_t \nabla_z
                                ! XYAVG_DOC: s\right>_{xy}$ \quad(turbulent
                                ! XYAVG_DOC: heat flux)
  integer :: idiag_fturbtz=0    ! XYAVG_DOC: $\left<\varrho T \chi_t0 \nabla_z
                                ! XYAVG_DOC: s\right>_{xy}$ \quad(turbulent
                                ! XYAVG_DOC: heat flux)
  integer :: idiag_fturbmz=0    ! XYAVG_DOC: $\left<\varrho T \chi_t0 \nabla_z
                                ! XYAVG_DOC: \overline{s}\right>_{xy}$ 
                                ! XYAVG_DOC: \quad(turbulent heat flux)
  integer :: idiag_fturbfz=0    ! XYAVG_DOC: $\left<\varrho T \chi_t0 \nabla_z
                                ! XYAVG_DOC: s'\right>_{xy}$ \quad(turbulent
                                ! XYAVG_DOC: heat flux)
  integer :: idiag_dcoolz=0     ! XYAVG_DOC: surface cooling flux
  integer :: idiag_heatmz=0     ! XYAVG_DOC: heating
  integer :: idiag_Kkramersmz=0 ! XYAVG_DOC: $\left< K_0 T^{(3-b)}/\rho^{(a+1)} \right>_{xy}$
  integer :: idiag_ethmz=0      ! XYAVG_DOC: $\left<\varrho e\right>_{xy}$
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
  integer :: idiag_ss2mx=0      ! YZAVG_DOC: $\left< s^2 \right>_{yz}$
  integer :: idiag_ppmx=0       ! YZAVG_DOC: $\left< p \right>_{yz}$
  integer :: idiag_TTmx=0       ! YZAVG_DOC: $\left< T \right>_{yz}$
  integer :: idiag_TT2mx=0      ! YZAVG_DOC: $\left< T^2 \right>_{yz}$
  integer :: idiag_uxTTmx=0     ! YZAVG_DOC: $\left< u_x T \right>_{yz}$
  integer :: idiag_uyTTmx=0     ! YZAVG_DOC: $\left< u_y T \right>_{yz}$
  integer :: idiag_uzTTmx=0     ! YZAVG_DOC: $\left< u_z T \right>_{yz}$
  integer :: idiag_fconvxmx=0   ! YZAVG_DOC: $\left< c_p \varrho u_x T \right>_{yz}$
  integer :: idiag_fradmx=0     ! YZAVG_DOC: $\left<F_{\rm rad}\right>_{yz}$
  integer :: idiag_fturbmx=0    ! YZAVG_DOC: $\left<\varrho T \chi_t \nabla_x
                                ! YZAVG_DOC: s\right>_{yz}$ \quad(turbulent
                                ! YZAVG_DOC: heat flux)
  integer :: idiag_Kkramersmx=0 ! YZAVG_DOC: $\left< K_0 T^(3-b)/rho^(a+1) \right>_{yz}$
  integer :: idiag_dcoolx=0     ! YZAVG_DOC: surface cooling flux
  integer :: idiag_fradx_kramers=0 ! YZAVG_DOC: $F_{\rm rad}$ (from Kramers'
                                   ! YZAVG_DOC: opacity)
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
  integer :: idiag_gTxmxy=0    ! ZAVG_DOC: $\left<\nabla_x T\right>_{z}$
  integer :: idiag_gTymxy=0    ! ZAVG_DOC: $\left<\nabla_y T\right>_{z}$
  integer :: idiag_gTzmxy=0    ! ZAVG_DOC: $\left<\nabla_z T\right>_{z}$
  integer :: idiag_gsxmxy=0    ! ZAVG_DOC: $\left<\nabla_x s\right>_{z}$
  integer :: idiag_gsymxy=0    ! ZAVG_DOC: $\left<\nabla_y s\right>_{z}$
  integer :: idiag_gszmxy=0    ! ZAVG_DOC: $\left<\nabla_z s\right>_{z}$
  integer :: idiag_gTxgsxmxy=0 ! ZAVG_DOC: $\left<\left(\nabla T\times\nabla s\right)_x\right>_{z}$
  integer :: idiag_gTxgsymxy=0 ! ZAVG_DOC: $\left<\left(\nabla T\times\nabla s\right)_y\right>_{z}$
  integer :: idiag_gTxgszmxy=0 ! ZAVG_DOC: $\left<\left(\nabla T\times\nabla s\right)_z\right>_{z}$
  integer :: idiag_gTxgsx2mxy=0 ! ZAVG_DOC: $\left<\left(\nabla T\times\nabla s\right)_x^2\right>_{z}$
  integer :: idiag_gTxgsy2mxy=0 ! ZAVG_DOC: $\left<\left(\nabla T\times\nabla s\right)_y^2\right>_{z}$
  integer :: idiag_gTxgsz2mxy=0 ! ZAVG_DOC: $\left<\left(\nabla T\times\nabla s\right)_z^2\right>_{z}$
  integer :: idiag_fconvxy=0   ! ZAVG_DOC: $\left<c_p \varrho u_x T \right>_{z}$
  integer :: idiag_fconvyxy=0  ! ZAVG_DOC: $\left<c_p \varrho u_y T \right>_{z}$
  integer :: idiag_fconvzxy=0  ! ZAVG_DOC: $\left<c_p \varrho u_z T \right>_{z}$
  integer :: idiag_fradxy_Kprof=0 ! ZAVG_DOC: $F^{\rm rad}_x$ ($x$-component of radiative flux, $z$-averaged, from Kprof)
  integer :: idiag_fradymxy_Kprof=0 ! ZAVG_DOC: $F^{\rm rad}_y$ ($y$-component of radiative flux, $z$-averaged, from Kprof)
  integer :: idiag_fradxy_kramers=0 ! ZAVG_DOC: $F_{\rm rad}$ ($z$-averaged,
                                    ! ZAVG_DOC: from Kramers' opacity)
  integer :: idiag_fturbxy=0   ! ZAVG_DOC: $\left<\varrho T \chi_t \nabla_x
                               ! ZAVG_DOC: s\right>_{z}$
  integer :: idiag_fturbymxy=0 ! ZAVG_DOC: $\left<\varrho T \chi_t \nabla_y
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
! Auxiliaries
!
  real, dimension(:,:), pointer :: reference_state
  real, dimension (nx) :: Hmax,ss0,diffus_chi,diffus_chi3
!
  contains
!
!***********************************************************************
    subroutine register_energy
!
!  Initialise variables which should know that we solve an entropy
!  equation: iss, etc; increase nvar accordingly.
!
!  6-nov-01/wolf: coded
!
      use FArrayManager, only: farray_register_pde, farray_register_auxiliary, &
                               farray_index_append
!
      call farray_register_pde('ss',iss)
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
      if (lenergy_slope_limited) then
        if (iFF_diff==0) then
          call farray_register_auxiliary('Flux_diff',iFF_diff,vector=dimensionality)
          iFF_diff1=iFF_diff; iFF_diff2=iFF_diff+dimensionality-1
        endif
        call farray_register_auxiliary('Div_flux_diff_ss',iFF_div_ss)
        iFF_char_c=iFF_div_ss
        if (iFF_div_aa>0) iFF_char_c=max(iFF_char_c,iFF_div_aa+2)
        if (iFF_div_uu>0) iFF_char_c=max(iFF_char_c,iFF_div_uu+2)
        if (iFF_div_rho>0) iFF_char_c=max(iFF_char_c,iFF_div_rho)
      endif
!
!  Possible to calculate pressure gradient directly from the pressure, in which
!  case we must open an auxiliary slot in f to store the pressure. This method
!  is not compatible with non-blocking communication of ghost zones, so we turn
!  this behaviour off.
!
      if (lfpres_from_pressure) then
        if (ippaux==0) & 
          call farray_register_auxiliary('ppaux',ippaux)
        test_nonblocking=.true.
      endif
!
!  Heat conductivity and its gradient.
!
      if (lhcond_global) then
        if (iglobal_hcond==0) & 
          call farray_register_auxiliary('hcond',iglobal_hcond)
        if (iglobal_glhc==0) & 
          call farray_register_auxiliary('glhc',iglobal_glhc,vector=3)
      endif
!
!  Running average of entropy
!
      if (lss_running_aver_as_aux) then
        if (iss_run_aver==0) then
          call farray_register_auxiliary('ss_run_aver',iss_run_aver)
        else
          if (lroot) print*, 'register_energy: iss_run_aver = ', iss_run_aver
          call farray_index_append('iss_run_aver',iss_run_aver)
        endif
        if (lroot) write(15,*) 'ss_run_aver = fltarr(mx,my,mz)*one'
        aux_var(aux_count)=',ss_run_aver'
        if (naux+naux_com <  maux+maux_com) aux_var(aux_count)=trim(aux_var(aux_count))//' $'
        aux_count=aux_count+1
      endif
!
!  Enthalpy flux
!
      if (lFenth_as_aux) then
        if (iFenth==0) then
          call farray_register_auxiliary('Fenth',iFenth)
        else
          if (lroot) print*, 'register_energy: iFenth = ', iFenth
          call farray_index_append('iFenth',iFenth)
        endif
        if (lroot) write(15,*) 'Fenth = fltarr(mx,my,mz)*one'
        aux_var(aux_count)=',Fenth'
        if (naux+naux_com <  maux+maux_com) aux_var(aux_count)=trim(aux_var(aux_count))//' $'
        aux_count=aux_count+1
      endif
!
!  Fluctuating entropy = s - \mean_xy(s)
!
      if (lss_flucz_as_aux) then
        if (iss_flucz==0) then
          call farray_register_auxiliary('ss_flucz',iss_flucz)
        else
          if (lroot) print*, 'register_energy: iss_run_aver = ', iss_flucz
          call farray_index_append('iss_flucz',iss_flucz)
        endif
        if (lroot) write(15,*) 'ss_flucz = fltarr(mx,my,mz)*one'
        aux_var(aux_count)=',ss_flucz'
        if (naux+naux_com <  maux+maux_com) aux_var(aux_count)=trim(aux_var(aux_count))//' $'
        aux_count=aux_count+1
      endif
!
!  Fluctuating squared sound speed = TT - \mean_xy(TT)
!
      if (lTT_flucz_as_aux) then
        if (iTT_flucz==0) then
          call farray_register_auxiliary('TT_flucz',iTT_flucz)
        else
          if (lroot) print*, 'register_energy: iTT_run_aver = ', iTT_flucz
          call farray_index_append('iTT_flucz',iTT_flucz)
        endif
        if (lroot) write(15,*) 'TT_flucz = fltarr(mx,my,mz)*one'
        aux_var(aux_count)=',TT_flucz'
        if (naux+naux_com <  maux+maux_com) aux_var(aux_count)=trim(aux_var(aux_count))//' $'
        aux_count=aux_count+1
      endif
!
    endsubroutine register_energy
!***********************************************************************
    subroutine initialize_energy(f)
!
!  Called by run.f90 after reading parameters, but before the time loop.
!
!  21-jul-02/wolf: coded
!  28-mar-13/axel: reinitialize_ss added
!  26-feb-13/MR  : added temperature rescaling
!  15-nov-16/fred: option to use z-profile for reinitialize_ss
!
      use BorderProfiles, only: request_border_driving
      use EquationOfState, only: cs0, get_soundspeed, get_cp1, &
                                 select_eos_variable,gamma,gamma_m1
      use Gravity, only: gravz, g0, compute_gravity_star
      use Initcond
      use HDF5_IO, only: input_profile
      use Mpicomm, only: stop_it
      use SharedVariables, only: put_shared_variable, get_shared_variable
      use Sub, only: blob, write_zprof
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (nz) :: tmpz
      real, dimension(5) :: star_params
      real :: beta1, cp1, beta0, TT_bcz, star_cte
      integer :: i, j, q, n, m, stat
      logical :: lnothing, exist
      real, pointer :: gravx
!
!  Check any module dependencies.
!
      if (.not. leos) then
        call fatal_error('initialize_energy', &
            'EOS=noeos but entropy requires an EQUATION OF STATE for the fluid')
      endif
!
!  Logical variable lpressuregradient_gas shared with hydro modules.
!
      call get_shared_variable('lpressuregradient_gas',lpressuregradient_gas, &
                               caller='register_energy')
!
!  Tell the equation of state that we are here and what f variable we use.
!
      if (pretend_lnTT) then
        call select_eos_variable('lnTT',iss)
      else
        call select_eos_variable('ss',iss)
      endif
      if (ldensity.and..not.lstratz) then
        call get_shared_variable('mpoly',mpoly,caller='initialize_energy')
      else
        if (lroot) call warning('initialize_eos','mpoly not obtained from density,'// &
                                'set impossible')
        allocate(mpoly); mpoly=impossible
      endif
!
!  Radiative diffusion: initialize flux etc.
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
          call warning('initialize_energy', &
              'You should not set Kbot and hcond0 at the same time')
        endif
      endif
!
!  hcond0 is given by mixinglength_flux in the MLT case
!  hcond1 and hcond2 follow below in the block 'if (lmultilayer) etc...'
!
      if (mixinglength_flux/=0.0) then
        Fbot=mixinglength_flux
        hcond0=-mixinglength_flux*(mpoly0+1.)*gamma_m1*gamma1/gravz
        hcond1 = (mpoly1+1.)/(mpoly0+1.)
        hcond2 = (mpoly2+1.)/(mpoly0+1.)
        Kbot=hcond0
        lmultilayer=.true.  ! just to be sure...
        if (lroot) print*, &
            'initialize_energy: hcond0 given by mixinglength_flux=', &
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
          if (lroot) print*,'initialize_energy: cs2top,cs2cool=',cs2top,cs2cool
          if (cs2cool /= 0.) then ! cs2cool is the value to go for
            if (lroot) print*,'initialize_energy: now set cs2top=cs2cool'
            cs2top=cs2cool
          else                  ! cs2cool=0, so go for cs2top
            if (lroot) print*,'initialize_energy: now set cs2cool=cs2top'
            cs2cool=cs2top
          endif
        endif
      endif
!
      if ((lgravz.or.lgravr) .and. lrun) then
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
            if (bcz12(iss,1)=='c1') then
              Fbot=-gamma/(gamma-1)*hcond0*gravz/(mpoly0+1)
              if (lroot) print*, &
                      'initialize_energy: Calculated Fbot = ', Fbot
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
            if (bcz12(iss,2)=='c1') then
              Ftop=-gamma/(gamma-1)*hcond0*gravz/(mpoly0+1)
              if (lroot) print*, &
                      'initialize_energy: Calculated Ftop = ',Ftop
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
            if (bcz12(iss,1)=='c1') then
              Fbot=-gamma/(gamma-1)*hcond0*gravz/(mpoly+1)
              if (lroot) print*, &
                  'initialize_energy: Calculated Fbot = ', Fbot
!
              Kbot=gamma_m1/gamma*(mpoly+1.)*Fbot
              FbotKbot=gamma/gamma_m1/(mpoly+1.)
              if (lroot) print*, 'initialize_energy: Calculated Fbot,Kbot=', &
                  Fbot, Kbot
            ! else
            !! Don't need Fbot in this case (?)
            !  Fbot=-gamma/(gamma-1)*hcond0*gravz/(mpoly+1)
            !  if (lroot) print*, &
            !       'initialize_energy: Calculated Fbot = ', Fbot
            endif
          else
!
!  Define hcond0 from the given value of Fbot.
!
            if (bcz12(iss,1)=='c1') then
              hcond0=-gamma_m1/gamma*(mpoly0+1.)*Fbot/gravz
              Kbot=hcond0
              FbotKbot=-gravz*gamma/gamma_m1/(mpoly0+1.)
              if (lroot) print*, &
                   'initialize_energy: Calculated hcond0 = ', hcond0
            endif
          endif
!
!  Calculate Ftop if it has not been set in run.in.
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
!
! FbotKbot is required here for boundcond
!
      elseif (lgravx) then
        if ((coord_system=='spherical'.or.lconvection_gravx).and.(.not.lread_hcond)) then
          if (iheatcond(1)=='K-const') then
            hcondxbot=hcond0
            hcondxtop=hcond0
            FbotKbot=Fbot/hcond0
            if (lroot) &
              print*,'initialize_energy: hcondxbot, hcondxtop, FbotKbot =', &
                hcondxbot, hcondxtop, FbotKbot
          else
            if (lroot) call warning('initialize_energy', &
              'setting hcondxbot, hcondxtop, and FbotKbot is not implemented for iheatcond\=K-const')
          endif
        endif
      endif
!
!  Make sure all relevant parameters are set for spherical shell problems.
!
      select case (initss(1))
        case ('geo-kws','geo-benchmark','shell_layers')
          if (lroot) then
            print*, 'initialize_energy: '// &
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
            beta0=cp1*g0/(mpoly0+1)*gamma/gamma_m1
            beta1=cp1*g0/(mpoly1+1)*gamma/gamma_m1
            T0=cs20/gamma_m1       ! T0 defined from cs20
            TT_ext=T0
            TT_bcz=TT_ext+beta0*(1./r_bcz-1./r_ext)
            TT_int=TT_bcz+beta1*(1./r_int-1./r_bcz)
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
            endif
            TT_int=TT_ext+beta1*(1./r_int-1./r_ext)
          endif
!
!  Set up cooling parameters for spherical shell in terms of sound speeds
!
          call get_soundspeed(TT_ext,cs2_ext)
          call get_soundspeed(TT_int,cs2_int)
          cs2cool=cs2_ext
          if (lroot) then
            print*,'initialize_energy: g0,mpoly,beta1',g0,mpoly,beta1
            print*,'initialize_energy: TT_int, TT_ext=',TT_int,TT_ext
            print*,'initialize_energy: cs2_ext, cs2_int=',cs2_ext, cs2_int
          endif
!
        case ('star_heat')
          if (hcond1==impossible) hcond1=(mpoly1+1.)/(mpoly0+1.)
          if (hcond2==impossible) hcond2=(mpoly2+1.)/(mpoly0+1.)
          if (lroot) print*,'initialize_energy: set cs2cool=cs20'
          cs2cool=cs20
          cs2_ext=cs20
          if (rcool==0.) rcool=r_ext
!
!  Compute the gravity profile inside the star.
!
          star_cte=(mpoly0+1.)/hcond0*(1.-gamma1)
          call compute_gravity_star(f, wheat, luminosity, star_cte)
!
        case ('piecew_poly_cylind')
          if (bcx12(iss,1)=='c1') then
            call get_shared_variable('gravx', gravx, caller='initialize_energy')
            Fbot=-gamma/gamma_m1*hcond0*gravx/(mpoly0+1.)
            FbotKbot=-gamma/gamma_m1*gravx/(mpoly0+1.)
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
        if (lroot) print*, 'initialize_energy: Global density gradient '// &
            'with beta_glnrho_global=', beta_glnrho_global
      endif
!
!  Turn off pressure gradient term and advection for 0-D runs.
!
      if (nwgrid==1) then
        lpressuregradient_gas=.false.
        ladvection_entropy=.false.
        print*, 'initialize_energy: 0-D run, turned off pressure gradient term'
        print*, 'initialize_energy: 0-D run, turned off advection of entropy'
      endif
!
!  reinitialize entropy
!
      if (reinitialize_ss) then
        do j=1,ninit
          select case (initss(j))
          case ('blob')
            call blob(ampl_ss(j),f,iss,radius_ss(j),center1_x(j),center1_y(j),center1_z(j),radius_ss_x(j))
          case ('TTrescl')
            call rescale_TT_in_ss(f)
          case ('zprofile')
            inquire(file='zprof.txt',exist=exist)
            if (exist) then
              open(31,file='zprof.txt')
            else
              inquire(file=trim(directory)//'/zprof.ascii',exist=exist)
              if (exist) then
                open(31,file=trim(directory)//'/zprof.ascii')
              else
                call fatal_error('reinitialize_rho','error - no zprof.txt input file')
              endif
            endif
            do n=1,nzgrid
              read(31,*,iostat=stat) tmpz(n)
              if (stat<0) exit
            enddo
            do n=n1,n2
              f(:,:,n,iss)=f(:,:,n,iss)+tmpz(n-nghost+nz*ipz)
            enddo
          case default
          endselect
        enddo
      endif
!
!  Read entropy profile (used for cooling to reference profile)
!
      if (lcooling_ss_mz .and. lrun) &
        call input_profile('ss_mz', 'z', ss_mz, mz, lhas_ghost=.true.)
!
!  Initialize heat conduction.
!
      lheatc_Kprof=.false.
      lheatc_Kconst=.false.
      lheatc_sfluct=.false.
      lheatc_chiconst=.false.
      lheatc_tensordiffusion=.false.
      lheatc_spitzer=.false.
      lheatc_hubeny=.false.
      lheatc_kramers=.false.
      lheatc_corona=.false.
      lheatc_shock=.false.
      lheatc_shock2=.false.
      lheatc_shock_profr=.false.
      lheatc_hyper3ss=.false.
      lheatc_hyper3ss_polar=.false.
      lheatc_hyper3ss_mesh=.false.
      lheatc_hyper3ss_aniso=.false.
      lheatc_smagorinsky=.false.
      lheatc_chit=.false.
!
      lnothing=.false.
!
!  Select which radiative heating we are using.
!
      if (lroot) print*,'initialize_energy: nheatc_max,iheatcond=',nheatc_max,iheatcond(1:nheatc_max)
      do i=1,nheatc_max
        select case (iheatcond(i))
        case ('K-profile')
          lheatc_Kprof=.true.
          if (lroot) print*, 'heat conduction: K-profile'
        case ('K-const')
          lheatc_Kconst=.true.
          if (lroot) print*, 'heat conduction: K=cte'
        case ('sfluct')
          lheatc_sfluct=.true.
          if (lroot) print*, 'heat conduction: work purely on s fluctuations'
        case ('chi-const')
          lheatc_chiconst=.true.
          if (lroot) print*, 'heat conduction: constant chi'
        case ('chi-cspeed')
          lheatc_chi_cspeed=.true.
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
        case ('chit')
          lheatc_chit=.true.
          if (lroot) print*, 'heat conduction: chit'
        case ('corona')
          lheatc_corona=.true.
          if (lroot) print*, 'heat conduction: corona'
        case ('shock')
          lheatc_shock=.true.
          if (lroot) print*, 'heat conduction: shock'
          if (.not. lshock) &
            call stop_it('initialize_energy: shock diffusity'// &
                           ' but module setting SHOCK=noshock')
        case ('shock2')
          lheatc_shock2=.true.
          if (lroot) print*, 'heat conduction: shock'
          if (.not. lshock) &
            call stop_it('initialize_energy: shock diffusity'// &
                           ' but module setting SHOCK=noshock')
        case ('chi-shock-profr')
          lheatc_shock_profr=.true.
          if (lroot) print*, 'initialize_energy: shock with a radial profile'
          if (.not. lshock) &
            call stop_it('heat conduction: shock diffusity'// &
                           ' but module setting SHOCK=noshock')
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
        case ('smagorinsky')
          lheatc_smagorinsky=.true.
          if (lroot) print*, 'heat conduction: smagorinsky'
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
      if (hcond0==impossible) hcond0=0.

      if (lroot) then
        if (lheatc_Kprof .and. hcond0==0..and..not.lread_hcond) then
          call warning('initialize_energy', 'hcond0 is zero and no profile read in')
        endif
        if (lheatc_sfluct .and. (chi_t==0.0)) then
          call warning('initialize_energy','chi_t is zero')
        endif
        if (lheatc_chiconst .and. (chi==0.0 .and. chi_t==0.0)) then
          call warning('initialize_energy','chi and chi_t are zero')
        endif
        if (lheatc_chi_cspeed .and. (chi==0.0 .and. chi_t==0.0)) then
          call warning('initialize_energy','chi and chi_t are zero')
        endif
        if (lheatc_sqrtrhochiconst .and. (chi_rho==0.0 .and. chi_t==0.0)) then
          call warning('initialize_energy','chi_rho and chi_t are zero')
        endif
        if (chi_t==0.0.and.chit_aniso/=0.0) then
          call warning('initialize_energy','nonzero chit_aniso makes no sense with zero chi_t! Set to zero.')
          chit_aniso=0.
        endif
        if (all(iheatcond=='nothing') .and. hcond0/=0.) then
          call warning('initialize_energy', 'No heat conduction, but hcond0 /= 0')
        endif
        if (lheatc_Kconst .and. Kbot==0.0) then
          call warning('initialize_energy','Kbot is zero')
        endif
        if ((lheatc_spitzer.or.lheatc_corona) .and. (Kgpara==0.0 .or. Kgperp==0.0) ) then
          call warning('initialize_energy','Kgperp or Kgpara is zero')
        endif
        if (lheatc_kramers .and. hcond0_kramers==0.0) then
          call warning('initialize_energy','hcond0_kramers is zero')
        endif
        if (lheatc_smagorinsky .and. Pr_smag1==0.0) then
          call warning('initialize_energy','Pr_smag1 is zero')
        endif
        if (lheatc_hyper3ss .and. chi_hyper3==0.0) &
             call warning('initialize_energy','chi_hyper3 is zero')
        if (lheatc_hyper3ss_polar .and. chi_hyper3==0.0) &
             call warning('initialize_energy','chi_hyper3 is zero')
        if (lheatc_hyper3ss_mesh .and. chi_hyper3_mesh==0.0) &
             call warning('initialize_energy','chi_hyper3_mesh is zero')
        if (lheatc_shock .and. chi_shock==0.0) then
          call warning('initialize_energy','chi_shock is zero')
        endif
        if (lheatc_shock2 .and. chi_shock2==0.0) then
          call warning('initialize_energy','chi_shock2 is zero')
        endif
        if (lheatc_shock_profr .and. chi_shock==0.0) then
          call warning('initialize_energy','chi_shock is zero')
        endif
      endif

      if ( (lheatc_hyper3ss_aniso) .and.  &
           ((chi_hyper3_aniso(1)==0. .and. nxgrid/=1 ).or. &
            (chi_hyper3_aniso(2)==0. .and. nygrid/=1 ).or. &
            (chi_hyper3_aniso(3)==0. .and. nzgrid/=1 )) ) &
           call fatal_error('initialize_energy', &
           'A diffusivity coefficient of chi_hyper3 is zero')
!
!  Dynamical hyper-diffusivity operates only for mesh formulation of hyper-diffusion
!
      if (ldynamical_diffusion.and..not.lheatc_hyper3ss_mesh) &
        call fatal_error("initialize_energy",&
             "Dynamical diffusion requires mesh hyper heat conductivity, switch iheatcond='hyper3-mesh'")
!
! Compute profiles of hcond and corresponding grad(log(hcond)).
!
      if (lheatc_Kprof) then
!
        if (lread_hcond) then
          if (coord_system=='spherical' .or. lconvection_gravx) then   ! .or. lgravx ?
            call read_hcond(nx,nxgrid,ipx,hcondxtop,hcondxbot)   ! no scaling by hcond0!
          elseif (lgravz) then
            call read_hcond(nz,nzgrid,ipz,hcondztop,hcondzbot)   ! no scaling by hcond0!
          endif
        elseif (hcond0/=0.) then
          if (lgravz.and.lmultilayer) then
            call get_gravz_heatcond      ! scaling with hcond0!
            call write_zprof('hcond',hcond_prof)
            call write_zprof('gloghcond',dlnhcond_prof)
          elseif (lgravx.and.lmultilayer) then           ! i.or. coord_system=='spherical' .or. lconvection_gravx  ?
            call get_gravx_heatcond      ! no scaling by hcond0!
          elseif (lhcond_global) then
!
!  Heat conductivity stored globally when not depending on only one coordinate:
!  compute hcond and glhc and put the results in f-array.
!
!  i.e. sphere or cylinder in a box !?
!
            do n=n1,n2; do m=m1,m2
              call get_prof_pencil(f(l1:l2,m,n,iglobal_hcond),f(l1:l2,m,n,iglobal_glhc:iglobal_glhc+2), &
                                   lsphere_in_a_box.or.lcylinder_in_a_box, &
                                   hcond0,hcond1,hcond2,r_int,r_ext,llog=.true.)  ! scaling by hcond0!
            enddo; enddo
          endif

          if (chi_t/=0.) &
            call warning('initialize_energy', &
                'certain combinations of hcond0 and chi_t can be unphysical')
          if (chit_aniso/=0.) &
            call warning('initialize_energy', &
                'certain combinations of hcond0 and chit_aniso can be unphysical')
        endif

        if (chi_t/=0. .and. .not.lgravz) then
          if (chit_aniso/=0.0) call chit_aniso_profile
        endif

      endif

      if (chi_t/=0..or.chi_t0/=0.) then
        if (lheatc_chiconst.or.lheatc_Kprof.or.lheatc_kramers.or. &
            (lheatc_chit.and.(lchit_total.or.lchit_mean))) call chit_profile
      endif
!
!  Compute profiles. Diffusion of total and mean use the same profile.
!
      if (chi_t1/=0.and.lheatc_chit.and.lchit_fluct) call chit_profile_fluct
!
!  Tell the BorderProfiles module if we intend to use border driving, so
!  that the module can request the right pencils.
!
      if (borderss/='nothing') call request_border_driving(borderss)
!
!  Shared variables.
!
      call put_shared_variable('hcond0',hcond0,caller='initialize_energy')
      call put_shared_variable('hcond1',hcond1)
      call put_shared_variable('Kbot',Kbot)
      call put_shared_variable('hcondxbot',hcondxbot)
      call put_shared_variable('hcondxtop',hcondxtop)
      call put_shared_variable('hcondzbot',hcondzbot)
      call put_shared_variable('hcondztop',hcondztop)
      call put_shared_variable('Fbot',Fbot)
      call put_shared_variable('Ftop',Ftop)
      call put_shared_variable('FbotKbot',FbotKbot)
      call put_shared_variable('FtopKtop',FtopKtop)
      call put_shared_variable('chi',chi)
      call put_shared_variable('chi_t',chi_t)
      call put_shared_variable('chit_prof1',chit_prof1)
      call put_shared_variable('chit_prof2',chit_prof2)
      call put_shared_variable('lmultilayer',lmultilayer)
      call put_shared_variable('lheatc_chiconst',lheatc_chiconst)
      call put_shared_variable('lviscosity_heat',lviscosity_heat)
      call put_shared_variable('lheatc_kramers',lheatc_kramers)
      call put_shared_variable('isothtop',isothtop)
      call put_shared_variable('cs2cool',cs2cool)
      call put_shared_variable('mpoly0',mpoly0)
      call put_shared_variable('mpoly1',mpoly1)
      call put_shared_variable('mpoly2',mpoly2)
!      call put_shared_variable('lheatc_chit',lheatc_chit)

      if (lheatc_kramers) then
        call put_shared_variable('hcond0_kramers',hcond0_kramers)
        call put_shared_variable('nkramers',nkramers)
      else
        idiag_Kkramersm=0; idiag_Kkramersmx=0; idiag_Kkramersmz=0
        idiag_fradz_kramers=0; idiag_fradx_kramers=0; idiag_fradxy_kramers=0
      endif

      if (.not.(lheatc_Kprof.or.lheatc_chiconst.or.lheatc_kramers.or.lheatc_smagorinsky)) &
        idiag_fturbz=0
      if (.not.lheatc_chiconst) idiag_fradz_constchi=0
      if (.not.(lheatc_Kprof.or.lheatc_Kconst)) idiag_fradmx=0

      if (.not.lheatc_Kprof) then
        idiag_fradz_Kprof=0; idiag_fturbmx=0; idiag_fradxy_Kprof=0
        idiag_fradymxy_Kprof=0; idiag_fturbxy=0; idiag_fturbymxy=0
        idiag_fturbrxy=0; idiag_fturbthxy=0
      endif

      if (.not.lheatc_chit) then
        idiag_fturbtz=0; idiag_fturbmz=0; idiag_fturbfz=0
      endif

      if (.not.lstart) then
        if (lFenth_as_aux.and..not. lcalc_cs2mz_mean_diag) &
             call stop_it('energy_after_boundary: Need to set lcalc_cs2mz_mean_diag=T'// &
                          ' in entropy_run_pars for enthalpy diagnostics.')
!
        if ((idiag_Fenthz/=0 .or. idiag_Fenthupz/=0 .or. idiag_Fenthdownz/=0) &
            .and..not. lFenth_as_aux) then
            call stop_it('initialize_energy: Need to set lFenth_as_aux=T'// &
                         ' in entropy_run_pars for enthalpy diagnostics.')
        endif
!
        if ((idiag_ssf2mz/=0 .or. idiag_ssf2upmz/=0 .or. idiag_ssf2downmz/=0) &
            .and..not. lss_flucz_as_aux) then
            call stop_it('initialize_energy: Need to set lss_flucz_as_aux=T'// &
                         ' in entropy_run_pars for entropy fluctuation diagnostics.')
        endif
!
        if ((idiag_TTf2mz/=0 .or. idiag_TTf2upmz/=0 .or. idiag_TTf2downmz/=0) &
            .and..not. lTT_flucz_as_aux) then
            call stop_it('initialize_energy: Need to set lTT_flucz_as_aux=T'// &
                         ' in entropy_run_pars for temperature fluctuation diagnostics.')
        endif
      endif

      star_params=(/wheat,luminosity,r_bcz,widthss,alpha_MLT/)
      call put_shared_variable('star_params',star_params)
!
!  Check if reduced sound speed is used
!
      if (ldensity) then
        call get_shared_variable('lreduced_sound_speed',&
             lreduced_sound_speed)
        if (lreduced_sound_speed) then
          call get_shared_variable('reduce_cs2',reduce_cs2)
          call get_shared_variable('lscale_to_cs2top',lscale_to_cs2top)
        endif
      endif
!
      if (llocal_iso) &
           call fatal_error('initialize_energy', &
           'llocal_iso switches on the local isothermal approximation. ' // &
           'Use ENERGY=noenergy in src/Makefile.local')
!
!  Get reference_state. Requires that density is initialized before energy.
!
      if (lreference_state) &
        call get_shared_variable('reference_state',reference_state)
!
      lslope_limit_diff=lslope_limit_diff .or. lenergy_slope_limited
!
    endsubroutine initialize_energy
!***********************************************************************
    subroutine rescale_TT_in_ss(f)
!
! rescales implicitly temperature according to
! S := c_P*alog(rescale_TTmeanxy)/gamma + <S>_z + (S-<S>_z)/rescale_TTmeanxy
!
!  26-feb-13/MR: coded following Axel's recipe
!
      use EquationOfState, only: get_cp1
!
      real, dimension(mx,my,mz,mfarray), intent(INOUT) :: f
!
      real, dimension(mx,my) :: ssmxy
      real :: cp1_
      integer :: n
      real :: fac
!
      if (rescale_TTmeanxy==1.) return
      if (rescale_TTmeanxy<=0.) &
        call fatal_error('rescale_TT_in_ss', &
                         'rescale_TTmeanxy<=0')
!
      call calc_ssmeanxy(f,ssmxy)
!
      call get_cp1(cp1_)
!
      fac = alog(rescale_TTmeanxy)/(cp1_*gamma)
!
      do n=n1,n2
        f(:,:,n,iss) = (f(:,:,n,iss)-ssmxy)/rescale_TTmeanxy + ssmxy + fac
      enddo
!
    endsubroutine rescale_TT_in_ss
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
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=entropy_run_pars, IOSTAT=iostat)
!
    endsubroutine read_energy_run_pars
!***********************************************************************
    subroutine write_energy_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=entropy_run_pars)
!
    endsubroutine write_energy_run_pars
!***********************************************************************
    subroutine init_energy(f)
!
!  Initial condition for entropy.
!
!  07-nov-2001/wolf: coded
!  24-nov-2002/tony: renamed for consistancy (i.e. init_[variable name])
!  20-jan-2015/MR: changes for use of reference state
!
      use SharedVariables, only: get_shared_variable
      use EquationOfState, only: get_cp1, cs0, &
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
      intent(inout) :: f
!
      real, dimension (nx) :: tmp,pot
      real, dimension (nx) :: pp,lnrho,ss,r_mn
      real, dimension (mx) :: ss_mx
      real :: cs2int,ss0,ssint,ztop,ss_ext,pot0,pot_ext,cp1
      real, pointer :: fac_cs
      integer, pointer :: isothmid
      integer :: j,n,m
      logical :: lnothing, save_pretend_lnTT
!
!  If pretend_lnTT is set then turn it off so that initial conditions are
!  correctly generated in terms of entropy, but then restore it later
!  when we convert the data back again.
!
      save_pretend_lnTT=pretend_lnTT
      pretend_lnTT=.false.
      lnothing=.true.
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
          case ('gaussian-noise'); call gaunoise(ampl_ss(j),f,iss,iss)
          case ('blob')
            call blob(ampl_ss(j),f,iss,radius_ss(j),center1_x(j),center1_y(j),center1_z(j),radius_ss_x(j))
          case ('blob_radeq')
            call blob_radeq(ampl_ss(j),f,iss,radius_ss(j),center1_x(j),center1_y(j),center1_z(j))
          case ('isothermal'); call isothermal_entropy(f,T0)
          case ('isothermal_lnrho_ss')
            print*, 'init_energy: Isothermal density and entropy stratification'
            call isothermal_lnrho_ss(f,T0,rho0)
          case ('hydrostatic-isentropic')
            call hydrostatic_isentropic(f,lnrho_bot,ss_const)
          case ('wave')
            do n=n1,n2; do m=m1,m2
              f(l1:l2,m,n,iss)=f(l1:l2,m,n,iss)+ss_const+ampl_ss(j)*sin(kx_ss*x(l1:l2)+pi)
            enddo; enddo
          case ('wave-pressure-equil')
            call get_cp1(cp1)
            do n=n1,n2; do m=m1,m2
              tmp=ampl_ss(j)*cos(kx_ss*x(l1:l2))*cos(ky_ss*y(m))*cos(kz_ss*z(n))
              f(l1:l2,m,n,iss)=f(l1:l2,m,n,iss)+ss_const+tmp
              if (ldensity_nolog) then
                tmp=1./exp(cp1*tmp)
                f(l1:l2,m,n,irho)=f(l1:l2,m,n,irho)*tmp
                if (lreference_state) &
                  f(l1:l2,m,n,irho)=f(l1:l2,m,n,irho)-(1.-tmp)*reference_state(:,iref_rho)
              else
                f(l1:l2,m,n,ilnrho)=f(l1:l2,m,n,ilnrho)-cp1*tmp
              endif
            enddo; enddo
          case('Ferriere'); call ferriere(f)
          case('Ferriere-hs'); call ferriere_hs(f,rho0hs)
          case('Galactic-hs'); call galactic_hs(f,rho0hs,cs0hs,H0hs)
          case('xjump'); call jump(f,iss,ss_left,ss_right,widthss,'x')
          case('yjump'); call jump(f,iss,ss_left,ss_right,widthss,'y')
          case('zjump'); call jump(f,iss,ss_left,ss_right,widthss,'z')
          case('xyjump'); call jump(f,iss,ss_left,ss_right,widthss,'xy')
          case('x-y-jump'); call jump(f,iss,ss_left,ss_right,widthss,'x-y')
          case('sinxsinz'); call sinxsinz(ampl_ss(j),f,iss,kx_ss,ky_ss,kz_ss)
          case('hor-fluxtube')
            call htube(ampl_ss(j),f,iss,iss,radius_ss(j),epsilon_ss, &
                center1_x(j),center1_z(j))
          case ('hor-tube')
            call htube2(ampl_ss(j),f,iss,iss,radius_ss(j),epsilon_ss)
          case ('const_chit'); call strat_const_chit(f)
          case ('read_arr_file'); call read_outside_scal_array(f, "ss.arr", iss)
          case ('mixinglength')
             call mixinglength(mixinglength_flux,f)
             hcond0=-mixinglength_flux*(mpoly0+1.)*gamma_m1*gamma1/gravz
             print*,'init_energy : Fbot, hcond0=', Fbot, hcond0
          case ('sedov')
            if (lroot) print*,'init_energy: sedov - thermal background with gaussian energy burst'
            call blob(thermal_peak,f,iss,radius_ss(j), &
                center1_x(j),center1_y(j),center1_z(j))
          case ('sedov-dual')
            if (lroot) print*,'init_energy: sedov - thermal background with gaussian energy burst'
            call blob(thermal_peak,f,iss,radius_ss(j), &
                center1_x(j),center1_y(j),center1_z(j))
            call blob(thermal_peak,f,iss,radius_ss(j), &
                center2_x,center2_y,center2_z)
          case ('shock2d')
            if (ldensity_nolog) &
              call fatal_error('init_energy','shock2d only applicable for logarithmic density')
            call shock2d(f)
          case ('isobaric')
!
!  ss = - ln(rho/rho0)
!
            if (lroot) print*,'init_energy: isobaric stratification'
            if (pp_const==0.) then
              if (ldensity_nolog) then
                if (lreference_state) then
                  do n=n1,n2; do m=m1,m2
                    f(:,m,n,iss) = -(alog(f(:,m,n,irho)+reference_state(:,iref_rho))-lnrho0)
                  enddo; enddo
                else
                  f(:,:,:,iss) = -(alog(f(:,:,:,irho))-lnrho0)
                endif
              else
                f(:,:,:,iss) = -(f(:,:,:,ilnrho)-lnrho0)
              endif
            else
              pp=pp_const
              do n=n1,n2; do m=m1,m2
                call getlnrho(f(:,m,n,ilnrho),lnrho)
                call eoscalc(ilnrho_pp,lnrho,pp,ss=ss)
                f(l1:l2,m,n,iss)=ss
              enddo; enddo
            endif
          case ('isentropic', '1')
!
!  ss = const.
!
            if (lroot) print*,'init_energy: isentropic stratification'
            f(:,:,:,iss)=0.
            if (ampl_ss(j)/=0.) then
              print*, 'init_energy: put bubble: radius_ss(j),ampl_ss=', &
                  radius_ss(j), ampl_ss(j)
              do n=n1,n2; do m=m1,m2
                tmp=x(l1:l2)**2+y(m)**2+z(n)**2
                f(l1:l2,m,n,iss)=f(l1:l2,m,n,iss)+ &
                    ampl_ss(j)*exp(-tmp/max(radius_ss(j)**2-tmp,1e-20))
              enddo; enddo
            endif
          case ('linprof', '2')
!
!  Linear profile of ss, centered around ss=0.
!
            if (lroot) print*,'init_energy: linear entropy profile'
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
                  if (cs2cool==0.0) call fatal_error('init_energy', &
                      'inconsistency - cs2cool ca not be 0')
                  ss_ext = 0.0 + log(cs2cool/cs2_ext)
                  where (pot <= pot_ext) ! isentropic for r<r_ext
                    f(l1:l2,m,n,iss) = 0.0
                  elsewhere              ! isothermal for r>r_ext
                    f(l1:l2,m,n,iss) = ss_ext + gamma_m1*(pot-pot_ext)/cs2cool
                  endwhere
                else         ! gamma=1 --> simply isothermal (I guess [wd])
                  ! [NB: Never tested this..]
                  call getlnrho(f(:,m,n,ilnrho),lnrho)
                  f(l1:l2,m,n,iss) = -gamma_m1/gamma*(lnrho-lnrho0)
                endif
              enddo; enddo
            endif
          case ('piecew-poly', '4')
!
!  Piecewise polytropic convection setup.
!  cs0, rho0 and ss0=0 refer to height z=zref
!
            if (lroot) print*, &
                'init_energy: piecewise polytropic vertical stratification (ss)'
!
            cs2int = cs0**2
            ss0 = 0.              ! reference value ss0 is zero
            ssint = ss0
            f(:,:,:,iss) = 0.    ! just in case
!
            call get_shared_variable('isothmid', isothmid, caller='init_energy')
            call get_shared_variable('fac_cs', fac_cs)
!
!  Top layer.
            call polytropic_ss_z(f,mpoly2,zref,z2,z0+2*Lz, &
                                 isothtop,cs2int,ssint,fac_cs)
!  Unstable layer.
            call polytropic_ss_z(f,mpoly0,z2,z1,z2,isothmid,cs2int,ssint)
!  Stable layer.
            call polytropic_ss_z(f,mpoly1,z1,z0,z1,0,cs2int,ssint)

          case ('piecew-disc', '41')
!
!  Piecewise polytropic convective disc.
!  cs0, rho0 and ss0=0 refer to height z=zref.
!
            if (lroot) print*, 'init_energy: piecewise polytropic disc'
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
            if (lroot) print*,'init_energy: polytropic vertical stratification'
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
            if (lroot) print*,'init_energy: kws temperature in spherical shell'
            call shell_ss(f)
          case ('geo-benchmark')
!
!  Radial temperature profiles for spherical shell problem.
!
            if (lroot) print*, 'init_energy: '// &
                'benchmark temperature in spherical shell'
            call shell_ss(f)
          case ('shell_layers')
!
!  Radial temperature profiles for spherical shell problem.
!
            call information('init_energy', &
                'two polytropic layers in a spherical shell')
            call shell_ss_layers(f)
          case ('star_heat')
!
!  Radial temperature profiles for spherical shell problem.
!
            call information('init_energy', &
                'star_heat: now done in initial_condition_ss')
          case ('piecew_poly_cylind')
            call piecew_poly_cylind(f)
          case ('polytropic_simple')
!
!  Vertical temperature profiles for convective layer problem.
!
            call layer_ss(f)
          case ('blob_hs')
            print*, 'init_energy: put blob in hydrostatic equilibrium: '// &
                'radius_ss(j),ampl_ss(j)=', radius_ss(j), ampl_ss(j)
            call blob(ampl_ss(j),f,iss,radius_ss(j),center1_x(j),center1_y(j),center1_z(j))
            call blob(-ampl_ss(j),f,ilnrho,radius_ss(j),center1_x(j),center1_y(j),center1_z(j))
            if (ldensity_nolog) then
              f(:,:,:,irho) = exp(f(:,:,:,ilnrho))
              if (lreference_state) &
                f(l1:l2,:,:,irho) = f(l1:l2,:,:,irho) - spread(spread(reference_state(:,iref_rho),2,my),3,mz)
            endif
          case ('single_polytrope')
            call single_polytrope(f)
          case default
!
!  Catch unknown values.
!
            write(unit=errormsg,fmt=*) 'No such value for initss(' &
                //trim(iinit_str)//'): ',trim(initss(j))
            call fatal_error('init_energy',errormsg)
        endselect
!
        if (lroot) print*, 'init_energy: initss('//trim(iinit_str)//') = ', &
            trim(initss(j))
!
      enddo
!
      if (lnothing.and.lroot) print*, 'init_energy: nothing'
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
        if (lroot) print*, 'init_energy: adding hexagonal perturbation to ss'
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,iss) = f(l1:l2,m,n,iss) &
              + ampl_ss(j)*(2*cos(sqrt(3.)*0.5*khor_ss*x(l1:l2)) &
                          *cos(0.5*khor_ss*y(m)) &
                         + cos(khor_ss*y(m))) * cos(pi*z(n))
        enddo; enddo
      case default
!
!  Catch unknown values.
!
        write (unit=errormsg,fmt=*) 'No such value for pertss:', pertss
        call fatal_error('init_energy',errormsg)
!
      endselect
!
!  Interface for user's own initial condition.
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

      if (lreference_state) then    ! always meaningful?
        do n=n1,n2; do m=m1,m2
          f(:,m,n,iss) = f(:,m,n,iss) - reference_state(:,iref_s)
        enddo; enddo
      endif
!
    endsubroutine init_energy
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
         f,mpoly,zint,zbot,zblend,isoth,cs2int,ssint,fac_cs)
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
      use EquationOfState, only: get_cp1
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mz) :: stp
      real :: tmp,mpoly,zint,zbot,zblend,beta1,cs2int,ssint
      real :: cp1
      integer :: isoth
!
      real,intent(in),optional    :: fac_cs
      intent(inout) :: cs2int,ssint,f
!
      integer :: n,m
!
!  Warning: beta1 is here not dT/dz, but dcs2/dz = (gamma-1)c_p dT/dz
!
      call get_cp1(cp1)
      if (headt .and. isoth/=0.0) print*,'ssint=',ssint
      stp = step(z,zblend,widthss)
!
      do n=n1,n2; do m=m1,m2
        if (isoth/=0.0) then ! isothermal layer
          beta1 = 0.0
          tmp = ssint - gamma_m1*gravz*(z(n)-zint)/cs2int/cp1
        else
          beta1 = gamma*gravz/(mpoly+1)
          tmp = (1.0 + beta1*(z(n)-zint)/cs2int)
          ! Abort if args of log() are negative
          if ( (tmp<=0.0) .and. (z(n)<=zblend) ) then
            call fatal_error_local('polytropic_ss_z', &
                'Imaginary entropy values -- your z_inf is too low.')
          endif
          tmp = max(tmp,epsi)  ! ensure arg to log is positive
          tmp = ssint + (1-mpoly*gamma_m1)/gamma*log(tmp)/cp1
        endif
!
!  Smoothly blend the old value (above zblend) and the new one (below
!  zblend) for the two regions.
!
        f(l1:l2,m,n,iss) = stp(n)*f(l1:l2,m,n,iss)  + (1-stp(n))*tmp
!
      enddo; enddo
!
      if (isoth/=0) then
        if (present(fac_cs)) then
          ssint = ssint - gamma_m1*gravz*(zbot-zint)/cs2int/cp1+2*log(fac_cs)/gamma/cp1
        else
          ssint = ssint - gamma_m1*gravz*(zbot-zint)/cs2int/cp1
        endif
      else
        ssint = ssint + (1-mpoly*gamma_m1)/gamma &
                      * log(1 + beta1*(zbot-zint)/cs2int)/cp1
      endif
      if (isoth.ne.0 .and. present(fac_cs)) then
        cs2int = fac_cs**2*cs2int ! cs2 at layer interface (bottom)
      else
        cs2int = cs2int + beta1*(zbot-zint) ! cs2 at layer interface (bottom)
      endif
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
      integer :: n,m
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
            call fatal_error_local('polytropic_ss_disc', &
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
      call putlnrho(f(:,:,n1,ilnrho),lnrho)

      do n=n1+1,n2

        call eoscalc(ilnrho_ss,lnrho,ss_const,cs2=cs2_)
        lnrho_m=lnrho+dz*gravz/cs2_/2
        call eoscalc(ilnrho_ss,lnrho_m,ss_const,cs2=cs2_)
        lnrho=lnrho+dz*gravz/cs2_

        call putlnrho(f(:,:,n,ilnrho),lnrho)

      enddo
!
!  Entropy is simply constant.
!
      f(:,:,:,iss)=ss_const
!
    endsubroutine hydrostatic_isentropic
!***********************************************************************
    subroutine strat_const_chit(f)
!
!  Solve:
!  ds/dz=-(F/rho*chit)
!  dlnrho/dz=-ds/dz-|gravz|/cs2 using 2nd order Runge-Kutta.
!
!   2-mar-13/axel: coded
!
      use EquationOfState, only: gamma, gamma_m1, lnrho0, cs20
      use Gravity, only: gravz
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real :: zz,lnrho,ss,rho,cs2,dlnrho,dss
      integer :: iz
!
      if (headtt) print*,'init_energy : mixinglength stratification'
      if (.not.lgravz) then
        call fatal_error('mixinglength','works only for vertical gravity')
      endif
!
!  Integrate downward.
!
      if (lroot) open(11,file='data/zprof_const_chit.dat',status='unknown')
      lnrho=lnrho0
      ss=0.
      zz=xyz0(3)+Lxyz(3)
      do iz=nzgrid+nghost,0,-1
        n=iz-ipz*nz
!
!  compute right-hand side
!
        rho=exp(lnrho)
        cs2=cs20*exp(gamma_m1*(lnrho-lnrho0)+gamma*ss)
!
!  write data
!
        if (lroot) write(11,'(4(2x,1pe12.4))') zz,lnrho,ss,cs2
        if (n>=1.and.n<=mz) then
          call putlnrho(f(:,:,n,ilnrho),lnrho)
          f(:,:,n,iss)=ss
        endif
!
!  right-hand sides
!
        dss=-entropy_flux/rho
        dlnrho=-dss+gravz/cs2
!
!  integrate
!
        zz=zz-dz
        lnrho=lnrho-dlnrho*dz
        ss=ss-dss*dz
      enddo
      if (lroot) close(11)
!
    endsubroutine strat_const_chit
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
      if (headtt) print*,'init_energy : mixinglength stratification'
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
      rbot=rho0
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
      do iter=1,10
!
!  New estimate.
!
        rhotop=rt_old+(rt_new-rt_old)/(rb_new-rb_old)*(rbot-rb_old)
!
!  Check convergence.
!
        crit=abs(rhotop/rt_new-1.)
        if (crit<=1e-4) exit
!
        call strat_MLT (rhotop,mixinglength_flux,lnrhom,tempm,rhobot)
!
!  Update new estimates.
!
        rt_old=rt_new
        rb_old=rb_new
        rt_new=rhotop
        rb_new=rhobot
      enddo
      if (lfirst_proc_z) print*,'- iteration completed: rhotop,crit=',rhotop,crit
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
        if (ldensity_nolog) then
          f(:,:,n+nghost,irho)=exp(lnrho)   ! reference state?
        else
          f(:,:,n+nghost,ilnrho)=lnrho
        endif
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
      use EquationOfState, only: eoscalc, ilnrho_lnTT, get_cp1
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
          call getlnrho(f(:,m,n,ilnrho),lnrho)
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
            TT = TT_ext+beta1*(1./r_mn-1./r_ext)+pert_TT
          where (r_mn <= r_int) TT = TT_int
!
          call getlnrho(f(:,m,n,ilnrho),lnrho)
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
      use EquationOfState, only: eoscalc, ilnrho_lnTT, get_cp1
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

        call getlnrho(f(:,m,n,ilnrho),lnrho)
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
      use Mpicomm, only: mpibcast_real, MPI_COMM_WORLD
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(nx) :: absz
      double precision, dimension(nx) :: n_c,n_w,n_i,n_h
!  T in K, k_B s.t. pp is in code units ( = 9.59e-15 erg/cm/s^2)
!  (i.e. k_B = 1.381e-16 (erg/K) / 9.59e-15 (erg/cm/s^2) )
      real, parameter :: T_c_cgs=500.0,T_w_cgs=8.0e3,T_i_cgs=8.0e3,T_h_cgs=1.0e6
      real :: T_c,T_w,T_i,T_h
      real, dimension(nx) :: rho,pp
      real :: kpc,fmpi1
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
        call putrho(f(:,m,n,ilnrho),rho)
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
          fmpi1=cs2bot
          call mpibcast_real(fmpi1,0,comm=MPI_COMM_WORLD)
          cs2bot=fmpi1
          fmpi1=cs2top
          call mpibcast_real(fmpi1,ncpus-1,comm=MPI_COMM_WORLD)
          cs2top=fmpi1
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
!  20-jan-15/MR: changes for use of reference state.
!
      use EquationOfState , only: eosperturb, getmu
      use Mpicomm, only: mpibcast_real, MPI_COMM_WORLD
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(nx) :: rho,pp
      real :: rho0hs,muhs,fmpi1
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
        call fatal_error('ferriere_hs', &
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
        call putrho(f(:,m,n,ilnrho),rho)

        if (lentropy) then
!  Isothermal
          pp=rho*gamma_m1/gamma*T0
          call eosperturb(f,nx,pp=pp)
!
          fmpi1=cs2bot
          call mpibcast_real(fmpi1,0,comm=MPI_COMM_WORLD)
          cs2bot=fmpi1
          fmpi1=cs2top
          call mpibcast_real(fmpi1,ncpus-1,comm=MPI_COMM_WORLD)
          cs2top=fmpi1
!
         endif
       enddo
     enddo
!
      if (lroot) print*, 'Ferriere-hs: cs2bot=',cs2bot, ' cs2top=',cs2top
!
    endsubroutine ferriere_hs
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
!   20-jan-15/MR: changes for use of reference state.
!
      use EquationOfState, only: eosperturb
      use Mpicomm, only: mpibcast_real, MPI_COMM_WORLD
!
      real, dimension (mx,my,mz,mfarray) :: f

      real, dimension(nx) :: rho,pp,ss
      real :: rho0hs,cs0hs,H0hs,fmpi1
!
      if (lroot) print*, &
         'Galactic-hs: hydrostatic equilibrium density and entropy profiles'
!
      do n=n1,n2
      do m=m1,m2
!
        rho=rho0hs*exp(1 - sqrt(1 + (z(n)/H0hs)**2))
        call putrho(f(:,m,n,ilnrho),rho)
!
        if (lentropy) then
!
!  Isothermal
!
          pp=rho*cs0hs**2
          call eosperturb(f,nx,pp=pp)
          if (ldensity_nolog) then
            if (lreference_state) then
              ss=log(f(l1:l2,m,n,irho)+reference_state(:,iref_rho))
            else
              ss=log(f(l1:l2,m,n,irho))
            endif
          else
            ss=f(l1:l2,m,n,ilnrho)
          endif

          fmpi1=cs2bot
          call mpibcast_real(fmpi1,0,comm=MPI_COMM_WORLD)
          cs2bot=fmpi1
          fmpi1=cs2top
          call mpibcast_real(fmpi1,ncpus-1,comm=MPI_COMM_WORLD)
          cs2top=fmpi1
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
      integer :: l,ind,n,m
!
      if (lroot) print*,'shock2d: initial condition, gamma=',gamma
!
!  First quadrant:
!
        rpp(1) = 1.5
        rpr(1) = 1.5
        rpu(1) = 0.
        rpv(1) = 0.
!
!  Second quadrant:
!
        rpp(2) = 0.3
        rpr(2) = 0.532258064516129
        rpu(2) = 1.206045378311055
        rpv(2) = 0.0
!
!  Third quadrant:
!
        rpp(3) = 0.029032258064516
        rpr(3) = 0.137992831541219
        rpu(3) = 1.206045378311055
        rpv(3) = 1.206045378311055
!
!  Fourth quadrant:
!
        rpp(4) = 0.3
        rpr(4) = 0.532258064516129
        rpu(4) = 0.0
        rpv(4) = 1.206045378311055
!
!  s=-lnrho+log(gamma*p)/gamma
!
        do n=n1,n2; do m=m1,m2; do l=l1,l2
          if ( x(l)>=0.0 .and. y(m)>=0.0 ) then
            ind=1
          elseif ( x(l)<0.0 .and. y(m)>=0.0 ) then
            ind=2
          elseif ( x(l)<0.0 .and. y(m)<0.0 ) then
            ind=3
          elseif ( x(l)>=0.0 .and. y(m)<0.0 ) then
            ind=4
          else
            cycle
          endif
          if (ldensity_nolog) then
            f(l,m,n,irho)=rpr(ind)
            if (lreference_state) f(l,m,n,irho)=f(l,m,n,irho)-reference_state(l-l1+1,iref_rho)
            f(l,m,n,iss)=log(gamma*rpp(ind))/gamma-log(rpr(ind))
          else
            f(l,m,n,ilnrho)=log(rpr(ind))
            f(l,m,n,iss)=log(gamma*rpp(ind))/gamma-f(l,m,n,ilnrho)
          endif
          f(l,m,n,iux)=rpu(ind)
          f(l,m,n,iuy)=rpv(ind)
        enddo; enddo; enddo
!
    endsubroutine shock2d
!***********************************************************************
    subroutine pencil_criteria_energy
!
!  All pencils that the Entropy module depends on are specified here.
!
!  20-nov-04/anders: coded
!
      integer :: i
!
      if (lheatc_Kconst .or. lheatc_chiconst .or. lheatc_Kprof .or. &
          tau_cor>0 .or. &
          lheatc_sqrtrhochiconst) lpenc_requested(i_cp1)=.true.
      if (ldt) then
        lpenc_requested(i_cs2)=.true.
        lpenc_requested(i_ee)=.true.
      endif
      do i=1,3
        if (gradS0_imposed(i)/=0.) lpenc_requested(i_uu)=.true.
      enddo
      if (lpressuregradient_gas) lpenc_requested(i_fpres)=.true.
      if (ladvection_entropy) then
        if (lweno_transport) then
          lpenc_requested(i_rho1)=.true.
          lpenc_requested(i_transprho)=.true.
          lpenc_requested(i_transprhos)=.true.
        elseif (lfargo_advection) then
          lpenc_requested(i_uuadvec_gss)=.true.
        else
          lpenc_requested(i_ugss)=.true.
        endif
      endif
      if (lviscosity.and.lviscosity_heat) then
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_visc_heat)=.true.
        if (pretend_lnTT) lpenc_requested(i_cv1)=.true.
      endif
      if (lthdiff_Hmax) lpenc_diagnos(i_cv1)=.true.
      if (lcalc_cs2mean) lpenc_requested(i_cv1)=.true.
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
      if (cool/=0.0 .or. cool_ext/=0.0 .or. cool_int/=0.0) then
        lpenc_requested(i_cs2)=.true.
        if (cooltype=='rho_cs2') lpenc_requested(i_rho)=.true.
        if (cooltype=='two-layer') lpenc_requested(i_rho)=.true.
        if (cooltype=='corona') then
          lpenc_requested(i_cv)=.true.
           lpenc_requested(i_rho)=.true.
        endif
        if (cooling_profile=='surface_pp') lpenc_requested(i_pp)=.true.
      endif
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
      if (lheatc_sfluct) then
        lpenc_requested(i_del2ss)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_gss)=.true.
!needed?lpenc_requested(i_cp1)=.true.
!needed?lpenc_requested(i_rho1)=.true.
!needed?lpenc_requested(i_del2lnTT)=.true.
      endif
      if (lheatc_sqrtrhochiconst) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
      endif
      if (lheatc_Kprof) then
        if (hcond0/=0..or.lread_hcond) then
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
        !lpenc_requested(i_lnrho)=.true.
        !lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_cv1)=.true.
        lpenc_requested(i_del2lnTT)=.true.
      endif
      if (lheatc_chit) then
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_gss)=.true.
        lpenc_requested(i_del2ss)=.true.
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
        lpenc_requested(i_cp)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_glnrho)=.true.
        if (chi_t/=0.) then
           lpenc_requested(i_gss)=.true.
           lpenc_requested(i_del2ss)=.true.
           if (chiB/=0.0.and.(.not.lchiB_simplified)) then
             lpenc_requested(i_bij)=.true.
           endif
        endif
      endif
      if (lheatc_chi_cspeed) then
        lpenc_requested(i_cp)=.true.
        lpenc_requested(i_lnTT)=.true.
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
      if (lheatc_shock.or.lheatc_shock2) then
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_gss)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_gshock)=.true.
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_glnTT)=.true.
        if (pretend_lnTT) lpenc_requested(i_del2lnrho)=.true.
      endif
      if (lheatc_shock_profr) then
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
      if (lheatc_smagorinsky) then
        lpenc_requested(i_gss)=.true.
        lpenc_requested(i_del2ss)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_nu_smag)=.true.
      endif
      if (tau_cool/=0.0 .or. cool_uniform/=0.0) then
        lpenc_requested(i_cp)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_rho)=.true.
!
!  extra pencils for variable cooling time
!
        if (ltau_cool_variable) then
          if (lcartesian_coords.or.lcylindrical_coords) then
            lpenc_requested(i_rcyl_mn1)=.true.
          elseif (lspherical_coords) then
            lpenc_requested(i_r_mn1)=.true.
          endif
        endif
!
!  for photoelectric dust heating in debris disks
!
        if (lphotoelectric_heating) lpenc_requested(i_rhop)=.true.
        if (lphotoelectric_heating_radius) lpenc_requested(i_peh)=.true.
      endif
!
      if (tau_cool2/=0.0) lpenc_requested(i_rho)=.true.
!
! To be used to make the radiative cooling term propto exp(-P/P0), where P is the gas pressure
!
      if (Pres_cutoff/=impossible) then 
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_cs2)=.true.
      endif
        
      if (cool_newton/=0.0) then
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_cp)=.true.
      endif
!
      if (maxval(abs(beta_glnrho_scaled))/=0.0) lpenc_requested(i_cs2)=.true.
!
!  magnetic chi-quenching
!
      if (chiB/=0.0) lpenc_requested(i_b2)=.true.
!
      if (lfargo_advection) then
        lpenc_requested(i_uu_advec)=.true.
        lpenc_requested(i_gss)=.true.
        lpenc_requested(i_uuadvec_gss)=.true.
      endif
!      
      if (borderss == 'initial-temperature') lpenc_requested(i_lnrho)=.true.
!
      lpenc_diagnos2d(i_ss)=.true.
!
      if (idiag_dtchi/=0) lpenc_diagnos(i_rho1)=.true.
      if (idiag_dtH/=0) lpenc_diagnos(i_ee)=.true.
      if (idiag_Hmax/=0) lpenc_diagnos(i_ee)=.true.
      if (idiag_tauhmin/=0) lpenc_diagnos(i_cv1)=.true.
      if (idiag_ethdivum/=0) lpenc_diagnos(i_divu)=.true.
      if (idiag_cgam/=0) then
        lpenc_diagnos(i_TT)=.true.
        lpenc_diagnos(i_cp1)=.true.
        lpenc_diagnos(i_rho1)=.true.
      endif
      if (idiag_csm/=0 .or. idiag_csmax/=0) lpenc_diagnos(i_cs2)=.true.
      if (idiag_eem/=0) lpenc_diagnos(i_ee)=.true.
      if (idiag_ppm/=0) lpenc_diagnos(i_pp)=.true.
      if (idiag_pdivum/=0) then
        lpenc_diagnos(i_pp)=.true.
        lpenc_diagnos(i_divu)=.true.
      endif
      if (idiag_ssruzm/=0 .or. idiag_ssuzm/=0 .or. idiag_ssm/=0 .or. &
          idiag_ss2m/=0 .or. idiag_ssmz/=0 .or. idiag_ss2mz/=0 .or. & 
          idiag_ssupmz/=0 .or. idiag_ssdownmz/=0 .or. &
          idiag_ss2upmz/=0 .or. idiag_ss2downmz/=0 .or. & 
          idiag_ssf2mz/=0 .or. idiag_ssf2upmz/=0 .or. idiag_ssf2downmz/=0 .or. & 
          idiag_ssmy/=0 .or. idiag_ssmx/=0 .or. idiag_ss2mx/=0 .or. & 
          idiag_ssmr/=0) & 
          lpenc_diagnos(i_ss)=.true.
      if (idiag_ppmx/=0 .or. idiag_ppmy/=0 .or. idiag_ppmz/=0) &
         lpenc_diagnos(i_pp)=.true.
      lpenc_diagnos(i_rho)=.true.
      lpenc_diagnos(i_ee)=.true.
      if (idiag_ethm/=0 .or. idiag_ethtot/=0 .or. idiag_ethdivum/=0 .or. &
          idiag_ethmz/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_ee)=.true.
      endif
      if (idiag_ssuzm/=0) lpenc_diagnos(i_uu)=.true.
      if (idiag_ssruzm/=0 .or. idiag_fconvm/=0 .or. idiag_fconvz/=0 .or. &
          idiag_Fenthz/=0 .or. idiag_Fenthupz/=0 .or. idiag_Fenthdownz/=0 .or. &
          idiag_fconvxmx/=0) then
        lpenc_diagnos(i_cp)=.true.
        lpenc_diagnos(i_uu)=.true.
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_TT)=.true.  !(to be replaced by enthalpy)
      endif
      if (idiag_fradz/=0) lpenc_diagnos(i_gTT)=.true.
      if (idiag_fconvxy/=0 .or. idiag_fconvyxy/=0 .or. idiag_fconvzxy/=0) then
        lpenc_diagnos2d(i_cp)=.true.
        lpenc_diagnos2d(i_uu)=.true.
        lpenc_diagnos2d(i_rho)=.true.
        lpenc_diagnos2d(i_TT)=.true.
      endif
      if (idiag_fturbz/=0 .or. idiag_fturbtz/=0 .or. idiag_fturbmz/=0 .or. &
          idiag_fturbfz/=0 .or. idiag_fturbmx/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_TT)=.true.
        lpenc_diagnos(i_gss)=.true.
      endif
      if (idiag_fturbxy/=0 .or. idiag_fturbrxy/=0 .or. &
          idiag_fturbthxy/=0 .or. idiag_fturbymxy/=0) then
        lpenc_diagnos2d(i_rho)=.true.
        lpenc_diagnos2d(i_TT)=.true.
        lpenc_diagnos2d(i_gss)=.true.
      endif
      if (idiag_fradxy_Kprof/=0 .or. idiag_fradymxy_Kprof/=0) then
        lpenc_diagnos2d(i_TT)=.true.
        lpenc_diagnos2d(i_glnTT)=.true.
      endif
      if (idiag_fradxy_kramers/=0) then
        lpenc_diagnos2d(i_rho)=.true.
        lpenc_diagnos2d(i_TT)=.true.
        lpenc_diagnos2d(i_glnTT)=.true.
      endif
      if (idiag_TTm/=0 .or. idiag_TTmx/=0 .or. idiag_TTmy/=0 .or. &
          idiag_TTmz/=0 .or. idiag_TTmr/=0 .or. idiag_TTmax/=0 .or. &
          idiag_TTdownmz/=0 .or. idiag_TTupmz/=0 .or. &
          idiag_TTmin/=0 .or. idiag_uxTTmz/=0 .or.idiag_uyTTmz/=0 .or. &
          idiag_uzTTmz/=0 .or. idiag_TT2mx/=0 .or. idiag_TT2mz/=0 .or. &
          idiag_uzTTupmz/=0 .or. idiag_uzTTdownmz/=0 .or. &
          idiag_TT2downmz/=0 .or. idiag_TT2upmz/=0 .or. &
          idiag_uxTTmx/=0 .or. idiag_uyTTmx/=0 .or. idiag_uzTTmx/=0) &
          lpenc_diagnos(i_TT)=.true.
      if (idiag_TTupmz/=0 .or. idiag_TTdownmz/=0 .or. &
          idiag_TT2upmz/=0 .or. idiag_TT2downmz/=0 .or. &
          idiag_TTf2mz/=0 .or. idiag_TTf2upmz/=0 .or. idiag_TTf2downmz/=0 .or. &
          idiag_ssupmz/=0 .or. idiag_ssdownmz/=0 .or. &
          idiag_ss2upmz/=0 .or. idiag_ss2downmz/=0 .or. &
          idiag_ssf2mz/=0 .or. idiag_ssf2upmz/=0 .or. idiag_ssf2downmz/=0 .or. &
          idiag_uzTTupmz/=0 .or. idiag_uzTTdownmz/=0) &
        lpenc_diagnos(i_uu)=.true.
      if (idiag_gTmax/=0) then
        lpenc_diagnos(i_glnTT) =.true.
        lpenc_diagnos(i_TT) =.true.
      endif
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
      if (idiag_fradbot/=0 .or. idiag_fradtop/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_cp)=.true.
        lpenc_diagnos(i_TT)=.true.
        lpenc_diagnos(i_glnTT)=.true.
      endif
      if (idiag_cs2mphi/=0) lpenc_diagnos(i_cs2)=.true.
!
!  diagnostics for baroclinic term
!
      if (idiag_gTrms/=0.or.idiag_gTxgsrms/=0) &
        lpenc_diagnos(i_gTT)=.true.
      if (idiag_gTxmxy/=0 .or. idiag_gTymxy/=0 .or. idiag_gTzmxy/=0 .or. &
          idiag_gTxgsxmz/=0 .or. idiag_gTxgsymz/=0 .or. idiag_gTxgszmz/=0 .or. &
          idiag_gTxgsx2mz/=0 .or. idiag_gTxgsy2mz/=0 .or. idiag_gTxgsz2mz/=0 .or. &
          idiag_gTxgsxmxy/=0 .or. idiag_gTxgsymxy/=0 .or. idiag_gTxgszmxy/=0 .or. &
          idiag_gTxgsx2mxy/=0 .or. idiag_gTxgsy2mxy/=0 .or. idiag_gTxgsz2mxy/=0) &
        lpenc_diagnos2d(i_gTT)=.true.
      if (idiag_gsrms/=0.or.idiag_gTxgsrms/=0) lpenc_diagnos(i_gss)=.true.
      if (idiag_gsxmxy/=0 .or. idiag_gsymxy/=0 .or. idiag_gszmxy/=0 .or. &
          idiag_gTxgsxmz/=0 .or. idiag_gTxgsymz/=0 .or. idiag_gTxgszmz/=0 .or. &
          idiag_gTxgsx2mz/=0 .or. idiag_gTxgsy2mz/=0 .or. idiag_gTxgsz2mz/=0 .or. &
          idiag_gTxgsxmxy/=0 .or. idiag_gTxgsymxy/=0 .or. idiag_gTxgszmxy/=0 .or. &
          idiag_gTxgsx2mxy/=0 .or. idiag_gTxgsy2mxy/=0 .or. idiag_gTxgsz2mxy/=0) &
        lpenc_diagnos2d(i_gss)=.true.
!
!  Cooling for cold core collapse
!
      if (lprestellar_cool_iso) then
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_divu)=.true.
        lpenc_requested(i_pp)=.true.
      endif
!
! Store initial stratification as pencils if f-array space allocated
!
      if (iglobal_lnrho0/=0) lpenc_requested(i_initlnrho)=.true.
      if (iglobal_ss0/=0) lpenc_requested(i_initss)=.true.
!
    endsubroutine pencil_criteria_energy
!***********************************************************************
    subroutine pencil_interdep_energy(lpencil_in)
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
      if (lpencil_in(i_sglnTT)) then
        lpencil_in(i_sij)=.true.
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
      if (lpencil_in(i_uuadvec_gss)) then
        lpencil_in(i_uu_advec)=.true.
        lpencil_in(i_gss)=.true.
      endif
!
    endsubroutine pencil_interdep_energy
!***********************************************************************
    subroutine calc_pencils_energy(f,p)
!
!  Calculate Entropy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-nov-04/anders: coded
!  15-mar-15/MR: changes for use of reference state.
!
      use EquationOfState, only: gamma1
      use Sub, only: u_dot_grad, grad, multmv, h_dot_grad
      use WENO_transport, only: weno_transp
!
      real, dimension(mx,my,mz,mfarray), intent(IN) :: f
      type(pencil_case),                 intent(OUT):: p
!
      integer :: j
!
! Ma2
      if (lpencil(i_Ma2)) p%Ma2=p%u2/p%cs2
! ugss
      if (lpencil(i_ugss)) then
        call u_dot_grad(f,iss,p%gss,p%uu,p%ugss,UPWIND=lupw_ss)
        if (lreference_state) p%ugss = p%ugss + p%uu(:,1)*reference_state(:,iref_gs)  ! suppress if grad(s)=0
      endif
! for pretend_lnTT
      if (lpencil(i_uglnTT)) &
          call u_dot_grad(f,iss,p%glnTT,p%uu,p%uglnTT,UPWIND=lupw_ss)
! sglnTT
      if (lpencil(i_sglnTT)) call multmv(p%sij,p%glnTT,p%sglnTT)
! dsdr
     !if (lpencil(i_dsdr)) then
     !    call grad(f,iss,gradS)
     !    p%dsdr=gradS(:,1)
     !endif
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
              p%fpres(:,j)=-p%cs2*(p%glnrho(:,j) + p%glnTT(:,j))*gamma1
            enddo
!  TH: The following would work if one uncomments the intrinsic operator
!  extensions in sub.f90. Please Test.
!          p%fpres      =-p%cs2*(p%glnrho + p%glnTT)*gamma1
!if(iproc==3) print*, 'n,m,fpres=', n,m,maxval(p%fpres),minval(p%fpres)
          else
            do j=1,3
              p%fpres(:,j)=-p%cs2*(p%glnrho(:,j) + p%cp1tilde*p%gss(:,j))
            enddo
          endif
        endif
      endif
!  transprhos
      if (lpencil(i_transprhos)) then
        if (lreference_state) then
          call weno_transp(f,m,n,iss,irho,iux,iuy,iuz,p%transprhos,dx_1,dy_1,dz_1, &
                           reference_state(:,iref_s), reference_state(:,iref_rho))
        else
          call weno_transp(f,m,n,iss,irho,iux,iuy,iuz,p%transprhos,dx_1,dy_1,dz_1)
        endif
      endif
! initlnrho
      if (lpencil(i_initlnrho).and.iglobal_lnrho0/=0) &
          p%initlnrho=f(l1:l2,m,n,iglobal_lnrho0)
! initss
      if (lpencil(i_initss).and.iglobal_ss0/=0) &
          p%initss=f(l1:l2,m,n,iglobal_ss0)
!
      if (lpencil(i_uuadvec_gss)) call h_dot_grad(p%uu_advec,p%gss,p%uuadvec_gss)
!
    endsubroutine calc_pencils_energy
!***********************************************************************
    subroutine denergy_dt(f,df,p)
!
!  Calculate right hand side of entropy equation,
!  ds/dt = -u.grads + [H-C + div(K*gradT) + mu0*eta*J^2 + ...]
!
!  17-sep-01/axel: coded
!   9-jun-02/axel: pressure gradient added to du/dt already here
!   2-feb-03/axel: added possibility of ionization
!  21-oct-15/MR: added timestep adaptation for slope-limited diffusion
! 
      use Diagnostics
      use EquationOfState, only: gamma1
      use Interstellar, only: calc_heat_cool_interstellar
      use Special, only: special_calc_energy
      use Sub
      use Viscosity, only: calc_viscous_heat
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(inout)  :: f,p
      intent(inout) :: df
!
      real :: ztop,xi,profile_cor
      real, dimension(nx) :: uduu
      integer :: j,i
!
      Hmax = 0.
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'denergy_dt: SOLVE denergy_dt'
      if (headtt) call identify_bcs('ss',iss)
      if (headtt) print*,'denergy_dt: lnTT,cs2,cp1=', p%lnTT(1), p%cs2(1), p%cp1(1)
!
!  ``cs2/dx^2'' for timestep
!
      if (ldensity) then
        if (lhydro.and.lfirst.and.ldt) then
          if (lreduced_sound_speed) then
            if (lscale_to_cs2top) then
              advec_cs2=reduce_cs2*cs2top*dxyz_2
            else
              advec_cs2=reduce_cs2*p%cs2*dxyz_2
            endif
          else
            advec_cs2=p%cs2*dxyz_2
          endif
        endif
        if (headtt.or.ldebug) &
          print*, 'denergy_dt: max(advec_cs2) =', maxval(advec_cs2)
      endif
!
!  Pressure term in momentum equation (setting lpressuregradient_gas to
!  .false. allows suppressing pressure term for test purposes).
!
      if (lhydro) then
        if (lpressuregradient_gas) then
          df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + p%fpres
!
!  If reference state is used, -grad(p')/rho is needed in momentum equation, hence fpres -> fpres + grad(p0)/rho.
!
          if (lreference_state) &
            df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + reference_state(:,iref_gp)*p%rho1
        endif
!
!  Add pressure force from global density gradient.
!
        if (maxval(abs(beta_glnrho_global))/=0.0) then
          if (headtt) print*, 'denergy_dt: adding global pressure gradient force'
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
!  NOTE: in the entropy module it is lnTT that is advanced, so
!  there are additional cv1 terms on the right hand side.
!
      if (ladvection_entropy) then
        if (pretend_lnTT) then
          df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) - p%divu*gamma_m1-p%uglnTT
        else
          if (lweno_transport) then
            df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
                - p%transprhos*p%rho1 + p%ss*p%rho1*p%transprho
          elseif (lfargo_advection) then
            df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) - p%uuadvec_gss
          else
            df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) - p%ugss
          endif
        endif
      endif
!
!  Add advection term from an imposed spatially constant gradient of S.
!  This makes sense really only for periodic boundary conditions.
!
        do j=1,3
          if (gradS0_imposed(j)/=0.) &
            df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)-gradS0_imposed(j)*p%uu(:,j)
        enddo
!
!  Calculate viscous contribution to entropy.
!
      diffus_chi=0.; diffus_chi3=0.
 
      if (lviscosity .and. lviscosity_heat) call calc_viscous_heat(df,p,Hmax)
!
!  Entry possibility for "personal" entries.
!  In that case you'd need to provide your own "special" routine.
!
      if (lspecial) call special_calc_energy(f,df,p)
!
!  Thermal conduction delegated to different subroutines.
!
      if (lheatc_Kprof)    call calc_heatcond(f,df,p)
      if (lheatc_Kconst)   call calc_heatcond_constK(df,p)
      if (lheatc_sfluct)   call calc_heatcond_sfluct(df,p)
      if (lheatc_chiconst) call calc_heatcond_constchi(f,df,p)
      if (lheatc_chi_cspeed) call calc_heatcond_cspeed_chi(df,p)
      if (lheatc_sqrtrhochiconst) call calc_heatcond_sqrtrhochi(df,p)
      if (lheatc_shock.or.lheatc_shock2)    call calc_heatcond_shock(df,p)
      if (lheatc_shock_profr)    call calc_heatcond_shock_profr(df,p)
      if (lheatc_hyper3ss) call calc_heatcond_hyper3(df,p)
      if (lheatc_spitzer)  call calc_heatcond_spitzer(df,p)
      if (lheatc_hubeny)   call calc_heatcond_hubeny(df,p)
      if (lheatc_kramers)  call calc_heatcond_kramers(f,df,p)
      if (lheatc_chit)     call calc_heatcond_chit(f,df,p)
      if (lheatc_smagorinsky)  call calc_heatcond_smagorinsky(f,df,p)
      if (lheatc_corona) then
        call calc_heatcond_spitzer(df,p)
        call newton_cool(df,p)
        call calc_heat_cool_RTV(df,p)
      endif
      if (lheatc_tensordiffusion) call calc_heatcond_tensor(df,p)
      if (lheatc_hyper3ss_polar) call calc_heatcond_hyper3_polar(f,df)
      if (lheatc_hyper3ss_mesh)  call calc_heatcond_hyper3_mesh(f,df)
      if (lheatc_hyper3ss_aniso) call calc_heatcond_hyper3_aniso(f,df)

      if (lfirst.and.ldt) then
        maxdiffus=max(maxdiffus,diffus_chi)
        maxdiffus3=max(maxdiffus3,diffus_chi3)
      endif
      !!!if (lenergy_slope_limited.and.lfirst) &
      !!!  df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)-f(l1:l2,m,n,iFF_div_ss)
!
      !if(lfirst .and. ldiagnos) print*,'DIV:iproc,m,f=',iproc,m,f(l1:l2,m,n,iFF_div_ss)
!
!  Explicit heating/cooling terms.
!
      if ((luminosity/=0.0) .or. (cool/=0.0) .or. &
          (tau_cor/=0.0) .or. (tauheat_buffer/=0.0) .or. &
          heat_uniform/=0.0 .or. heat_gaussianz/=0.0 .or. tau_cool/=0.0 .or. &
          tau_cool_ss/=0.0 .or. cool_uniform/=0.0 .or. cool_newton/=0.0 .or. &
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
!  Enforce maximum heating rate timestep constraint
!
      if (lfirst.and.ldt) then
        if (lthdiff_Hmax) then 
          ss0 = abs(df(l1:l2,m,n,iss))
          dt1_max=max(dt1_max,ss0*p%cv1/cdts)
        else
          dt1_max=max(dt1_max,abs(Hmax/p%ee/cdts))
        endif
      endif
!
!  Calculate entropy related diagnostics.
!
      if (ldiagnos) then
        if (idiag_uduum/=0) then
            uduu=0.
            do i = 1, 3
              uduu=uduu+p%uu(:,i)*df(l1:l2,m,n,iux-1+i)
            enddo
            call sum_mn_name(p%rho*uduu,idiag_uduum)
        endif
      endif
!
      call calc_diagnostics_energy(f,p)

    endsubroutine denergy_dt
!***********************************************************************
    subroutine calc_0d_diagnostics_energy(p)
!
!  Calculate entropy related diagnostics.
!
      use Diagnostics
      use Sub, only: cross, dot2

      type(pencil_case) :: p

      real, dimension(nx) :: ufpres, glnTT2, Ktmp
      real, dimension(nx) :: gT2,gs2,gTxgs2
      real, dimension(nx,3) :: gTxgs
      real :: uT,fradz,TTtop
      integer :: i

      if (ldiagnos) then

        !uT=unit_temperature !(define shorthand to avoid long lines below)
        uT=1. !(AB: for the time being; to keep compatible with auto-test
        
        if (.not.lgpu) then
          if (idiag_dtchi/=0) &
              call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
          if (idiag_dtH/=0) then
            if (lthdiff_Hmax) then
              call max_mn_name(ss0*p%cv1/cdts,idiag_dtH,l_dt=.true.)
            else
              call max_mn_name(Hmax/p%ee/cdts,idiag_dtH,l_dt=.true.)
            endif
          endif
          if (idiag_Hmax/=0) call max_mn_name(Hmax/p%ee,idiag_Hmax)
          if (idiag_tauhmin/=0) then
            if (lthdiff_Hmax) then
              call max_mn_name(ss0*p%cv1,idiag_tauhmin,lreciprocal=.true.)
            else
              call max_mn_name(Hmax/p%ee,idiag_tauhmin,lreciprocal=.true.)
            endif
          endif
          if (idiag_dtc/=0) &
            call max_mn_name(sqrt(advec_cs2)/cdt,idiag_dtc,l_dt=.true.)
        endif
        if (idiag_ssmax/=0) call max_mn_name(p%ss*uT,idiag_ssmax)
        if (idiag_ssmin/=0) call max_mn_name(-p%ss*uT,idiag_ssmin,lneg=.true.)
        if (idiag_TTmax/=0) call max_mn_name(p%TT*uT,idiag_TTmax)
        if (idiag_TTmin/=0) call max_mn_name(-p%TT*uT,idiag_TTmin,lneg=.true.)
        if (idiag_gTmax/=0) then
          call dot2(p%glnTT,glnTT2)
          call max_mn_name(p%TT*sqrt(glnTT2),idiag_gTmax)
        endif
        if (idiag_TTm/=0)    call sum_mn_name(p%TT*uT,idiag_TTm)
        if (idiag_pdivum/=0) call sum_mn_name(p%pp*p%divu,idiag_pdivum)
        call max_mn_name(p%yH,idiag_yHmax)
        call sum_mn_name(p%yH,idiag_yHm)
        if (idiag_ethm/=0)   call sum_mn_name(p%rho*p%ee,idiag_ethm)
        if (idiag_ethtot/=0) call integrate_mn_name(p%rho*p%ee,idiag_ethtot)
        if (idiag_ethdivum/=0) &
            call sum_mn_name(p%rho*p%ee*p%divu,idiag_ethdivum)
        if (idiag_ssruzm/=0) call sum_mn_name(p%ss*p%rho*p%uu(:,3),idiag_ssruzm)
        if (idiag_ssuzm/=0) call sum_mn_name(p%ss*p%uu(:,3),idiag_ssuzm)
        call sum_mn_name(p%ss,idiag_ssm)
        if (idiag_ss2m/=0) call sum_mn_name(p%ss**2,idiag_ss2m)
        call sum_mn_name(p%ee,idiag_eem)
        call sum_mn_name(p%pp,idiag_ppm)
        call sum_mn_name(p%cs2,idiag_csm,lsqrt=.true.)
        call max_mn_name(p%cs2,idiag_csmax,lsqrt=.true.)
        if (idiag_cgam/=0) call sum_mn_name(16.*real(sigmaSB)*p%TT**3*p%cp1*p%rho1,idiag_cgam)
        if (idiag_ugradpm/=0) &
            call sum_mn_name(p%cs2*(p%uglnrho+p%ugss),idiag_ugradpm)
        if (idiag_fconvm/=0) &
            call sum_mn_name(p%cp*p%rho*p%uu(:,3)*p%TT,idiag_fconvm)
        if (idiag_ufpresm/=0) then
            ufpres=0.
            do i = 1,3
              ufpres=ufpres+p%uu(:,i)*p%fpres(:,i)
            enddo
            call sum_mn_name(ufpres,idiag_ufpresm)
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
        if (idiag_fradbot/=0.and.lfirst_proc_z.and.n==n1) then
          if (hcond0==0.) then
            Ktmp=chi*p%rho*p%cp
          else
            Ktmp=hcond0
          endif
          fradz=sum(-Ktmp*p%TT*p%glnTT(:,3)*dsurfxy)
          call surf_mn_name(fradz,idiag_fradbot)
        endif
!
!  Radiative heat flux at the top (assume here that hcond=hcond0=const).
!
        if (idiag_fradtop/=0.and.llast_proc_z.and.n==n2) then
          if (hcond0==0.) then
            Ktmp=chi*p%rho*p%cp
          else
            Ktmp=hcond0
          endif
          fradz=sum(-Ktmp*p%TT*p%glnTT(:,3)*dsurfxy)
          call surf_mn_name(fradz,idiag_fradtop)
        endif
!
!  Mean temperature at the top.
!
        if (idiag_TTtop/=0.and.llast_proc_z.and.n==n2) then
          TTtop=sum(p%TT*dsurfxy)
          call surf_mn_name(TTtop,idiag_TTtop)
        endif
!
!  Calculate integrated temperature in in limited radial range.
!
        if (idiag_TTp/=0) call sum_lim_mn_name(p%rho*p%cs2*gamma1,idiag_TTp,p)

      endif

    endsubroutine calc_0d_diagnostics_energy
!***********************************************************************
    subroutine calc_1d_diagnostics_energy(f,p)
!
      use Diagnostics
      use Sub, only: cross
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
      type(pencil_case) :: p
!
      real, dimension (nx,3) :: gTxgs
      real, dimension (nx) :: uzmask
!
!  1-D averages.
!
      if (l1davgfirst) then
!
        if (idiag_ethmz/=0) call xysum_mn_name_z(p%rho*p%ee,idiag_ethmz)
        if (idiag_fradz/=0) call xysum_mn_name_z(-hcond0*p%gTT(:,3),idiag_fradz)
        if (idiag_fconvz/=0) call xysum_mn_name_z(p%cp*p%rho*p%uu(:,3)*p%TT,idiag_fconvz)
        if (idiag_fconvxmx/=0) call yzsum_mn_name_x(p%cp*p%rho*p%uu(:,1)*p%TT,idiag_fconvxmx)
        call yzsum_mn_name_x(p%ss,idiag_ssmx)
        if (idiag_ss2mx/=0) call yzsum_mn_name_x(p%ss**2,idiag_ss2mx)
        call xzsum_mn_name_y(p%ss,idiag_ssmy)
        call xysum_mn_name_z(p%ss,idiag_ssmz)
        if (idiag_ss2mz/=0) call xysum_mn_name_z(p%ss**2,idiag_ss2mz)
        call yzsum_mn_name_x(p%pp,idiag_ppmx)
        call xzsum_mn_name_y(p%pp,idiag_ppmy)
        call xysum_mn_name_z(p%pp,idiag_ppmz)
        call yzsum_mn_name_x(p%TT,idiag_TTmx)
        call xzsum_mn_name_y(p%TT,idiag_TTmy)
        call xysum_mn_name_z(p%TT,idiag_TTmz)
        if (idiag_fradz_constchi/=0) call xysum_mn_name_z(-chi*p%rho*p%TT*p%glnTT(:,3)/p%cp1,idiag_fradz_constchi)
        if (idiag_TT2mx/=0) call yzsum_mn_name_x(p%TT**2,idiag_TT2mx)
        if (idiag_TT2mz/=0) call xysum_mn_name_z(p%TT**2,idiag_TT2mz)
        if (idiag_uxTTmz/=0) call xysum_mn_name_z(p%uu(:,1)*p%TT,idiag_uxTTmz)
        if (idiag_uxTTmx/=0) call yzsum_mn_name_x(p%uu(:,1)*p%TT,idiag_uxTTmx)
        if (idiag_uyTTmx/=0) call yzsum_mn_name_x(p%uu(:,2)*p%TT,idiag_uyTTmx)
        if (idiag_uzTTmx/=0) call yzsum_mn_name_x(p%uu(:,3)*p%TT,idiag_uzTTmx)
        if (idiag_uyTTmz/=0) call xysum_mn_name_z(p%uu(:,2)*p%TT,idiag_uyTTmz)
        if (idiag_uzTTmz/=0) call xysum_mn_name_z(p%uu(:,3)*p%TT,idiag_uzTTmz)
        call phizsum_mn_name_r(p%ss,idiag_ssmr)
        call phizsum_mn_name_r(p%TT,idiag_TTmr)
        if (lFenth_as_aux) &
            call xysum_mn_name_z(f(l1:l2,m,n,iFenth),idiag_Fenthz)
        if (lss_flucz_as_aux) then
            if (idiag_ssf2mz/=0) call xysum_mn_name_z(f(l1:l2,m,n,iss_flucz)**2,idiag_ssf2mz)
        endif
        if (lTT_flucz_as_aux) then
          if (idiag_TTf2mz/=0) call xysum_mn_name_z(f(l1:l2,m,n,iTT_flucz)**2,idiag_TTf2mz)
        endif
!
        if (idiag_Fenthupz/=0 .or. idiag_ssupmz/=0 .or. idiag_ss2upmz/=0 .or. &
            idiag_ssf2upmz/=0 .or. idiag_TTupmz/=0 .or. idiag_TT2upmz/=0 .or. &
            idiag_TTf2upmz/=0 .or. idiag_uzTTupmz/=0) then
          where (p%uu(:,3) > 0.)
            uzmask = p%uu(:,3)/abs(p%uu(:,3))
          elsewhere
            uzmask = 0.
          endwhere
!
          if (lFenth_as_aux) then
            if (idiag_Fenthupz/=0) call xysum_mn_name_z(uzmask*f(l1:l2,m,n,iFenth),idiag_Fenthupz)
          endif
          if (idiag_ssupmz  /=0) call xysum_mn_name_z(uzmask*p%ss,idiag_ssupmz)
          if (idiag_ss2upmz /=0) call xysum_mn_name_z(uzmask*p%ss**2,idiag_ss2upmz)
          if (lss_flucz_as_aux) then
            if (idiag_ssf2upmz/=0) call xysum_mn_name_z(uzmask*f(l1:l2,m,n,iss_flucz)**2,idiag_ssf2upmz)
          endif
          if (idiag_TTupmz  /=0) call xysum_mn_name_z(uzmask*p%TT,idiag_TTupmz)
          if (idiag_TT2upmz /=0) call xysum_mn_name_z(uzmask*p%TT**2,idiag_TT2upmz)
          if (lTT_flucz_as_aux) then
            if (idiag_TTf2upmz/=0) call xysum_mn_name_z(uzmask*f(l1:l2,m,n,iTT_flucz)**2,idiag_TTf2upmz)
          endif
          if (idiag_uzTTupmz/=0) call xysum_mn_name_z(uzmask*p%uu(:,3)*p%TT,idiag_uzTTupmz)
        endif
!
        if (idiag_Fenthdownz/=0 .or. idiag_ssdownmz/=0 .or. idiag_ss2downmz/=0 .or. &
            idiag_ssf2downmz/=0 .or. idiag_TTdownmz/=0 .or. idiag_TT2downmz/=0 .or. &
            idiag_TTf2downmz/=0 .or. idiag_uzTTdownmz/=0) then
          where (p%uu(:,3) < 0.)
            uzmask = -p%uu(:,3)/abs(p%uu(:,3))
          elsewhere
            uzmask = 0.
          endwhere
          if (lFenth_as_aux) then
            if (idiag_Fenthdownz/=0) call xysum_mn_name_z(uzmask*f(l1:l2,m,n,iFenth),idiag_Fenthdownz)
          endif
          if (idiag_ssdownmz  /=0) call xysum_mn_name_z(uzmask*p%ss,idiag_ssdownmz)
          if (idiag_ss2downmz /=0) call xysum_mn_name_z(uzmask*p%ss**2,idiag_ss2downmz)
          if (lss_flucz_as_aux) then
            if (idiag_ssf2downmz/=0) call xysum_mn_name_z(uzmask*f(l1:l2,m,n,iss_flucz)**2,idiag_ssf2downmz)
          endif
          if (idiag_TTdownmz  /=0) call xysum_mn_name_z(uzmask*p%TT,idiag_TTdownmz)
          if (idiag_TT2downmz /=0) call xysum_mn_name_z(uzmask*p%TT**2,idiag_TT2downmz)
          if (lTT_flucz_as_aux) then
            if (idiag_TTf2downmz/=0) call xysum_mn_name_z(uzmask*f(l1:l2,m,n,iTT_flucz)**2,idiag_TTf2downmz)
          endif
          if (idiag_uzTTdownmz/=0) call xysum_mn_name_z(uzmask*p%uu(:,3)*p%TT,idiag_uzTTdownmz)
        endif
!
!  For the 1D averages of the baroclinic term
!
        if (idiag_gTxgsxmz/=0 .or. idiag_gTxgsx2mz/=0 .or. &
            idiag_gTxgsymz/=0 .or. idiag_gTxgsy2mz/=0 .or. &
            idiag_gTxgszmz/=0 .or. idiag_gTxgsz2mz/=0) then
          call cross(p%gTT,p%gss,gTxgs)
          call xysum_mn_name_z(gTxgs(:,1),idiag_gTxgsxmz)
          call xysum_mn_name_z(gTxgs(:,2),idiag_gTxgsymz)
          call xysum_mn_name_z(gTxgs(:,3),idiag_gTxgszmz)
          if (idiag_gTxgsx2mz/=0) &
              call xysum_mn_name_z(gTxgs(:,1)**2,idiag_gTxgsx2mz)
          if (idiag_gTxgsy2mz/=0) &
              call xysum_mn_name_z(gTxgs(:,2)**2,idiag_gTxgsy2mz)
          if (idiag_gTxgsz2mz/=0) &
              call xysum_mn_name_z(gTxgs(:,3)**2,idiag_gTxgsz2mz)
        endif
      endif
!
    endsubroutine calc_1d_diagnostics_energy
!***********************************************************************
    subroutine calc_2d_diagnostics_energy(p)
!
!  2-D averages.
!
      use Diagnostics
      use Sub, only: cross

      type(pencil_case) :: p

      real, dimension (nx,3) :: gTxgs, tmpvec

      if (l2davgfirst) then
!
!  Phi-averages
!
        call phisum_mn_name_rz(p%ss,idiag_ssmphi)
        call phisum_mn_name_rz(p%cs2,idiag_cs2mphi)

        call zsum_mn_name_xy(p%TT,idiag_TTmxy)
        call ysum_mn_name_xz(p%TT,idiag_TTmxz)
        call zsum_mn_name_xy(p%ss,idiag_ssmxy)
        call ysum_mn_name_xz(p%ss,idiag_ssmxz)
        if (idiag_uxTTmxy/=0) call zsum_mn_name_xy(p%uu(:,1)*p%TT,idiag_uxTTmxy)
        call zsum_mn_name_xy(p%uu,idiag_uyTTmxy,(/0,1,0/),p%TT)
        call zsum_mn_name_xy(p%uu,idiag_uzTTmxy,(/0,0,1/),p%TT)
        if (idiag_fconvxy/=0) &
            call zsum_mn_name_xy(p%cp*p%rho*p%uu(:,1)*p%TT,idiag_fconvxy)
        if (idiag_fconvyxy/=0) &
            call zsum_mn_name_xy(p%uu,idiag_fconvyxy,(/0,1,0/),p%cp*p%rho*p%TT)
        if (idiag_fconvzxy/=0) &
            call zsum_mn_name_xy(p%uu,idiag_fconvzxy,(/0,0,1/),p%cp*p%rho*p%TT)
        call zsum_mn_name_xy(p%gTT(:,1),idiag_gTxmxy)
        call zsum_mn_name_xy(p%gTT,idiag_gTymxy,(/0,1,0/))
        call zsum_mn_name_xy(p%gTT,idiag_gTzmxy,(/0,0,1/))
        call zsum_mn_name_xy(p%gss(:,1),idiag_gsxmxy)
        call zsum_mn_name_xy(p%gss,idiag_gsymxy,(/0,1,0/))
        call zsum_mn_name_xy(p%gss,idiag_gszmxy,(/0,0,1/))

        if (chi_t/=0..and.chit_aniso/=0.) then
          if (idiag_fturbrxy/=0) &
            call zsum_mn_name_xy(-chit_aniso_prof*p%rho*p%TT* &
                                 (costh(m)**2*p%gss(:,1)-sinth(m)*costh(m)*p%gss(:,2)),idiag_fturbrxy)
          if (idiag_fturbthxy/=0) then
            tmpvec(:,(/1,3/))=0.
            tmpvec(:,2)= -chit_aniso_prof*p%rho*p%TT* &
                         (-sinth(m)*costh(m)*p%gss(:,1)+sinth(m)**2*p%gss(:,2))
            call zsum_mn_name_xy(tmpvec,idiag_fturbthxy,(/0,1,0/))   ! not correct for Yin-Yang: phi component in tmpvec missing

          endif
        endif
!
!  2D averages of the baroclinic term.
!
         if (idiag_gTxgsxmxy/=0 .or. idiag_gTxgsx2mxy/=0 .or. &
             idiag_gTxgsymxy/=0 .or. idiag_gTxgsy2mxy/=0 .or. &
             idiag_gTxgszmxy/=0 .or. idiag_gTxgsz2mxy/=0) then
           call cross(p%gTT,p%gss,gTxgs)
           call zsum_mn_name_xy(gTxgs(:,1),idiag_gTxgsxmxy)
           call zsum_mn_name_xy(gTxgs,idiag_gTxgsymxy,(/0,1,0/))
           call zsum_mn_name_xy(gTxgs,idiag_gTxgszmxy,(/0,0,1/))
           if (idiag_gTxgsx2mxy/=0) &
               call zsum_mn_name_xy(gTxgs(:,1)**2,idiag_gTxgsxmxy)
           if (idiag_gTxgsy2mxy/=0) &
               call zsum_mn_name_xy(gTxgs**2,idiag_gTxgsymxy,(/0,1,0/))
           if (idiag_gTxgsz2mxy/=0) &
               call zsum_mn_name_xy(gTxgs**2,idiag_gTxgszmxy,(/0,0,1/))
         endif
       endif

    endsubroutine calc_2d_diagnostics_energy
!***********************************************************************
    subroutine calc_diagnostics_energy(f,p)

      real, dimension (mx,my,mz,mfarray) :: f
      type(pencil_case) :: p

      call calc_2d_diagnostics_energy(p)
      call calc_1d_diagnostics_energy(f,p)
      call calc_0d_diagnostics_energy(p)

    endsubroutine calc_diagnostics_energy
!***********************************************************************
    subroutine calc_ssmeanxy(f,ssmxy)
!
!  Calculate azimuthal (z) average of entropy
!
!  26-feb-14/MR: outsourced from energy_after_boundary
!
      use Sub, only: finalize_aver
!
      real, dimension (mx,my,mz,mfarray), intent(IN) :: f
      real, dimension (mx,my),            intent(OUT):: ssmxy
!
      integer :: l,m
      real :: fact
!
      fact=1./nzgrid
      do l=1,mx
        do m=1,my
          ssmxy(l,m)=fact*sum(f(l,m,n1:n2,iss))
        enddo
      enddo
!
      call finalize_aver(nprocz,3,ssmxy)
!
    endsubroutine calc_ssmeanxy
!***********************************************************************
    subroutine energy_after_boundary(f)
!
!  Calculate <s>, which is needed for diffusion with respect to xy-flucts.
!
!  17-apr-10/axel: adapted from magnetic_after_boundary
!  12-feb-15/MR  : changed for reference state; not yet done in averages of entropy.
!
      use Deriv, only: der_x, der2_x, der_z, der2_z
      use Mpicomm, only: mpiallreduce_sum, stop_it
      use Sub, only: finalize_aver,calc_all_diff_fluxes,div
      use EquationOfState, only : lnrho0, cs20, get_cv1, get_cp1
!
      real, dimension (mx,my,mz,mfarray),intent(INOUT) :: f
!
      real, dimension (mx,my,mz) :: cs2p, ruzp
      real, dimension (mz) :: ruzmz
!
      integer :: l,m,n,lf
      real :: fact, cv1, cp1, tmp1
!
!  Compute horizontal average of entropy. Include the ghost zones,
!  because they have just been set.
!
      if (lcalc_ssmean) then
        fact=1./nxygrid
        do n=1,mz
          ssmz(n)=fact*sum(f(l1:l2,m1:m2,n,iss))
        enddo
        call finalize_aver(nprocxy,12,ssmz)
!
!  Entropy fluctuations as auxilliary array
!
        if (lss_flucz_as_aux) then
        do n=1,mz
           f(l1:l2,m1:m2,n,iss_flucz) = f(l1:l2,m1:m2,n,iss) - ssmz(n)
        enddo
       endif
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
        fact=1./nyzgrid
        do l=1,mx
          ssmx(l)=fact*sum(f(l,m1:m2,n1:n2,iss))
        enddo
        call finalize_aver(nprocyz,23,ssmx)
!
        gssmx(:,2:3)=0.
        call der_x(ssmx,gssmx(:,1))
        call der2_x(ssmx,del2ssmx)
        if (lspherical_coords) del2ssmx=del2ssmx+2.*gssmx(:,1)/x(l1:l2)
!
!  Azimuthal (z) average
!
        fact=1./nzgrid
        do l=1,nx
          lf=l+nghost
          do m=1,my
            ssmxy(l,m)=fact*sum(f(lf,m,n1:n2,iss))
          enddo
        enddo
        call finalize_aver(nprocz,3,ssmxy)
!
      endif
!
!  Compute average sound speed cs2(x) and cs2(x,y)
!
      if (lcalc_cs2mean) then
        call get_cv1(cv1)
        if (lcartesian_coords) then
          fact=1./nyzgrid
          do l=1,nx
            lf=l+nghost
            if (ldensity_nolog) then
              if (lreference_state) then
                cs2mx(l)= fact*sum(cs20*exp(gamma_m1*(alog(f(lf,m1:m2,n1:n2,irho)+reference_state(l,iref_rho)) &
                         -lnrho0)+cv1*(f(lf,m1:m2,n1:n2,iss)+reference_state(l,iref_s)))) 
              else
                cs2mx(l)= fact*sum(cs20*exp(gamma_m1*(alog(f(lf,m1:m2,n1:n2,irho)) &
                         -lnrho0)+cv1*f(lf,m1:m2,n1:n2,iss)))
              endif
            else
              cs2mx(l)= fact*sum(cs20*exp(gamma_m1*(f(lf,m1:m2,n1:n2,ilnrho) &
                       -lnrho0)+cv1*f(lf,m1:m2,n1:n2,iss)))
            endif
          enddo
        elseif (lspherical_coords) then
          fact=1./((cos(y0)-cos(y0+Lxyz(2)))*Lxyz(3))
          do l=1,nx
            lf=l+nghost
            tmp1=0.
            do n=n1,n2
              if (ldensity_nolog) then
                if (lreference_state) then
                  tmp1= tmp1+sum(cs20*exp(gamma_m1*(alog(f(lf,m1:m2,n,irho)+reference_state(l,iref_rho)) &
                       -lnrho0)+cv1*(f(lf,m1:m2,n,iss)+reference_state(l,iref_s)))*dVol_y(m1:m2))*dVol_z(n)
                else
                  tmp1= tmp1+sum(cs20*exp(gamma_m1*(alog(f(lf,m1:m2,n,irho)) &
                       -lnrho0)+cv1*f(lf,m1:m2,n,iss))*dVol_y(m1:m2))*dVol_z(n)
                endif
              else
                tmp1= tmp1+sum(cs20*exp(gamma_m1*(f(lf,m1:m2,n,ilnrho) &
                     -lnrho0)+cv1*f(lf,m1:m2,n,iss))*dVol_y(m1:m2))*dVol_z(n)
              endif
            enddo
            cs2mx(l)=tmp1
          enddo
          cs2mx=fact*cs2mx
        elseif (lcylindrical_coords) then
          call fatal_error('energy_after_boundary','calculation of mean c_s not implemented for cylidrical coordinates')
        endif
        call finalize_aver(nprocyz,23,cs2mx)
!
!  Do 2D averages at the same time
!
        fact=1./nzgrid
        cs2mxy = 0.
!
        do l=1,nx
          lf=l+nghost
          do m=1,my
            if (ldensity_nolog) then
              if (lreference_state) then
                cs2mxy(l,m)= fact*sum(cs20*exp(gamma_m1*(alog(f(lf,m,n1:n2,irho)+reference_state(l,iref_rho)) &
                            -lnrho0)+cv1*(f(lf,m,n1:n2,iss)+reference_state(l,iref_s))))
              else
                cs2mxy(l,m)= fact*sum(cs20*exp(gamma_m1*(alog(f(lf,m,n1:n2,irho)) &
                            -lnrho0)+cv1*f(lf,m,n1:n2,iss)))
              endif
            else
              cs2mxy(l,m)= fact*sum(cs20*exp(gamma_m1*(f(lf,m,n1:n2,ilnrho) &
                          -lnrho0)+cv1*f(lf,m,n1:n2,iss)))
            endif
          enddo
        enddo
        call finalize_aver(nprocz,3,cs2mxy)
!
      endif
!
!  Compute average sound speed cs2(z)
!
      if (lcalc_cs2mz_mean .or. lcalc_cs2mz_mean_diag) then
        call get_cv1(cv1)
!
        fact=1./nxygrid
        cs2mz=0.
        if (ldensity_nolog) then
          if (lreference_state) then
            do n=1,mz
              cs2mz(n)= fact*sum(cs20*exp(gamma_m1*(alog(f(l1:l2,m1:m2,n,irho)+spread(reference_state(:,iref_rho),2,ny)) &
                       -lnrho0)+cv1*(f(l1:l2,m1:m2,n,iss)+reference_state(l-l1+1,iref_s))))
            enddo
          else
            do n=1,mz
              cs2mz(n)= fact*sum(cs20*exp(gamma_m1*(alog(f(l1:l2,m1:m2,n,irho)) &
                       -lnrho0)+cv1*f(l1:l2,m1:m2,n,iss)))
            enddo
          endif
        else
          do n=1,mz
            cs2mz(n)= fact*sum(cs20*exp(gamma_m1*(f(l1:l2,m1:m2,n,ilnrho) &
                     -lnrho0)+cv1*f(l1:l2,m1:m2,n,iss)))
          enddo
        endif
        call finalize_aver(nprocxy,12,cs2mz)
!
!  Sound speed fluctuations as auxilliary array
!
        if (lTT_flucz_as_aux) then
          call get_cp1(cp1)
          tmp1=cp1/gamma_m1
          do n=1,mz
            f(l1:l2,m1:m2,n,iTT_flucz) = tmp1* &
                         (cs20*exp(gamma_m1*(f(l1:l2,m1:m2,n,ilnrho)-lnrho0)+cv1*f(l1:l2,m1:m2,n,iss))-cs2mz(n))
          enddo
        endif
      endif
!
!  Compute volume average of entropy.
!
      if (lcalc_ss_volaverage) then
        tmp1=sum(f(l1:l2,m1:m2,n1:n2,iss))/nwgrid
        call mpiallreduce_sum(tmp1,ss_volaverage)
      endif
!
!  Slope limited diffusion following Rempel (2014)
!  First calculating the flux in a subroutine below
!  using a slope limiting procedure then storing in the
!  auxilaries variables in the f array. Finally the divergence
!  of the flux is calculated and stored in the f array.
!
      if (lenergy_slope_limited.and.lfirst) then
  
        f(:,:,:,iFF_diff1:iFF_diff2)=0.
        call calc_all_diff_fluxes(f,iss,islope_limiter,h_slope_limited)
!if (ldiagnos.and.iproc==0) print'(a,i2,22(1x,e14.8))','iss,IPROC=', iproc, f(4,10,:,iss)
!if (ldiagnos.and.iproc==0) print'(a,i2,22(1x,e14.8)/)','iFF_diff1,IPROC=', iproc, f(4,10,:,iFF_diff1)
        do n=n1,n2; do m=m1,m2
          call div(f,iFF_diff,f(l1:l2,m,n,iFF_div_ss),.true.)
        enddo; enddo

      endif
!
!  Compute running average of entropy
!
      if (lss_running_aver_as_aux) then
        if (t.lt.dt) f(:,:,:,iss_run_aver)=f(:,:,:,iss)
        f(:,:,:,iss_run_aver)=(1.-dt/tau_aver1)*f(:,:,:,iss_run_aver)+dt/tau_aver1*f(:,:,:,iss)
      endif
!
!  Enthalpy flux as auxilliary array
!
      if (lFenth_as_aux) then
!
        call get_cp1(cp1)
        tmp1=cp1/gamma_m1
!
        f(:,:,:,iFenth)=0.
!
        fact=1./nxygrid
        cs2p=0.
        ruzp=0.
        ruzmz=0.
        do n=1,mz
           ruzmz(n)= fact*sum(exp(f(l1:l2,m1:m2,n,ilnrho))*f(l1:l2,m1:m2,n,iuz))
        enddo
        call finalize_aver(nprocxy,12,ruzmz)
!
        do n=1,mz
           cs2p(l1:l2,m1:m2,n) = cs20*exp(gamma_m1*(f(l1:l2,m1:m2,n,ilnrho)-lnrho0)+cv1*f(l1:l2,m1:m2,n,iss))-cs2mz(n)
           ruzp(l1:l2,m1:m2,n) = exp(f(l1:l2,m1:m2,n,ilnrho))*f(l1:l2,m1:m2,n,iuz)-ruzmz(n)
        enddo
!
        f(l1:l2,m1:m2,:,iFenth)=tmp1*cs2p(l1:l2,m1:m2,:)*ruzp(l1:l2,m1:m2,:)
      endif
!
    endsubroutine energy_after_boundary
!***********************************************************************
    subroutine update_char_vel_energy(f)
!
!  Updates characteristic velocity for slope-limited diffusion.
!
!  25-sep-15/MR+joern: coded
!   9-oct-15/MR: added updating of characteristic velocity by sound speed
!  29-dec-15/joern: changed to staggered_max_scale
!
      use EquationOfState, only: eoscalc
!      use General, only: staggered_mean_scal
      use General, only: staggered_max_scal
!
      real, dimension(mx,my,mz,mfarray), intent(INOUT) :: f
!
      real, dimension(mx) :: cs2
!
      if (lslope_limit_diff) then
!
!  Calculate sound speed and store temporarily in first slot of diffusive fluxes.
!
        do n=1,mz; do m=1,my
          call eoscalc(f,mx,cs2=cs2)
          f(:,m,n,iFF_diff) = sqrt(cs2)   ! sqrt needed as we need the speed.
        enddo; enddo
!
!        call staggered_mean_scal(f,iFF_diff,iFF_char_c,w_sldchar_ene)
        call staggered_max_scal(f,iFF_diff,iFF_char_c,w_sldchar_ene)
!
      endif
!
    endsubroutine update_char_vel_energy
!***********************************************************************
    subroutine set_border_entropy(f,df,p)
!
!  Calculates the driving term for the border profile
!  of the ss variable.
!
!  28-jul-06/wlad: coded
!  28-apr-16/wlad: added case initial-temperature
!
      use BorderProfiles, only: border_driving, set_border_initcond
      use EquationOfState, only: gamma,gamma_m1,get_cp1
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(nx) :: f_target
      real, dimension(nx) :: lnrho_init,ss_init
      real :: cp1
!
      select case (borderss)
!
      case ('zero','0')
        f_target=0.0
      case ('constant')
        f_target=ss_const
      case ('initial-condition')
        call set_border_initcond(f,iss,f_target)
      case ('initial-temperature')
!
!  This boundary condition drives the entropy back not to the initial entropy,
!  but to the initial _temperature_. It is useful for situations when the density
!  in the boundary goes very low and the initial entropy corresponds to a very high         
!  temperature, making the simulation break. The sounds speed is 
!
!    cs2=cs20*exp(ss/cv+gamma_m1*(lnrho-lnrho0))
!
!  So the entropy it corresponds to is 
!
!    ss = cv*(log(cs2/cs20)-gamma_m1*(lnrho-lnrho0))
!         
!  Keep the density and replace cs2 by the initial cs2 to drive ss to the initial temperature. 
!
        call set_border_initcond(f,ilnrho,lnrho_init)
        call set_border_initcond(f,iss,ss_init)
        call get_cp1(cp1)
        !cs2_init = cs20*exp(gamma*cp1*ss_init+gamma_m1*(lnrho_init-lnrho0))
        !f_target = 1./(gamma*cp1)*(log(cs2_init/cs20)-gamma_m1*(p%lnrho-lnrho0))
!
!  The two lines above reduce to the one below        
!        
        f_target = ss_init - gamma_m1/(gamma*cp1)*(p%lnrho-lnrho_init)
        
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
    subroutine calc_heatcond_constchi(f,df,p)
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
      use Diagnostics!, only: max_mn_name
      use Sub, only: dot, multmv_transp, multsv_mn
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(inout) :: df
      intent(in) :: p
!
      real, dimension(nx) :: thdiff, g2, chit_prof
      real, dimension(nx,3) :: gradchit_prof
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
      call dot(p%glnrho+p%glnTT,p%glnTT,g2)
      if (pretend_lnTT) then
        thdiff=gamma*chi*(p%del2lnTT+g2)
      else
        thdiff=chi*p%cp*(p%del2lnTT+g2)
      endif
      if (chi_t/=0.) then
        call get_prof_pencil(chit_prof,gradchit_prof,lsphere_in_a_box,1.,chit_prof1,chit_prof2,xbot,xtop,p,f, &
                             stored_prof=chit_prof_stored,stored_dprof=dchit_prof_stored)

        call dot(p%glnrho+p%glnTT,p%gss,g2)
        thdiff=thdiff+chi_t*chit_prof*(p%del2ss+g2)
        call dot(gradchit_prof,p%gss,g2)
        thdiff=thdiff+chi_t*g2
!
!  Diagnostics
!
        if (l1davgfirst) then
          if (idiag_fturbz/=0) &
            call xysum_mn_name_z(-chi_t*chit_prof*p%rho*p%TT*p%gss(:,3),idiag_fturbz)
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
          diffus_chi=diffus_chi+(gamma*chi)*dxyz_2
        else
          diffus_chi=diffus_chi+chi*dxyz_2
        endif
        if (chi_t/=0.) diffus_chi = diffus_chi+chi_t*chit_prof*dxyz_2
!        if (ldiagnos.and.idiag_dtchi/=0) &
!            call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
      endif
!
    endsubroutine calc_heatcond_constchi
!***********************************************************************
    subroutine calc_heatcond_cspeed_chi(df,p)
!
!  Adapted from Heat conduction for constant value to
!  include temperature dependence to handle high temperatures
!  in hot diffuse cores of SN remnants in interstellar chi propto sqrt(T)
!  This routine also adds in turbulent diffusion, if chi_t /= 0.
!  Ds/Dt = ... + 1/(rho*T) grad(flux), where
!  flux = chi*rho*gradT + chi_t*rho*T*grads
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
      intent(inout) :: df
!
!  Check that chi is ok.
!
      if (headtt) print*,'calc_heatcond_cspeed_chi: chi=',chi
!
!  Heat conduction
!  Note: these routines require revision when ionization turned on
!  The variable g2 is reused to calculate glnP.gss a few lines below.
!
!  diffusion of the form:
!  rho*T*Ds/Dt = ... + nab.(rho*cp*chi*gradT)
!        Ds/Dt = ... + cp*chi*[del2lnTT+(glnrho+glnTT).glnTT]
!
!  Now chi=chi_0 sqrt(T) so additional 0.5*glnTT.glnTT
!
!  with additional turbulent diffusion
!  rho*T*Ds/Dt = ... + nab.(rho*T*chit*grads)
!        Ds/Dt = ... + chit*[del2ss+(glnrho+glnTT).gss]
!
!  Note: need thermally sensitive diffusion without magnetic field
!  for interstellar hydro runs to contrain SNr core temp
!  fred: 23.9.17 replaced 0.5 with chi_cspeed so exponent can be generalised
!
      thchi=chi*exp(chi_cspeed*p%lnTT)
      if (pretend_lnTT) then
        call dot(p%glnrho+(1.+chi_cspeed)*p%glnTT,p%glnTT,g2)
        thdiff=gamma*thchi*(p%del2lnTT+g2)
        if (chi_t/=0.) then
          call dot(p%glnrho+p%glnTT,p%gss,g2)
          thdiff=thdiff+chi_t*(p%del2ss+g2)
        endif
      else
        call dot(p%glnrho+(1.+chi_cspeed)*p%glnTT,p%glnTT,g2)
        thdiff=thchi*(p%del2lnTT+g2)*p%cp
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
      if (headtt) print*,'calc_heatcond_cspeed_chi: added thdiff'
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
!        if (ldiagnos.and.idiag_dtchi/=0) then
!          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
!        endif
      endif
!
    endsubroutine calc_heatcond_cspeed_chi
!***********************************************************************
    subroutine calc_heatcond_sqrtrhochi(df,p)
!
!  Adapted from Heat conduction for constant value to
!  include temperature dependence to handle high temperatures
!  in hot diffuse cores of SN remnants in interstellar chi propto sqrt(rho)
!  This routine also adds in turbulent diffusion, if chi_t /= 0.
!  Ds/Dt = ... + 1/(rho*T) grad(flux), where
!  flux = chi*rho*gradT + chi_t*rho*T*grads
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
      intent(inout) :: df
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
!        if (ldiagnos.and.idiag_dtchi/=0) then
!          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
!        endif
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
      intent(inout) :: df
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
      real, dimension (mx,my,mz,mvar),intent(inout) :: df
!
      real, dimension (nx) :: thdiff
!
!  Check that chi_hyper3_aniso is ok.
!
      if (headtt) print*, &
           'calc_heatcond_hyper3_aniso: chi_hyper3_aniso=', &
           chi_hyper3_aniso
!
!  Heat conduction.
!
      call del6fj(f,chi_hyper3_aniso,iss,thdiff)
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
!
      intent(in)  :: f
      intent(inout) :: df
!
      integer :: j
      real, dimension (nx) :: thdiff,tmp
!
      if (headtt) print*, 'calc_heatcond_hyper3_polar: chi_hyper3=', chi_hyper3
!
      thdiff=0.
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
           diffus_chi3=diffus_chi3+chi_hyper3*pi4_1*dxmin_pencil**4
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
!
      intent(in)  :: f
      intent(inout) :: df
!
      real, dimension (nx) :: thdiff,tmp,advec_hypermesh_ss
      integer :: j
!
      if (headtt) print*, 'calc_heatcond_hyper3_mesh: chi_hyper3=', chi_hyper3
!
      thdiff=0.0
      do j=1,3
        call der6(f,iss,tmp,j,IGNOREDX=.true.)
        if (ldynamical_diffusion) then
          thdiff = thdiff + chi_hyper3_mesh * tmp * dline_1(:,j)
        else
          thdiff = thdiff + chi_hyper3_mesh*pi5_1/60.*tmp*dline_1(:,j)
        endif
      enddo
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
!
      if (headtt) print*,'calc_heatcond_hyper3: added thdiff'
!
      if (lfirst.and.ldt) then
        if (ldynamical_diffusion) then
          diffus_chi3 = diffus_chi3 + chi_hyper3_mesh * (abs(dline_1(:,1)) + abs(dline_1(:,2)) + abs(dline_1(:,3)))
          advec_hypermesh_ss = 0.0
        else
          advec_hypermesh_ss=chi_hyper3_mesh*pi5_1*sqrt(dxyz_2)
        endif
        advec2_hypermesh=advec2_hypermesh+advec_hypermesh_ss**2
      endif
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
!  22-aug-14/joern: removed chi_t from here.
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: thdiff,g2,gshockglnTT,gshockgss
!
      intent(in) :: p
      intent(inout) :: df
!
!  Check that chi is ok.
!
      if (headtt) print*,'calc_heatcond_shock: chi_shock=',chi_shock
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
      if (lheatc_shock) then
        if (pretend_lnTT) then
          thdiff=gamma*chi_shock*(p%shock*(p%del2lnrho+g2)+gshockglnTT)
        else
          if (lchi_shock_density_dep) then
            call dot(p%gshock,p%gss,gshockgss)
            call dot(twothird*p%glnrho+p%glnTT,p%gss,g2)
            thdiff=exp(-onethird*p%lnrho)*chi_shock* &
                   (p%shock*(exp(twothird*p%lnrho)*p%del2ss+g2)+gshockgss)
          else
            thdiff=chi_shock*(p%shock*(p%del2lnTT+g2)+gshockglnTT)
          endif
        endif
      endif
      if (lheatc_shock2) then
        if (pretend_lnTT) then
          thdiff=gamma*chi_shock2*(p%shock**2*(p%del2lnrho+g2)+gshockglnTT)
        else
          if (lchi_shock_density_dep) then
            call dot(p%gshock,p%gss,gshockgss)
            call dot(twothird*p%glnrho+p%glnTT,p%gss,g2)
            thdiff=exp(-onethird*p%lnrho)*chi_shock2* &
                   (p%shock**2*(exp(twothird*p%lnrho)*p%del2ss+g2)+2*p%shock*gshockgss)
          else
            thdiff=chi_shock2*(p%shock**2*(p%del2lnTT+g2)+2*p%shock*gshockglnTT)
          endif
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
          if (lchi_shock_density_dep) then
            if (lheatc_shock) &
                diffus_chi=diffus_chi+exp(-onethird*p%lnrho)*chi_shock*p%shock*p%cp1*dxyz_2
            if (lheatc_shock2) &
                diffus_chi=diffus_chi+exp(-onethird*p%lnrho)*chi_shock2*p%shock**2*p%cp1*dxyz_2
          else
            if (lheatc_shock) &
                diffus_chi=diffus_chi+(gamma*chi_shock*p%shock)*dxyz_2
            if (lheatc_shock2) &
                diffus_chi=diffus_chi+(gamma*chi_shock2*p%shock**2)*dxyz_2
          endif
        else
          if (lheatc_shock) &
              diffus_chi=diffus_chi+(chi_shock*p%shock)*dxyz_2
          if (lheatc_shock2) &
              diffus_chi=diffus_chi+(chi_shock2*p%shock**2)*dxyz_2
        endif
!        if (ldiagnos.and.idiag_dtchi/=0) then
!          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
!        endif
      endif
!
    endsubroutine calc_heatcond_shock
!***********************************************************************
    subroutine calc_heatcond_shock_profr(df,p)
!.
!  21-aug-14/joern: adapted from calc_heatcond_shock
!                   but no chit adding anymore
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot, step, der_step
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: thdiff,g2,gshockglnTT,gchiglnTT,pchi_shock
      real, dimension (nx,3) :: gradchi_shock
!
      intent(in) :: p
      intent(inout) :: df
!
!  Check that chi is ok.
!
      if (headtt) print*,'calc_heatcond_shock: chi_shock,chi_jump_shock=', &
                         chi_shock,chi_jump_shock
!
      pchi_shock = chi_shock + chi_shock*(chi_jump_shock-1.)* &
                      step(p%r_mn,xchi_shock,widthchi_shock)
!
          gradchi_shock(:,1) = chi_shock*(chi_jump_shock-1.)* &
                            der_step(p%r_mn,xchi_shock,widthchi_shock)
          gradchi_shock(:,2) = 0.
          gradchi_shock(:,3) = 0.

      call dot(p%gshock,p%glnTT,gshockglnTT)
      call dot(p%glnrho+p%glnTT,p%glnTT,g2)
!
!  Shock entropy diffusivity with a radial profile
!  Write: chi_shock = pchi_shock*shock, gshock=grad(shock)
!  and gradchi_shock=grad(pchi_shock) so
!  Ds/Dt = ... + pchi_shock*[shock*(del2ss+glnTT.gss+glnrho.gss) + gshock.gss]
!              + shock*gradchi_shock.gss
!
      if (headtt) print*,'calc_heatcond_shock_profr: use shock diffusion with radial profile'

      call dot(p%gshock,p%glnTT,gshockglnTT)
      call dot(p%glnrho+p%glnTT,p%glnTT,g2)
      call dot(gradchi_shock,p%glnTT,gchiglnTT)

!
      if (pretend_lnTT) then
        thdiff=gamma*(pchi_shock*(p%shock*(p%del2lnrho+g2)+gshockglnTT)+gchiglnTT*p%shock)
      else
        thdiff=pchi_shock*(p%shock*(p%del2lnTT+g2)+gshockglnTT)+gchiglnTT*p%shock
      endif
!
!  Add heat conduction to entropy equation.
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
      if (headtt) print*,'calc_heatcond_shock_profr: added thdiff'
!
!  Check maximum diffusion from thermal diffusion.
!  With heat conduction, the second-order term for entropy is
!  gamma*pchi*del2ss.
!
      if (lfirst.and.ldt) then
        if (leos_idealgas) then
          diffus_chi=diffus_chi+(gamma*pchi_shock*p%shock)*dxyz_2
        else
          diffus_chi=diffus_chi+(pchi_shock*p%shock)*dxyz_2
        endif
!        if (ldiagnos.and.idiag_dtchi/=0) then
!          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
!        endif
      endif
!
    endsubroutine calc_heatcond_shock_profr
!***********************************************************************
    subroutine calc_heatcond_constK(df,p)
!
!  Heat conduction.
!
!   8-jul-02/axel: adapted from Wolfgang's more complex version
!  30-mar-06/ngrs: simplified calculations using p%glnTT and p%del2lnTT
!
      use Diagnostics, only: max_mn_name, yzsum_mn_name_x
      use Sub, only: dot
!
      type (pencil_case) :: p
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: chix
      real, dimension (nx) :: thdiff,g2
      real, dimension (nx) :: hcond
!
      intent(in) :: p
      intent(inout) :: df
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
      chix = p%rho1*hcond*p%cp1
!
!  Put empirical heat transport suppression by the B-field.
!
      if (chiB/=0.) chix = chix/(1.+chiB*p%b2)
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
!  1-D averages.
!
      if (l1davgfirst) then
        if (idiag_fradmx/=0) call yzsum_mn_name_x(-hcond*p%TT*p%glnTT(:,1),idiag_fradmx)
      endif
!
!  Check maximum diffusion from thermal diffusion.
!  With heat conduction, the second-order term for entropy is
!  gamma*chix*del2ss.
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chix*dxyz_2
!        if (ldiagnos.and.idiag_dtchi/=0) then
!          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
!        endif
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
      use EquationOfState, only: gamma
      use Debug_IO, only: output_pencil
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
      intent(inout) :: df
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
      call multsv_mn(tmps,p%glnTT,gvKpara)!
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
!  to magnetic field lines.
!
!  24-aug-09/bing: moved from denergy_dt to here
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
!        if (ldiagnos.and.idiag_dtchi/=0) then
!          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
!        endif
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
    subroutine calc_heatcond_kramers(f,df,p)
!
!  Heat conduction using Kramers' opacity law
!
!  23-feb-11/pete: coded
!  24-aug-15/MR: bounds for chi introduced
!
      use Diagnostics
      use Sub, only: dot
      use General, only: notanumber
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: p
      intent(inout) :: df
!
      real, dimension(nx) :: thdiff, chix, g2
      real, dimension(nx) :: Krho1, chit_prof, del2ss1
      real, dimension(nx,3) :: gradchit_prof, gss1
      integer :: j
!
!  Diffusion of the form
!      rho*T*Ds/Dt = ... + nab.(K*gradT),
!  where
!      K = K_0*(T**6.5/rho**2)**n.
!  In reality n=1, but we may need to use n\=1 for numerical reasons.
!
      Krho1 = hcond0_kramers*p%rho1**(2.*nkramers+1.)*p%TT**(6.5*nkramers)   ! = K/rho
      !Krho1 = hcond0_kramers*exp(-p%lnrho*(2.*nkramers+1.)+p%lnTT*(6.5*nkramers))   ! = K/rho
      if (chimax_kramers>0.) &
        Krho1 = max(min(Krho1,chimax_kramers/p%cp1),chimin_kramers/p%cp1)
      call dot(-2.*nkramers*p%glnrho+(6.5*nkramers+1)*p%glnTT,p%glnTT,g2)
      thdiff = Krho1*(p%del2lnTT+g2)
!
      if (chi_t/=0.0) then
!
        if (lcalc_ssmean) then
          do j=1,3; gss1(:,j)=p%gss(:,j)-gssmz(n-n1+1,j); enddo
          del2ss1=p%del2ss-del2ssmz(n-n1+1)
        else if (lcalc_ssmeanxy) then
          do j=1,3; gss1(:,j)=p%gss(:,j)-gssmx(:,j); enddo
          del2ss1=p%del2ss-del2ssmx
        else
          do j=1,3; gss1(:,j)=p%gss(:,j) ; enddo
          del2ss1=p%del2ss
        endif
        call dot(p%glnrho+p%glnTT,gss1,g2)
        call get_prof_pencil(chit_prof,gradchit_prof,lsphere_in_a_box,1.,chit_prof1,chit_prof2,xbot,xtop,p,f, &
                             stored_prof=chit_prof_stored,stored_dprof=dchit_prof_stored)
        thdiff=thdiff+chi_t*chit_prof*(del2ss1+g2)
        call dot(gradchit_prof,gss1,g2)
        thdiff=thdiff+chi_t*g2
      endif
!
      if (pretend_lnTT) thdiff = p%cv1*thdiff
!
!  Here chix = K/(cp rho) is needed for diffus_chi calculation.
!
      chix = p%cp1*Krho1
    
      if (ldiagnos.or.l1davgfirst.or.l2davgfirst) Krho1 = Krho1*p%rho      ! now Krho1=K

      if (ldiagnos) call sum_mn_name(Krho1,idiag_Kkramersm)
!
!  Write radiative flux array.
!
      if (l1davgfirst) then
        if (idiag_fradz_kramers/=0) call xysum_mn_name_z(-Krho1*p%TT*p%glnTT(:,3),idiag_fradz_kramers)
        call xysum_mn_name_z( Krho1, idiag_Kkramersmz)
        call yzsum_mn_name_x( Krho1, idiag_Kkramersmx)
        if (idiag_fradx_kramers/=0) call yzsum_mn_name_x(-Krho1*p%rho*p%TT*p%glnTT(:,1),idiag_fradx_kramers)
        if (idiag_fturbz/=0) call xysum_mn_name_z(-chi_t*chit_prof*p%rho*p%TT*p%gss(:,3),idiag_fturbz)
      endif
!
!  2d-averages
!
      if (l2davgfirst) then
        if (idiag_fradxy_kramers/=0) call zsum_mn_name_xy(-Krho1*p%TT*p%glnTT(:,1),idiag_fradxy_kramers)
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
          call fatal_error_local('calc_heatcond_kramers','NaNs in thdiff')
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
!        if (ldiagnos.and.idiag_dtchi/=0) then
!          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
!        endif
      endif
!
    endsubroutine calc_heatcond_kramers
!***********************************************************************
    subroutine calc_heatcond_smagorinsky(f,df,p)
!
!  Heat conduction using Smagorinsky viscosity
!
!  21-apr-17/pete: coded
!
      use Diagnostics
      use Debug_IO, only: output_pencil
      use Sub, only: dot, grad
      use General, only: notanumber
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (nx) :: thdiff, chix, g2
!
      intent(in) :: p
      intent(inout) :: f,df
!
      real, dimension(nx) :: del2ss1
      real, dimension(nx,3) :: gss1, gradnu
      integer :: j
!
      thdiff=0.
!
!  Diffusion of the form
!     rho*T*Ds/Dt = ... + nab.((nu_smag/Pr_smag)*rho*T*grads),
!
      if (lcalc_ssmean) then
        do j=1,3; gss1(:,j)=p%gss(:,j)-spread(gssmz(n-n1+1,j), 1, l2-l1+1); enddo
        del2ss1=p%del2ss-spread(del2ssmz(n-n1+1), 1, l2-l1+1)
      else if (lcalc_ssmeanxy) then
        do j=1,3; gss1(:,j)=p%gss(:,j)-gssmx(:,j); enddo
        del2ss1=p%del2ss-del2ssmx
      else
        do j=1,3; gss1(:,j)=p%gss(:,j) ; enddo
        del2ss1=p%del2ss
      endif
      if (lchit_noT) then
        call dot(p%glnrho,gss1,g2)
      else
        call dot(p%glnrho+p%glnTT,gss1,g2)
      endif
      thdiff=Pr_smag1*p%nu_smag*(del2ss1+g2)
!
      call grad(f,inusmag,gradnu)
!
      call dot(gradnu,gss1,g2)
      thdiff=thdiff+Pr_smag1*g2
!
!  Here chix = nu_smag/Pr_smag is needed for diffus_chi calculation.
!
      chix = Pr_smag1*p%nu_smag
!
!  Write radiative flux array.
!
      if (l1davgfirst) then
        call xysum_mn_name_z(-Pr_smag1*p%nu_smag*p%rho*p%TT*gss1(:,3),idiag_fturbz)
      endif
!
!  2d-averages
!
!      if (l2davgfirst) then
!        if (idiag_fradxy_kramers/=0) call zsum_mn_name_xy(-Krho1*p%TT*p%glnTT(:,1),idiag_fradxy_kramers)
!      endif
!
!  Check for NaNs initially.
!
      if (headt .and. (Pr_smag1/=0.0)) then
        if (notanumber(p%rho1))    print*,'calc_heatcond_smagorinsky: NaNs in rho1'
        if (notanumber(chix))      print*,'calc_heatcond_smagorinsky: NaNs in chix'
        if (notanumber(p%del2ss))  print*,'calc_heatcond_smagorinsky: NaNs in del2ss'
        if (notanumber(p%TT))      print*,'calc_heatcond_smagorinsky: NaNs in TT'
        if (notanumber(p%glnTT))   print*,'calc_heatcond_smagorinsky: NaNs in glnT'
        if (notanumber(g2))        print*,'calc_heatcond_smagorinsky: NaNs in g2'
        if (notanumber(thdiff))    print*,'calc_heatcond_smagorinsky: NaNs in thdiff'
!
!  Most of these should trigger the following trap.
!
        if (notanumber(thdiff)) then
          print*, 'calc_heatcond_smagorinsky: m,n,y(m),z(n)=', m, n, y(m), z(n)
          call fatal_error_local('calc_heatcond_smagorinsky','NaNs in thdiff')
        endif
      endif
!
!  At the end of this routine, add all contribution to
!  thermal diffusion on the rhs of the entropy equation,
!  so Ds/Dt = ... + thdiff = ... + (...)/(rho*T)
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
!
      if (headtt) print*,'calc_heatcond_smagorinsky: added thdiff'
!
!  Check maximum diffusion from thermal diffusion.
!  NB: With heat conduction, the second-order term for entropy is
!    gamma*chix*del2ss.
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+(gamma*chix)*dxyz_2
      endif
!
    endsubroutine calc_heatcond_smagorinsky
!***********************************************************************
    subroutine calc_heatcond(f,df,p)
!
!  In this routine general heat conduction profiles are being provided
!  and applied to the entropy equation.
!
!  17-sep-01/axel: coded
!  14-jul-05/axel: corrected expression for chi_t diffusion.
!  30-mar-06/ngrs: simplified calculations using p%glnTT and p%del2lnTT
!
      use Diagnostics
      use Debug_IO, only: output_pencil
      use Sub, only: dot, g2ij
      use General, only: notanumber
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df

      real, dimension (nx,3) :: glnThcond,glhc,gss1
      real, dimension (nx) :: chix,chit_prof
      real, dimension (nx) :: thdiff,g2,del2ss1
      real, dimension (nx) :: glnrhoglnT
      real, dimension (nx) :: hcond
      real, dimension (nx,3) :: gradchit_prof
      real, dimension (nx,3,3) :: tmp
      !real, save :: z_prev=-1.23e20
      real :: s2,c2,sc
      integer :: j
!
!  Heat conduction.
!
      if (hcond0/=0..or.lread_hcond) then
        if (headtt) then
          print*,'calc_heatcond: hcond0=',hcond0
          print*,'calc_heatcond: lgravz=',lgravz
          if (lgravz) print*,'calc_heatcond: Fbot,Ftop=',Fbot,Ftop
        endif

        call get_prof_pencil(hcond,glhc,lsphere_in_a_box.or.lcylinder_in_a_box, &
                             hcond0,hcond1,hcond2,r_bcz,r_ext,p,f,hcond_prof,dlnhcond_prof,llog=.true.)
!
!  Diffusion of the form
!
!  rho*T*Ds/Dt = ... + nab.(K*gradT)
!        Ds/Dt = ... + K/rho*[del2lnTT+(glnTT+glnhcond).glnTT]
!
!  where chix = K/(cp rho) is needed for diffus_chi calculation.
!
        if (notanumber(p%glnTT)) &
          call fatal_error_local('calc_heatcond', 'NaNs in p%glnTT')

        chix = p%rho1*hcond*p%cp1
        glnThcond = p%glnTT + glhc    ! grad ln(T*hcond)
        call dot(p%glnTT,glnThcond,g2)
        if (pretend_lnTT) then
          thdiff = p%cv1*p%rho1*hcond * (p%del2lnTT + g2)
        else
          if (lhcond0_density_dep) then
            call dot(p%glnTT,p%glnrho,glnrhoglnT)
            thdiff = sqrt(p%rho1)*hcond * (p%del2lnTT + g2+0.5*glnrhoglnT)
            chix = sqrt(p%rho1)*hcond*p%cp1
          else
            thdiff = p%rho1*hcond * (p%del2lnTT + g2)
          endif
        endif
!
!DM+GG This is done earlier now. The profile writing done earlier in this
! code also includes the ghost zones. This commented line may stay longer
! than the ones above.
!      if (lgravz) call write_zprof('hcond',hcond)
!
!  Write radiative flux array.
!
        if (l1davgfirst) then
          if (idiag_fradz_Kprof/=0) call xysum_mn_name_z(-hcond*p%TT*p%glnTT(:,3),idiag_fradz_Kprof)
          if (idiag_fradmx/=0) call yzsum_mn_name_x(-hcond*p%TT*p%glnTT(:,1),idiag_fradmx)
        endif
!
!  2d-averages
!
        if (l2davgfirst) then
          if (idiag_fradxy_Kprof/=0) call zsum_mn_name_xy(-hcond*p%TT*p%glnTT(:,1),idiag_fradxy_Kprof)
          if (idiag_fradymxy_Kprof/=0) call zsum_mn_name_xy(p%glnTT,idiag_fradymxy_Kprof,(/0,1,0/),-hcond*p%TT)
        endif
      endif  ! hcond0/=0.
!
!  "Turbulent" entropy diffusion.
!
!  Traditionally only present if g.gradss > 0 (unstable stratification)
!  but this is not currently being checked.
!
      if (chi_t/=0.) then
        if (headtt) &
          print*,'calc_headcond: "turbulent" entropy diffusion: chi_t=',chi_t
        call get_prof_pencil(chit_prof,gradchit_prof,lsphere_in_a_box,1.,chit_prof1,chit_prof2,xbot,xtop,p,f, &
                             stored_prof=chit_prof_stored,stored_dprof=dchit_prof_stored)
!
!  'Standard' formulation where the flux contains temperature:
! 
!  ... + div(rho*T*chi_t*grads) = ... + chi_t*[del2s + (glnrho+glnTT+glnchi).grads]
!
!  Alternative (correct?) formulation where temperature is absent from the flux:
!
!  ... + div(rho*chi_t*grads) = ... + chi_t*[del2s + (glnrho+glnchi).grads]
!
!  If lcalc_ssmean=T or lcalc_ssmeanxy=T, mean stratification is subtracted.
!
        if (lcalc_ssmean .or. lcalc_ssmeanxy) then
          if (lcalc_ssmean) then
            do j=1,3; gss1(:,j)=p%gss(:,j)-gssmz(n-n1+1,j); enddo
            del2ss1=p%del2ss-del2ssmz(n-n1+1)
          else if (lcalc_ssmeanxy) then
            do j=1,3; gss1(:,j)=p%gss(:,j)-gssmx(:,j); enddo
            del2ss1=p%del2ss-del2ssmx
          endif
!  
          if (lchit_noT) then
            call dot(p%glnrho,gss1,g2)
          else
            call dot(p%glnrho+p%glnTT,gss1,g2)
          endif
!
          thdiff=thdiff+chi_t*chit_prof*(del2ss1+g2)
          call dot(gradchit_prof,gss1,g2)
          thdiff=thdiff+chi_t*g2
        else
          if (lchit_noT) then
            call dot(p%glnrho,p%gss,g2)
          else
            call dot(p%glnrho+p%glnTT,p%gss,g2)
          endif
!
          thdiff=thdiff+chi_t*chit_prof*(p%del2ss+g2)
          call dot(gradchit_prof,p%gss,g2)
          thdiff=thdiff+chi_t*g2
        endif
        if (l1davgfirst) then
          if (idiag_fturbz/=0) call xysum_mn_name_z(-chi_t*chit_prof*p%rho*p%TT*p%gss(:,3),idiag_fturbz)
          if (idiag_fturbmx/=0) call yzsum_mn_name_x(-chi_t*chit_prof*p%rho*p%TT*p%gss(:,1),idiag_fturbmx)
        endif
        if (l2davgfirst) then
          if (idiag_fturbxy/=0  ) call zsum_mn_name_xy(-chi_t*chit_prof*p%rho*p%TT*p%gss(:,1),idiag_fturbxy)
          if (idiag_fturbymxy/=0) call zsum_mn_name_xy(p%gss,idiag_fturbymxy,(/0,1,0/),-chi_t*chit_prof*p%rho*p%TT)
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
          if (headtt) &
            print*, 'calc_heatcond: '// &
                'anisotropic "turbulent" entropy diffusion: chit_aniso=', &
                 chit_aniso
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
            thdiff=thdiff+chit_aniso_prof* &
                ((1.-3.*c2)*r1_mn*p%gss(:,1) + &
                (-sc*tmp(:,1,2))+(p%glnrho(:,2)+p%glnTT(:,2))*(-sc*p%gss(:,1)))
          else
!
!  Otherwise use the full formulation.
!
            thdiff=thdiff + &
                ((dchit_aniso_prof*c2+chit_aniso_prof/x(l1:l2))*p%gss(:,1)+ &
                sc*(chit_aniso_prof/x(l1:l2)-dchit_aniso_prof)*p%gss(:,2))
!
            call g2ij(f,iss,tmp)
            thdiff=thdiff+chit_aniso_prof* &
                ((-sc*(tmp(:,1,2)+tmp(:,2,1))+c2*tmp(:,1,1)+s2*tmp(:,2,2))+ &
                ((p%glnrho(:,1)+p%glnTT(:,1))*(c2*p%gss(:,1)-sc*p%gss(:,2))+ &
                ( p%glnrho(:,2)+p%glnTT(:,2))*(s2*p%gss(:,2)-sc*p%gss(:,1))))
!
          endif
        endif
      endif
!
!  Check for NaNs initially.
!
      if (hcond0/=0..or.lread_hcond.or.chi_t/=0.) then
        if (headt) then
!
          if (notanumber(p%rho1))    print*,'calc_heatcond: NaNs in rho1'
          if (chi_t/=0.) then
            if (notanumber(p%del2ss))  print*,'calc_heatcond: NaNs in del2ss'
          endif
          if (hcond0/=0..or.lread_hcond) then
            if (notanumber(hcond))     print*,'calc_heatcond: NaNs in hcond'
            if (notanumber(1/hcond))   print*,'calc_heatcond: NaNs in 1/hcond'
            if (notanumber(glhc))      print*,'calc_heatcond: NaNs in glhc'
            if (notanumber(chix))      print*,'calc_heatcond: NaNs in chix'
            if (notanumber(glnThcond)) print*,'calc_heatcond: NaNs in glnThcond'
          endif
          if (notanumber(g2))        print*,'calc_heatcond: NaNs in g2'
!
!  Most of these should trigger the following trap.
!
          if (notanumber(thdiff)) then
            print*,'calc_heatcond: NaNs in thdiff'
            print*, 'calc_heatcond: m,n,y(m),z(n)=', m, n, y(m), z(n)
            call fatal_error_local('calc_heatcond','NaNs in thdiff')
          endif
        endif
        if ((hcond0/=0..or.lread_hcond).and.lwrite_prof .and. ip<=9) then
          call output_pencil('chi.dat',chix,1)
          call output_pencil('hcond.dat',hcond,1)
          call output_pencil('glhc.dat',glhc,3)
        endif

        if (headt .and. lfirst .and. ip == 13) &
           call output_pencil('heatcond.dat',thdiff,1)
!
!  At the end of this routine, add all contribution to
!  thermal diffusion on the rhs of the entropy equation,
!  so Ds/Dt = ... + thdiff = ... + (...)/(rho*T)
!
        df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
!
        if (headtt) print*,'calc_heatcond: added thdiff'
!
      endif
!
!  Check maximum diffusion from thermal diffusion.
!  NB: With heat conduction, the second-order term for entropy is
!    gamma*chix*del2ss.
!
      if (lfirst.and.ldt) then
        if (hcond0/=0..or.lread_hcond) &
          diffus_chi=diffus_chi+gamma*chix*dxyz_2
        if (chi_t/=0.) &
          diffus_chi=diffus_chi+chi_t*chit_prof*dxyz_2
!        if (ldiagnos.and.idiag_dtchi/=0) &
!          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
      endif
!
    endsubroutine calc_heatcond
!***********************************************************************
    subroutine calc_heatcond_sfluct(df,p)
!
!  In this routine general heat conduction profiles are being provided.
!
!   9-jun-16/axel: adapted from calc_heatcond
!
      use Sub, only: dot
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: gss1
      real, dimension (nx) :: thdiff,g2,del2ss1
      integer :: j
!
      intent(in) :: p
      intent(inout) :: df
!
!  Entropy diffusion on fluctuations only
!
      if (headtt) print*,'calc_heatcond_sfluct: chi_t=',chi_t
!
!  "Turbulent" entropy diffusion (operates on entropy fluctuations only).
!
      if (chi_t/=0.) then
        if (headtt) &
          print*,'calc_headcond: "turbulent" entropy diffusion: chi_t=',chi_t
!
!  ... + div(rho*T*chi*grads) = ... + chi*[del2s + (glnrho+glnTT+glnchi).grads]
!  If lcalc_ssmean=T or lcalc_ssmeanxy=T, mean stratification is subtracted.
!
        if (lcalc_ssmean .or. lcalc_ssmeanxy) then
          if (lcalc_ssmean) then
            do j=1,3; gss1(:,j)=p%gss(:,j)-gssmz(n-n1+1,j); enddo
            del2ss1=p%del2ss-del2ssmz(n-n1+1)
          else if (lcalc_ssmeanxy) then
            do j=1,3; gss1(:,j)=p%gss(:,j)-gssmx(:,j); enddo
            del2ss1=p%del2ss-del2ssmx
          endif
          call dot(p%glnrho+p%glnTT,gss1,g2)
          thdiff=chi_t*(del2ss1+g2)
        else
          call fatal_error("calc_heatcond_sfluct","lcalc_ssmean(xy) needed")
        endif
      else
        call fatal_error("calc_heatcond_sfluct","chi_t must not be 0")
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
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+chi_t*dxyz_2
      endif
!
    endsubroutine calc_heatcond_sfluct
!***********************************************************************
    subroutine calc_heatcond_chit(f,df,p)
!
!  Subgrid-scale ("turbulent") diffusion of entropy.
!
!  Three variants: 1) chi_t acts on total ss
!                  2)  -----||----- suitably averaged ss
!                  3)  -----||----- fluctuations of ss
!  All variants are standalone but 2) and 3) can also occur simultaneously
!  with different profiles and/or magnitudes for chi_t[0,1].
!
!   21-apr-17/pete: adapted from calc_heatcond_sfluct
!
      use Sub, only: dot, grad, del2
      use Diagnostics
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: gss0, gss1
      real, dimension (nx) :: thdiff, g2, del2ss0, del2ss1
      real, dimension (nx,3) :: gradchit_prof, gradchit_prof_fluct
      real, dimension (nx) :: chit_prof, chit_prof_fluct
      integer :: j
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  General turbulent entropy diffusion
!
      if (headtt) then
        print*,'calc_heatcond_chit: chi_t0=',chi_t0
        print*,'calc_heatcond_chit: chi_t1=',chi_t1
      endif

      thdiff=0.

      if (chi_t0/=0.) then
!
!  Compute profiles. Diffusion of total and mean use the same profile
!
        if (lchit_total .or. lchit_mean) &
          call get_prof_pencil(chit_prof,gradchit_prof,lsphere_in_a_box,1.,chit_prof1,chit_prof2,xbot,xtop,p,f, &
                               stored_prof=chit_prof_stored,stored_dprof=dchit_prof_stored)
!
!  chi_t acts on the total entropy 
!      
        if (lchit_total) then
          if (lchit_noT) then 
            call dot(p%glnrho,p%gss,g2)
          else
            call dot(p%glnrho+p%glnTT,p%gss,g2)
          endif
          thdiff=thdiff+chi_t0*chit_prof*(p%del2ss+g2)
          call dot(gradchit_prof,p%gss,g2)
          thdiff=thdiff+chi_t0*g2
        endif
!
!  chi_t acts on averaged entropy, need to have either lcalc_ssmean=T
!  or lcalc_ssmeanxy=T
!
        if (lchit_mean .and. (lcalc_ssmean .or. lcalc_ssmeanxy)) then
          if (lcalc_ssmean) then
            do j=1,3; gss0(:,j)=spread(gssmz(n-n1+1,j), 1, l2-l1+1); enddo
            del2ss0=spread(del2ssmz(n-n1+1), 1, l2-l1+1)
          else if (lcalc_ssmeanxy) then
            do j=1,3; gss0(:,j)=gssmx(:,j); enddo
            del2ss0=del2ssmx
          endif
          if (lchit_noT) then 
            call dot(p%glnrho,gss0,g2)
          else
            call dot(p%glnrho+p%glnTT,gss0,g2)
          endif
          thdiff=thdiff+chi_t0*chit_prof*(del2ss0+g2)
          call dot(gradchit_prof,gss0,g2)
          thdiff=thdiff+chi_t0*g2
        endif
        if (l1davgfirst) then
          if (idiag_fturbtz/=0) call xysum_mn_name_z(-chi_t0*chit_prof*p%rho*p%TT*p%gss(:,3),idiag_fturbtz)
          if (idiag_fturbmz/=0) call xysum_mn_name_z(-chi_t0*chit_prof*p%rho*p%TT*gss0(:,3),idiag_fturbmz)
        endif
      endif
!
!  chi_t acts on fluctuations of entropy, again need to have either 
!  lcalc_ssmean=T or lcalc_ssmeanxy=T
!
      if (lchit_fluct.and.chi_t1/=0.) then
        if (lcalc_ssmean .or. lcalc_ssmeanxy .or. lss_running_aver_as_aux) then
          if ((lgravr.or.lgravx.or.lgravz).and..not.(lsphere_in_a_box.or.lchi_t1_noprof)) then
            call get_prof_pencil(chit_prof_fluct,gradchit_prof_fluct,.false., &  ! no 2D/3D profiles of chit_fluct implemented
                                 stored_prof=chit_prof_fluct_stored,stored_dprof=dchit_prof_fluct_stored)
          else
            chit_prof_fluct=chi_t1; gradchit_prof_fluct=0.
          endif
        endif

        if (lcalc_ssmean .or. lcalc_ssmeanxy) then
          
          if (lcalc_ssmean) then
            do j=1,3; gss1(:,j)=p%gss(:,j)-gssmz(n-n1+1,j); enddo
            del2ss1=p%del2ss-del2ssmz(n-n1+1)
          else if (lcalc_ssmeanxy) then
            do j=1,3; gss1(:,j)=p%gss(:,j)-gssmx(:,j); enddo
            del2ss1=p%del2ss-del2ssmx
          endif
          if (lchit_noT) then 
            call dot(p%glnrho,gss1,g2)
          else
            call dot(p%glnrho+p%glnTT,gss1,g2)
          endif

          if ((lgravr.or.lgravx.or.lgravz).and..not.lsphere_in_a_box) then
            thdiff=thdiff+chit_prof_fluct*(del2ss1+g2)
            call dot(gradchit_prof_fluct,gss1,g2)
            thdiff=thdiff+g2
          endif
!
          if (l1davgfirst) then
            if (idiag_fturbfz/=0) call xysum_mn_name_z(-chit_prof_fluct*p%rho*p%TT*gss1(:,3),idiag_fturbfz)
          endif
        endif
!
! chi_t1 acting on deviations from a running 3D mean of entropy   
!
        if (lss_running_aver_as_aux) then
          call grad(f(:,:,:,iss_run_aver),gss1)
          do j=1,3
            gss0(:,j)=p%gss(:,j)-gss1(:,j)
          enddo
          call del2(f(:,:,:,iss_run_aver),del2ss1)
          del2ss0=p%del2ss-del2ss1
          if (lchit_noT) then 
            call dot(p%glnrho,gss0,g2)
          else
            call dot(p%glnrho+p%glnTT,gss0,g2)
          endif
          thdiff=thdiff+chit_prof_fluct*(del2ss0+g2)
          call dot(gradchit_prof_fluct,gss0,g2)
          thdiff=thdiff+g2
        endif
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
!
      if (lfirst.and.ldt) then
        if (chi_t0/=0..and.(lchit_total .or. lchit_mean)) &
          diffus_chi=diffus_chi+chi_t0*chit_prof*dxyz_2
        if (lcalc_ssmean .or. lcalc_ssmeanxy .or. lss_running_aver_as_aux) then
          if (chi_t1/=0..and.lchit_fluct) &
            diffus_chi=diffus_chi+chit_prof_fluct*dxyz_2
        endif
      endif
!
    endsubroutine calc_heatcond_chit
!***********************************************************************
    subroutine calc_heat_cool(df,p,Hmax)
!
!  Add combined heating and cooling to entropy equation.
!
!  02-jul-02/wolf: coded
!
      use Diagnostics, only: sum_mn_name, xysum_mn_name_z
      use EquationOfState, only: cs0, get_cp1, lnrho0, &
                                 gamma,gamma_m1
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: Hmax
!
      real, dimension (nx) :: heat, TT_drive
      real :: profile_buffer
      real :: xi,cp1,Rgas
!
      intent(in) :: p
      intent(inout) :: Hmax,df
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
      if (lgravz .and. (.not.lcooling_general) .and. &
          ( (luminosity/=0.) .or. (cool/=0.) ) ) &
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
!  Add heating with Gaussian profile.
!
      if (heat_gaussianz/=0.0) then
        heat=heat+heat_gaussianz*exp(-.5*z(n)**2/heat_gaussianz_sigma**2)/ &
          (sqrt(2*pi)*heat_gaussianz_sigma)
      endif
!
!  Add spatially uniform heating.
!
      if (heat_uniform/=0.0) then
        if (zheat_uniform_range/=0.) then
          if (abs(z(n)) <= zheat_uniform_range) heat=heat+heat_uniform
        else
          heat=heat+heat_uniform
        endif
      endif
!
!  Add spatially uniform cooling in Ds/Dt equation.
!
      if (cool_uniform/=0.0) heat=heat-cool_uniform*p%rho*p%cp*p%TT
      if (cool_newton/=0.0) then
        if (lchromospheric_cooling) then
          call get_cp1(cp1)
          Rgas=(1.0-gamma1)/cp1
          xi=(z(n)-z_cor)/(xyz1(3)-z_cor)
          profile_buffer=xi**2*(3-2*xi)
          TT_drive=(cs0**2*(exp(gamma*cp1*p%initss+&
                   gamma_m1*(p%initlnrho-lnrho0))))/(gamma*Rgas)
          if (z(n).gt.-1.0) &
          heat=heat-cool_newton*p%cp*p%rho*profile_buffer*(p%TT-TT_drive)
        else
          heat=heat-cool_newton*p%TT**4
        endif
      endif
!
!  Add cooling with constant time-scale to TTref_cool.
!
      if (tau_cool/=0.0) then
        if (.not.ltau_cool_variable) then
          if (lphotoelectric_heating) then
            TT_drive=TTref_cool*p%rhop
          else
            TT_drive=TTref_cool
          endif
          if (TT_floor /= impossible) call apply_floor(TT_drive)
          heat=heat-p%rho*p%cp*gamma1*(p%TT-TT_drive)/tau_cool
        else
          call calc_heat_cool_variable(heat,p)
        endif
      endif
!
      if (lphotoelectric_heating_radius) then
        if (lcylindrical_coords .and. nzgrid==1) heat=heat+peh_factor*p%peh*dx_1(l1:l2)*dy_1(m)/p%r_mn
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
!  AB: I guess this should be disabled soon, in favor of using
!      lcooling_general=T, cooling_profile='cubic_step', cooltype='corona'
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
      if (limpose_heat_ceiling) where(heat>heat_ceiling) heat=heat_ceiling
!
!  Add heating/cooling to entropy equation.
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + p%TT1*p%rho1*heat
      if (lfirst.and.ldt) Hmax=Hmax+heat*p%rho1
!
!  Volume heating/cooling term on the mean entropy with respect to ss_const.
!  Allow for local heating/cooling w.r.t. ss when lcalc_ss_volaverage=F
!
      if (tau_cool_ss/=0.) then
        if (lcalc_ss_volaverage) then
          df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
            -(ss_volaverage-ss_const)/tau_cool_ss
        elseif (lcooling_ss_mz) then
          df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)-(p%ss-ss_mz(n))/tau_cool_ss
        else
          df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)-(p%ss-ss_const)/tau_cool_ss
        endif
      endif
!
!  Heating/cooling related diagnostics.
!
      if (ldiagnos) call sum_mn_name(heat,idiag_heatm)
!
!  Write divergence of cooling flux.
!
      if (l1davgfirst) call xysum_mn_name_z(heat,idiag_heatmz)
!
    endsubroutine calc_heat_cool
!***********************************************************************
    subroutine apply_floor(a)
!
      real, dimension(nx), intent(inout) :: a
      integer :: i
!
      do i=1,nx
        if (a(i) < TT_floor) a(i)=TT_floor
      enddo
!
    endsubroutine apply_floor
!***********************************************************************
    subroutine calc_heat_cool_variable(heat,p)
!
! Thermal relaxation for radially stratified global Keplerian disks
!
      use EquationOfState, only: rho0,gamma1
      use Sub, only: quintic_step
!
      real, dimension(nx), intent(inout) :: heat
      real, dimension (nx) :: rr1,TT_drive, OO, pborder
      real, save :: tau1_cool,rho01
      logical, save :: lfirstcall=.true.
      type (pencil_case), intent(in) :: p
!
      if (lcartesian_coords.or.lcylindrical_coords) then
        rr1=p%rcyl_mn1
      elseif (lspherical_coords) then
        rr1=p%r_mn1
      endif
      !period=2*pi/omega
      !period_inv=.5*pi_1*rr1**1.5
      OO=rr1**1.5
!
! Set the constants
!
      if (lfirstcall) then
        tau1_cool=1./tau_cool
        if (lphotoelectric_heating) rho01=1./rho0
        lfirstcall=.false.
      endif
!
      if (lphotoelectric_heating) then
        TT_drive=TTref_cool*rr1**TT_powerlaw*p%rhop*rho01
      else
        TT_drive=TTref_cool*rr1**TT_powerlaw
      endif
!
      if (TT_floor /= impossible) call apply_floor(TT_drive)
!
!  Apply tapering function. The idea is to smooth it to zero near the
!  boundaries. Useful in connection with border profiles when the initial
!  condition is not the same as the power law.
!
      if (lborder_heat_variable) then
        pborder=quintic_step(x(l1:l2),r_int,widthss,SHIFT=1.) - &
                quintic_step(x(l1:l2),r_ext,widthss,SHIFT=-1.)
      else
        pborder=1.
      endif
!
!  tau_cool = tau_cool_0 * Omega0/Omega
!  tau_cool1 = 1/tau_cool_0 * Omega/Omega0
!
      heat=heat-p%rho*p%cp*gamma1*&
           (p%TT-TT_drive)*tau1_cool*OO*pborder
!
    endsubroutine calc_heat_cool_variable
!***********************************************************************
    subroutine get_heat_cool_general(heat,p)
!
!  Subroutine to do volume heating and cooling in a layer independent of
!  gravity.
!
!  17-may-10/dhruba: coded
!  16-nov-13/nishant: two-layer (jumps) at arb heights and arb temperatures
!  31-dec-14/axel: added cubic_step (as in corona of temperature_ionization)
!
      use Sub, only: erfunc, cubic_step
!
      real, dimension (nx) :: heat
      type (pencil_case) :: p
!
      real, dimension (nx) :: prof,prof2
      real :: ztop, zbot
!
      intent(in) :: p
!
!  Calculate profile.
!
      select case (cooling_profile)
      case ('gaussian-z')
        prof=spread(exp(-0.5*((zcool-z(n))/wcool)**2), 1, l2-l1+1)
!
!  Cooling with a profile linear in z (unstable).
!
      case ('lin-z')
        prof=spread(z(n)/wcool,1,l2-l1+1)
!
!  Sinusoidal cooling profile (periodic).
!
      case ('sin-z')
        prof=spread(sin(z(n)/wcool),1,l2-l1+1)
!
!  Error function cooling profile.
!
      case ('surface_z')
        prof=spread(.5*(1.+erfunc((z(n)-zcool)/wcool)),1,l2-l1+1)
!
!  Error function cooling profile (two-layer).
!
      case ('two-layer')
        prof=spread(.5*(1.+erfunc((z(n)-zcool)/wcool)),1,l2-l1+1)
        prof2=spread(.5*(1.+erfunc((z(n)-zcool2)/wcool2)),1,l2-l1+1)
!
!  add "coronal" heating (to simulate a hot corona; see temperature_ionization)
!
      case ('cubic_step')
        ztop=xyz0(3)+Lxyz(3)
        prof = cubic_step(z(n),(ztop+z_cor)/2,(ztop-z_cor)/2)
!
!  Similar to cubic_step, but with the same type of step at negative z.
!  This is useful for accretion discs.
!
      case ('cubic_step_topbot')
        zbot=xyz0(3)
        ztop=xyz0(3)+Lxyz(3)
        prof = cubic_step(z(n),(ztop+z_cor)/2,(ztop-z_cor)/2) &
              +cubic_step(z(n),(zbot-z_cor)/2,(zbot+z_cor)/2)
!
!  Error function cooling profile with respect to pressure.
!
      case ('surface_pp')
        prof=.5*(1.-erfunc((p%pp-ppcool)/wcool))
!
      endselect
!
!  Note: the cooltype 'Temp' used below was introduced by Axel for the
!  aerosol runs. Although this 'Temp' does not match with the cooltype
!  'Temp' used in other parts of this subroutine. I have introduced
!  'Temp2' which is the same as the 'Temp' elsewhere. Later runs
!  will clarify this. - Dhruba
!  AB: not sure; in general we need to multiply with cv, as in 'corona'
!
     select case (cooltype)
     case('constant')
       heat=heat-cool*prof
     case('corona')
       heat=heat-cool*prof*p%cv*p%rho*(p%TT-TT_cor)
     case ('Temp')
       if (headtt) print*, 'calc_heat_cool: cs20,cs2cool=', cs20, cs2cool
       heat=heat-cool*(p%cs2-(cs20-prof*cs2cool))/cs2cool
     case('Temp2')
       heat=heat-cool*prof*(p%cs2-cs2cool)/cs2cool
     case('rho_cs2')
       heat=heat-cool*prof*p%rho*(p%cs2-cs2cool)
     case ('two-layer')
       heat = heat - cool *prof *p%rho*(p%cs2-cs2cool) &
                   - cool2*prof2*p%rho*(p%cs2-cs2cool2)
     case('plain')
       heat=heat-cool*prof
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
      use HDF5_IO, only: output_profile
      use Sub, only: step, cubic_step
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
      if (luminosity/=0.0) then
        if (nygrid==1) then
          prof = spread(exp(-0.5*((z(n)-zbot)/wheat)**2), 1, l2-l1+1) &
               /(sqrt(pi/2.)*wheat*Lx)
        else
          prof = spread(exp(-0.5*((z(n)-zbot)/wheat)**2), 1, l2-l1+1) &
               /(sqrt(pi/2.)*wheat*Lx*Ly)
        endif
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
        prof = spread(step(z(n),z2,wcool),1,nx)
      case ('step2')
        prof = spread(step(z(n),zcool,wcool),1,nx)
      case ('cubic_step')
        prof = spread(cubic_step(z(n),z2,wcool),1,nx)
!
!  Cooling with a profile linear in z (unstable).
!
      case ('lin-z')
        prof=spread(z(n)/wcool, 1, l2-l1+1)
      case default
        call fatal_error('get_heat_cool_gravz','please select a cooltype')
      endselect
!
      if (lcalc_cs2mz_mean) then
        heat = heat - cool*prof*(cs2mz(n)-cs2cool)/cs2cool
      else 
        heat = heat - cool*prof*(p%cs2-cs2cool)/cs2cool
      endif
!
!  Write out cooling profile (during first time step only) and apply.
!  MR: Later to be moved to initialization!
!
      if (m==m1) call output_profile('cooling_profile',z(n:n),prof(1:1),'z', lsave_name=(n==n1))
!
!  Write divergence of cooling flux.
!
      if (l1davgfirst) call xysum_mn_name_z(heat,idiag_dcoolz)
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
      use Debug_IO, only: output_pencil
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
           call output_pencil('heat.dat',heat,1)
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
          theta_profile=(onethird-(p%rcyl_mn/p%z_mn)**2)*deltaT_poleq
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
      use Diagnostics, only: zsum_mn_name_xy, yzsum_mn_name_x
      use Debug_IO, only: output_pencil
      use Sub, only: step
!
      type (pencil_case) :: p
!
      real, dimension (nx) :: heat,prof,prof2,cs2_tmp,fac_tmp
      real, dimension (nx), save :: cs2cool_x
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
          call output_pencil('heat.dat',heat,1)
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
!  Heating/cooling with two cooling layers on top of each other
!
      case ('shell2')
        if (rcool==0.0) rcool=r_ext
        if (rcool2==0.0) rcool2=r_ext
        prof2= step(x(l1:l2),rcool2,wcool2)
        prof = step(x(l1:l2),rcool,wcool)-prof2
        heat = heat - cool*prof*(p%cs2-cs2cool)/cs2cool-cool2*prof2*(p%cs2-cs2cool2)/cs2cool2
!
!  Similar to shell2, it cools/heats to defined profile
!
      case ('shell3')
        if (rcool==0.0) rcool=r_ext
        if (rcool2==0.0) rcool2=r_ext
        prof = step(x(l1:l2),rcool,wcool)
        prof2 = cs2cool + abs(cs2cool2-cs2cool)*step(x(l1:l2),rcool2,wcool2)
        heat = heat - cool*prof*(p%cs2-prof2)/prof2
!
!  Cool the mean temperature toward a specified profile stored in a file
!
      case ('shell_mean_yz')
        if (it == 1) call read_cooling_profile_x(cs2cool_x)
        if (rcool==0.0) rcool=r_ext
        prof = step(x(l1:l2),rcool,wcool)
        if (lcalc_cs2mean) then
          heat = heat - cool*prof*(cs2mx-cs2cool_x)
        else
          heat = heat - cool*prof*(p%cs2-cs2cool_x)
        endif
!
!  Cool (and not heat) toward a specified profile stored in a file
!
      case ('shell_mean_yz2')
        if (it == 1) call read_cooling_profile_x(cs2cool_x)
        if (rcool==0.0) rcool=r_ext
        prof = step(x(l1:l2),rcool,wcool)
        cs2_tmp=p%cs2-cs2cool_x
        where (cs2_tmp < 0.) cs2_tmp=0.
        heat = heat - cool*prof*(p%cs2-cs2cool_x)
!
!  Cool only downflows toward a specified profile stored in a file
!
      case ('shell_mean_downflow')
        if (it == 1) call read_cooling_profile_x(cs2cool_x)
        if (rcool==0.0) rcool=r_ext
        prof = step(x(l1:l2),rcool,wcool)
        cs2_tmp=p%cs2-cs2cool_x
        where (p%uu(:,1) .gt. 0.)
          fac_tmp=0.
        elsewhere
          fac_tmp=1.
        endwhere
        heat = heat - cool*prof*fac_tmp*(p%cs2-downflow_cs2cool_fac*cs2cool_x)
!
!  Latitude dependent heating/cooling: imposes a latitudinal variation
!  of temperature proportional to cos(theta) at each depth. deltaT gives
!  the amplitude of the variation between theta_0 and the equator.
!
      case ('latheat')
        if (rcool==0.0) rcool=r_ext
        prof = step(x(l1:l2),rcool1,wcool)-step(x(l1:l2),rcool2,wcool)
        prof1= 1.+deltaT*cos(2.*pi*(y(m)-y0)/Ly)
        heat = heat - cool*prof*(cs2mxy(:,m)-prof1*cs2mx)
!
!  Latitude dependent heating/cooling (see above) plus additional cooling
!  layer on top.
!
      case ('shell+latheat')
        if (rcool==0.0) rcool=r_ext
        prof = step(x(l1:l2),rcool1,wcool)-step(x(l1:l2),rcool2,wcool)
        prof1= 1.+deltaT*cos(2.*pi*(y(m)-y0)/Ly)
        heat = heat - cool*prof*(cs2mxy(:,m)-prof1*cs2mx)
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
        heat = heat - cool*prof*(ssmxy(:,m)-prof1*ssmx(l1+nghost:l2+nghost))
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
      if (l1davgfirst) call yzsum_mn_name_x(heat,idiag_dcoolx)
      if (l2davgfirst) call zsum_mn_name_xy(heat,idiag_dcoolxy)
!
    endsubroutine get_heat_cool_gravx_spherical
!***********************************************************************
    subroutine get_heat_cool_gravx_cartesian(heat,p)
!
!  Subroutine to calculate the heat/cool term in cartesian coordinates
!  with gravity along x direction.
!
!  This is equivalent to the previous routine.
!
!   4-aug-10/gustavo: coded
!
      use Debug_IO, only: output_pencil
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
          call output_pencil('heat.dat',heat,1)
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
            'get_heat_cool_gravx_cartesian: No such value for cooltype: ', trim(cooltype)
        call fatal_error('calc_heat_cool',errormsg)
      endselect
!
    endsubroutine get_heat_cool_gravx_cartesian
!***********************************************************************
    subroutine get_heat_cool_corona(heat,p)
!
!  Subroutine to calculate the heat/cool term for hot corona
!  Assume a linearly increasing reference profile.
!  This 1/rho1 business is clumsy, but so would be obvious alternatives...
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
        if (lcalc_cs2mz_mean) then
          heat=heat+profile_cor*(TT_cor-cs2mz(n)/gamma_m1*p%cp1)/(p%rho1*tau_cor*p%cp1)
        else
          heat=heat+profile_cor*(TT_cor-1/p%TT1)/(p%rho1*tau_cor*p%cp1)
        endif
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
      use Debug_IO, only: output_pencil
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: lnQ,rtv_cool,lnTT_SI,lnneni,delta_lnTT,tmp
      integer :: i,imax
      real :: unit_lnQ
!
!  Parameters for subroutine cool_RTV in SI units (from Cook et al. 1989).
!
      double precision, parameter, dimension (10) :: &
          intlnT_1 =(/4.605, 8.959, 9.906, 10.534, 11.283, 12.434, 13.286, 14.541, 17.51, 20.723 /)
      double precision, parameter, dimension (9) :: &
          lnH_1 = (/ -542.398,  -228.833, -80.245, -101.314, -78.748, -53.88, -80.452, -70.758, -91.182/), &
       B_1   = (/     50.,      15.,      0.,      2.0,      0.,    -2.,      0., -twothird,    0.5 /)
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
      intent(inout) :: df
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
    select case (cool_type)
      case (1)
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
        if (Pres_cutoff/=impossible) then
          rtv_cool=exp(lnneni+lnQ-unit_lnQ-p%lnTT-p%lnrho)*exp(-p%rho*p%cs2*gamma1/Pres_cutoff)
        else
          rtv_cool=exp(lnneni+lnQ-unit_lnQ-p%lnTT-p%lnrho)
        endif
        tmp=maxval(rtv_cool*gamma)/(cdts)
      case (2)
!
!  Second set of parameters
!
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
        if (Pres_cutoff/=impossible) then
          rtv_cool=exp(lnneni+lnQ-unit_lnQ-p%lnTT-p%lnrho)*exp(-p%rho*p%cs2*gamma1/Pres_cutoff)
        else
          rtv_cool=exp(lnneni+lnQ-unit_lnQ-p%lnTT-p%lnrho)
        endif
        tmp=maxval(rtv_cool*gamma)/(cdts)
      case (3)
        rtv_cool=0.0
        if (z(n) > z_cor) then
          call get_lnQ(lnTT_SI, lnQ, delta_lnTT)
          rtv_cool = exp(lnQ-unit_lnQ+lnneni-p%lnTT-p%lnrho)
        endif
        tmp = max (rtv_cool/cdts, abs (rtv_cool/max (tini, delta_lnTT)))
      case default
        rtv_cool(:)=0.
        tmp=maxval(rtv_cool*gamma)/(cdts)
    end select
!
      rtv_cool=rtv_cool * cool_RTV  ! for adjusting by setting cool_RTV in run.in
!
!  Add to entropy equation.
!
      if (lfirst .and. ip == 13) &
           call output_pencil('rtv.dat',rtv_cool*exp(p%lnrho+p%lnTT),1)
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss)-rtv_cool
!
      if (lfirst.and.ldt) then
        dt1_max=max(dt1_max,tmp)
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
      intent(inout) :: df
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
!  12-nov-10/mvaisala: adapted from calc_heat_cool
!  12-feb-15/MR: adapted for reference state.
!
      use Debug_IO, only: output_pencil
      use EquationOfState, only: eoscalc,ilnrho_ss,irho_ss
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: p
      intent(inout) :: df
!
      real, dimension (nx) :: pp
!
!  Add heating/cooling to entropy equation.
!
      if (ldensity_nolog) then
        if (lreference_state) then
          call eoscalc(irho_ss,f(l1:l2,m,n,irho)+reference_state(:,iref_rho), &
                               f(l1:l2,m,n,iss)+reference_state(:,iref_s),pp=pp)
        else
          call eoscalc(irho_ss,f(l1:l2,m,n,irho),f(l1:l2,m,n,iss),pp=pp)
        endif
      else
        call eoscalc(ilnrho_ss,f(l1:l2,m,n,ilnrho),f(l1:l2,m,n,iss),pp=pp)
      endif
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + 0.5*pp*p%divu
!
!  0.5 Because according to the virial theorem a collapsing cloud loses approximately
!  half of its internal energy to radiation.
!
      !df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + p%pp*p%divu
!
    endsubroutine calc_heat_cool_prestellar
!***********************************************************************
    subroutine rprint_energy(lreset,lwrite)
!
!  Reads and registers print parameters relevant to entropy.
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Diagnostics, only: parse_name,set_type
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
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtc=0; idiag_ethm=0; idiag_ethdivum=0
        idiag_ssruzm=0; idiag_ssuzm=0; idiag_ssm=0; idiag_ss2m=0
        idiag_eem=0; idiag_ppm=0; idiag_csm=0; idiag_cgam=0; idiag_pdivum=0; idiag_heatm=0
        idiag_ugradpm=0; idiag_ethtot=0; idiag_dtchi=0; idiag_ssmphi=0
        idiag_fradbot=0; idiag_fradtop=0; idiag_TTtop=0
        idiag_yHmax=0; idiag_yHm=0; idiag_TTmax=0; idiag_TTmin=0; idiag_TTm=0
        idiag_ssmax=0; idiag_ssmin=0; idiag_gTmax=0; idiag_csmax=0
        idiag_gTrms=0; idiag_gsrms=0; idiag_gTxgsrms=0
        idiag_fconvm=0; idiag_fconvz=0; idiag_dcoolz=0; idiag_heatmz=0; idiag_fradz=0
        idiag_Fenthz=0; idiag_Fenthupz=0; idiag_Fenthdownz=0
        idiag_fturbz=0; idiag_ppmx=0; idiag_ppmy=0; idiag_ppmz=0
        idiag_fturbtz=0; idiag_fturbmz=0; idiag_fturbfz=0
        idiag_ssmx=0; idiag_ss2mx=0; idiag_ssmy=0; idiag_ssmz=0; idiag_ss2mz=0
        idiag_ssupmz=0; idiag_ssdownmz=0; idiag_ss2upmz=0; idiag_ss2downmz=0
        idiag_ssf2mz=0; idiag_ssf2upmz=0; idiag_ssf2downmz=0
        idiag_ssmr=0; idiag_TTmr=0
        idiag_TTmx=0; idiag_TTmy=0; idiag_TTmz=0; idiag_TTmxy=0; idiag_TTmxz=0
        idiag_TTupmz=0; idiag_TTdownmz=0
        idiag_uxTTmz=0; idiag_uyTTmz=0; idiag_uzTTmz=0; idiag_cs2mphi=0
        idiag_uzTTupmz=0; idiag_uzTTdownmz=0
        idiag_ssmxy=0; idiag_ssmxz=0; idiag_fradz_Kprof=0; idiag_uxTTmxy=0
        idiag_uyTTmxy=0; idiag_uzTTmxy=0; idiag_TT2mx=0; idiag_TT2mz=0
        idiag_TT2upmz=0; idiag_TT2downmz=0
        idiag_TTf2mz=0; idiag_TTf2upmz=0; idiag_TTf2downmz=0
        idiag_uxTTmx=0; idiag_uyTTmx=0; idiag_uzTTmx=0;
        idiag_fturbxy=0; idiag_fturbrxy=0; idiag_fturbthxy=0; idiag_fturbmx=0
        idiag_fradxy_Kprof=0; idiag_fconvxy=0; idiag_fradmx=0
        idiag_fradx_kramers=0; idiag_fradz_kramers=0; idiag_fradxy_kramers=0
        idiag_fconvyxy=0; idiag_fconvzxy=0; idiag_dcoolx=0; idiag_dcoolxy=0
        idiag_ufpresm=0; idiag_uduum=0; idiag_fradz_constchi=0
        idiag_gTxmxy=0; idiag_gTymxy=0; idiag_gTzmxy=0
        idiag_gsxmxy=0; idiag_gsymxy=0; idiag_gszmxy=0
        idiag_gTxgsxmxy=0;idiag_gTxgsymxy=0;idiag_gTxgszmxy=0
        idiag_fradymxy_Kprof=0; idiag_fturbymxy=0; idiag_fconvxmx=0
        idiag_Kkramersm=0; idiag_Kkramersmx=0; idiag_Kkramersmz=0
        idiag_Hmax=0; idiag_dtH=0; idiag_tauhmin=0; idiag_ethmz=0
      endif
!
!  iname runs through all possible names that may be listed in print.in.
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtc',idiag_dtc)
        call parse_name(iname,cname(iname),cform(iname),'dtchi',idiag_dtchi)
        call parse_name(iname,cname(iname),cform(iname),'dtH',idiag_dtH)
        call parse_name(iname,cname(iname),cform(iname),'Hmax',idiag_Hmax)
        call parse_name(iname,cname(iname),cform(iname),'tauhmin',idiag_tauhmin)
        call parse_name(iname,cname(iname),cform(iname),'ethtot',idiag_ethtot)
        call parse_name(iname,cname(iname),cform(iname),'ethdivum',idiag_ethdivum)
        call parse_name(iname,cname(iname),cform(iname),'ethm',idiag_ethm)
        call parse_name(iname,cname(iname),cform(iname),'ssruzm',idiag_ssruzm)
        call parse_name(iname,cname(iname),cform(iname),'ssuzm',idiag_ssuzm)
        call parse_name(iname,cname(iname),cform(iname),'ssm',idiag_ssm)
        call parse_name(iname,cname(iname),cform(iname),'ss2m',idiag_ss2m)
        call parse_name(iname,cname(iname),cform(iname),'eem',idiag_eem)
        call parse_name(iname,cname(iname),cform(iname),'ppm',idiag_ppm)
        call parse_name(iname,cname(iname),cform(iname),'pdivum',idiag_pdivum)
        call parse_name(iname,cname(iname),cform(iname),'heatm',idiag_heatm)
        call parse_name(iname,cname(iname),cform(iname),'csm',idiag_csm)
        call parse_name(iname,cname(iname),cform(iname),'csmax',idiag_csmax)
        call parse_name(iname,cname(iname),cform(iname),'cgam',idiag_cgam)
        call parse_name(iname,cname(iname),cform(iname),'ugradpm',idiag_ugradpm)
        call parse_name(iname,cname(iname),cform(iname),'fradbot',idiag_fradbot)
        call parse_name(iname,cname(iname),cform(iname),'fradtop',idiag_fradtop)
        call parse_name(iname,cname(iname),cform(iname),'TTtop',idiag_TTtop)
        call parse_name(iname,cname(iname),cform(iname),'yHm',idiag_yHm)
        call parse_name(iname,cname(iname),cform(iname),'yHmax',idiag_yHmax)
        call parse_name(iname,cname(iname),cform(iname),'TTm',idiag_TTm)
        call parse_name(iname,cname(iname),cform(iname),'TTmax',idiag_TTmax)
        call parse_name(iname,cname(iname),cform(iname),'TTmin',idiag_TTmin)
        call parse_name(iname,cname(iname),cform(iname),'gTmax',idiag_gTmax)
        call parse_name(iname,cname(iname),cform(iname),'ssmin',idiag_ssmin)
        call parse_name(iname,cname(iname),cform(iname),'ssmax',idiag_ssmax)
        call parse_name(iname,cname(iname),cform(iname),'gTrms',idiag_gTrms)
        call parse_name(iname,cname(iname),cform(iname),'gsrms',idiag_gsrms)
        call parse_name(iname,cname(iname),cform(iname),'gTxgsrms',idiag_gTxgsrms)
        call parse_name(iname,cname(iname),cform(iname),'TTp',idiag_TTp)
        call parse_name(iname,cname(iname),cform(iname),'fconvm',idiag_fconvm)
        call parse_name(iname,cname(iname),cform(iname),'ufpresm',idiag_ufpresm)
        call parse_name(iname,cname(iname),cform(iname),'uduum',idiag_uduum)
        call parse_name(iname,cname(iname),cform(iname),'Kkramersm',idiag_Kkramersm)
      enddo
!
      if (idiag_fradbot/=0) call set_type(idiag_fradbot,lsurf=.true.)
      if (idiag_fradtop/=0) call set_type(idiag_fradtop,lsurf=.true.)
      if (idiag_TTtop/=0) call set_type(idiag_TTtop,lsurf=.true.)
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ssmx',idiag_ssmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ss2mx',idiag_ss2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ppmx',idiag_ppmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'TTmx',idiag_TTmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'TT2mx',idiag_TT2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'fconvxmx',idiag_fconvxmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'uxTTmx',idiag_uxTTmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'uyTTmx',idiag_uyTTmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'uzTTmx',idiag_uzTTmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'fradmx',idiag_fradmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'fturbmx',idiag_fturbmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'Kkramersmx',idiag_Kkramersmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'dcoolx',idiag_dcoolx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'fradx_kramers',idiag_fradx_kramers)
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
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fturbtz',idiag_fturbtz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fturbmz',idiag_fturbmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fturbfz',idiag_fturbfz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fconvz',idiag_fconvz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'Fenthz',idiag_Fenthz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'Fenthupz',idiag_Fenthupz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'Fenthdownz',idiag_Fenthdownz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'dcoolz',idiag_dcoolz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'heatmz',idiag_heatmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fradz',idiag_fradz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fradz_Kprof',idiag_fradz_Kprof)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fradz_constchi',idiag_fradz_constchi)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fradz_kramers',idiag_fradz_kramers)
        call parse_name(inamex,cnamez(inamez),cformz(inamez),'Kkramersmz',idiag_Kkramersmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ssmz',idiag_ssmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ssupmz',idiag_ssupmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ssdownmz',idiag_ssdownmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ss2mz',idiag_ss2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ss2upmz',idiag_ss2upmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ss2downmz',idiag_ss2downmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ssf2mz',idiag_ssf2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ssf2upmz',idiag_ssf2upmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ssf2downmz',idiag_ssf2downmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TTmz',idiag_TTmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TTupmz',idiag_TTupmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TTdownmz',idiag_TTdownmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TT2mz',idiag_TT2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TT2upmz',idiag_TT2upmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TT2downmz',idiag_TT2downmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TTf2mz',idiag_TTf2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TTf2upmz',idiag_TTf2upmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TTf2downmz',idiag_TTf2downmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ppmz',idiag_ppmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uxTTmz',idiag_uxTTmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uyTTmz',idiag_uyTTmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzTTmz',idiag_uzTTmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzTTupmz',idiag_uzTTupmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzTTdownmz',idiag_uzTTdownmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'gTxgsxmz',idiag_gTxgsxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'gTxgsymz',idiag_gTxgsymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'gTxgszmz',idiag_gTxgszmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'gTxgsx2mz',idiag_gTxgsx2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'gTxgsy2mz',idiag_gTxgsy2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'gTxgsz2mz',idiag_gTxgsz2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ethmz',idiag_ethmz)
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
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'fturbymxy',idiag_fturbymxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'fturbrxy',idiag_fturbrxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'fturbthxy',idiag_fturbthxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'fradxy_Kprof',idiag_fradxy_Kprof)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'fradymxy_Kprof',idiag_fradymxy_Kprof)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'fradxy_kramers',idiag_fradxy_kramers)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'fconvxy',idiag_fconvxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'fconvyxy',idiag_fconvyxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'fconvzxy',idiag_fconvzxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'dcoolxy',idiag_dcoolxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'gTxmxy',idiag_gTxmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'gTymxy',idiag_gTymxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'gTzmxy',idiag_gTzmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'gsxmxy',idiag_gsxmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'gsymxy',idiag_gsymxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'gszmxy',idiag_gszmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'gTxgsxmxy',idiag_gTxgsxmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'gTxgsymxy',idiag_gTxgsymxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'gTxgszmxy',idiag_gTxgszmxy)
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
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then 
        where(cnamev=='pp'.or.cnamev=='ss'.or.cnamev=='TT'.or.cnamev=='lnTT') cformv='DEFINED'
      endif
!
!  Write column where which entropy variable is stored.
!
      if (lwr) then
        call farray_index_append('iyH',iyH)
        call farray_index_append('ilnTT',ilnTT)
        call farray_index_append('iTT',iTT)
      endif
!
    endsubroutine rprint_energy
!***********************************************************************
    subroutine get_slices_energy(f,slices)
!
!  Write slices for animation of Entropy variables.
!
!  26-jul-06/tony: coded
!  12-feb-15/MR  : changes for use of reference state.
!
      use EquationOfState, only: eoscalc, ilnrho_ss, irho_ss
      use Slices_methods, only: assign_slices_scal, addto_slices, process_slices, exp2d
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      character(LEN=labellen) :: sname
      real :: ssadd
      integer :: l,n,m
!
      ssadd=0.
      sname=trim(slices%name)
!
!  Loop over slices
!
      select case (sname)
!
!  Entropy.
!
        case ('ss')
          call assign_slices_scal(slices,f,iss)
          if (lfullvar_in_slices) call addto_slices(slices,reference_state(:,iref_s))
!
!  Pressure.
!
        case ('pp')
!
          if (lwrite_slice_yz) then
            if (lreference_state) ssadd=reference_state(ix_loc-l1+1,iref_s)
            do m=m1,m2; do n=n1,n2
              call eoscalc(irho_ss,getrho_s(f(ix_loc,m,n,ilnrho),ix_loc), &
                           f(ix_loc,m,n,iss)+ssadd,pp=slices%yz(m-m1+1,n-n1+1))
            enddo; enddo
          endif

          do l=l1,l2 
!
            if (lreference_state) ssadd=reference_state(l-l1+1,iref_s)

            do n=n1,n2
              if (lwrite_slice_xz) &
                call eoscalc(irho_ss,getrho_s(f(l,iy_loc,n,ilnrho),l), &
                             f(l,iy_loc,n,iss)+ssadd,pp=slices%xz(l-l1+1,n-n1+1))
              if (lwrite_slice_xz2) &
                call eoscalc(irho_ss,getrho_s(f(l,iy2_loc,n,ilnrho),l), &
                             f(l,iy2_loc,n,iss)+ssadd,pp=slices%xz2(l-l1+1,n-n1+1))
            enddo

            do m=m1,m2
              if (lwrite_slice_xy) &
                call eoscalc(irho_ss,getrho_s(f(l,m,iz_loc,ilnrho),l), &
                             f(l,m,iz_loc, iss)+ssadd,pp=slices%xy(l-l1+1,m-m1+1))
              if (lwrite_slice_xy2) &
                call eoscalc(irho_ss,getrho_s(f(l,m,iz2_loc,ilnrho),l), &
                             f(l,m,iz2_loc,iss)+ssadd,pp=slices%xy2(l-l1+1,m-m1+1))
              if (lwrite_slice_xy3) &
                call eoscalc(irho_ss,getrho_s(f(l,m,iz3_loc,ilnrho),l), &
                             f(l,m,iz3_loc,iss)+ssadd,pp=slices%xy3(l-l1+1,m-m1+1))
              if (lwrite_slice_xy4) &
                call eoscalc(irho_ss,getrho_s(f(l,m,iz4_loc,ilnrho),l), &
                             f(l,m,iz4_loc,iss)+ssadd,pp=slices%xy4(l-l1+1,m-m1+1))
            enddo 
          enddo
          slices%ready=.true.
!
! Temperature
!
        case ('TT','lnTT')
!
          if (lwrite_slice_yz) then
            if (lreference_state) ssadd=reference_state(ix_loc-l1+1,iref_s)
            do m=m1,m2; do n=n1,n2
              call eoscalc(irho_ss,getrho_s(f(ix_loc,m,n,ilnrho),ix_loc), &
                           f(ix_loc,m,n,iss)+ssadd,lnTT=slices%yz(m-m1+1,n-n1+1))
            enddo; enddo
          endif

          do l=l1,l2

            if (lreference_state) ssadd=reference_state(l-l1+1,iref_s)

            do n=n1,n2
              if (lwrite_slice_xz) &
                call eoscalc(irho_ss,getrho_s(f(l,iy_loc,n,ilnrho),l), &
                             f(l,iy_loc,n,iss)+ssadd,lnTT=slices%xz(l-l1+1,n-n1+1))
              if (lwrite_slice_xz2) &
                call eoscalc(irho_ss,getrho_s(f(l,iy2_loc,n,ilnrho),l), &
                             f(l,iy2_loc,n,iss)+ssadd,lnTT=slices%xz2(l-l1+1,n-n1+1))
            enddo

            do m=m1,m2
              if (lwrite_slice_xy) &
                call eoscalc(irho_ss,getrho_s(f(l,m,iz_loc,ilnrho),l), &
                             f(l,m,iz_loc,iss)+ssadd,lnTT=slices%xy(l-l1+1,m-m1+1))
              if (lwrite_slice_xy2) &
                call eoscalc(irho_ss,getrho_s(f(l,m,iz2_loc,ilnrho),l), &
                             f(l,m,iz2_loc,iss)+ssadd,lnTT=slices%xy2(l-l1+1,m-m1+1))
              if (lwrite_slice_xy3) &
                call eoscalc(irho_ss,getrho_s(f(l,m,iz3_loc,ilnrho),l), &
                             f(l,m,iz3_loc,iss)+ssadd,lnTT=slices%xy3(l-l1+1,m-m1+1))
              if (lwrite_slice_xy4) &
                call eoscalc(irho_ss,getrho_s(f(l,m,iz4_loc,ilnrho),l), &
                             f(l,m,iz4_loc,iss)+ssadd,lnTT=slices%xy4(l-l1+1,m-m1+1))
            enddo
          enddo
!
          if (sname=='TT') call process_slices(slices,exp2d)
          slices%ready=.true.
!
      endselect
!
    endsubroutine get_slices_energy
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
        zprof_hcond(n) = 1. + (hcond1-1.)*cubic_step(zpt,z1,-widthss) &
                            + (hcond2-1.)*cubic_step(zpt,z2,+widthss)
        zprof_hcond(n) = hcond0*zprof_hcond(n)
        zprof_glhc(n,1:2) = 0.
        zprof_glhc(n,3) = (hcond1-1.)*cubic_der_step(zpt,z1,-widthss) &
                        + (hcond2-1.)*cubic_der_step(zpt,z2,+widthss)
        zprof_glhc(n,3) = hcond0*zprof_glhc(n,3)
      enddo
!
    endsubroutine calc_heatcond_zprof
!***********************************************************************
    subroutine get_gravz_heatcond
!
! Calculates z dependent heat conductivity and its logarithmic gradient.
! Scaled with hcond0=Kbot.
!
!  25-feb-2011/gustavo+dhruba: stolen from heatcond below
!
      use Gravity, only: z1, z2
      use Sub, only: step,der_step
!
      if (.not.lmultilayer) call fatal_error('get_gravz_heatcond:',&
           "don't call if you have only one layer")
!
      allocate(hcond_prof(nz),dlnhcond_prof(nz))
!
      hcond_prof = 1. + (hcond1-1.)*step(z(n1:n2),z1,-widthss) &
                      + (hcond2-1.)*step(z(n1:n2),z2, widthss)
!
      dlnhcond_prof = (  (hcond1-1.)*der_step(z(n1:n2),z1,-widthss) &
                       + (hcond2-1.)*der_step(z(n1:n2),z2, widthss))&
                      /hcond_prof

      hcond_prof = hcond0*hcond_prof
!
    endsubroutine get_gravz_heatcond
!***********************************************************************
    subroutine get_gravx_heatcond
!
! Calculates x dependent heat conductivity and its logarithmic gradient.
! No scaling with hcond0!
!
!  10-march-2011/gustavo: Adapted from idl script strat_poly.pro
!
      use SharedVariables, only: get_shared_variable
      use Sub, only: step, der_step
      use Mpicomm, only: stop_it

      real, pointer :: cv, gravx        ! ,xc,xb   ! are not put shared!
      real :: xb=.6, xc=.8              ! to set anything unproblematic
      real, dimension (:), pointer :: gravx_xpencil
      real, dimension (nx) :: dTTdxc,mpoly_xprof,dmpoly_dx
      real :: Lum
!
      if (.not.lmultilayer) call fatal_error('get_gravx_heatcond:',&
           "don't call if you have only one layer")

      !!call get_shared_variable('xb',xb, caller='get_gravx_heatcond')
      !!call get_shared_variable('xc',xc)

      call get_shared_variable('cv', cv, caller='get_gravx_heatcond')
      call get_shared_variable('gravx_xpencil', gravx_xpencil)
      call get_shared_variable('gravx', gravx)
!
! Radial profile of the polytropic index
      mpoly_xprof = mpoly0 - (mpoly0-mpoly1)*step(x(l1:l2),xb,widthss) &
                           - (mpoly1-mpoly2)*step(x(l1:l2),xc,widthss)
      dmpoly_dx = -(mpoly0-mpoly1)*der_step(x(l1:l2),xb,widthss) &
                  -(mpoly1-mpoly2)*der_step(x(l1:l2),xc,widthss)
!
! Hydrostatic equilibrium relations
      dTTdxc = gravx_xpencil / (cv*gamma_m1*(mpoly_xprof+1.))
      Lum = Fbot * (4.*pi*xyz0(1)**2)
!
! Kappa and its gradient are computed here
!
      allocate(hcond_prof(nx),dlnhcond_prof(nx))
      hcond_prof = -Lum/(4.*pi*x(l1:l2)**2*dTTdxc)
      dlnhcond_prof = Lum*cv*gamma_m1/(4.*pi*gravx) * dmpoly_dx/hcond_prof
!                     MR: how can this be the derivative of hcond_prof?
      if (lroot) print*, &
           ' get_gravx_heatcond: hcond computed from gravity dependent ' //&
           ' hydrostatic equilibrium relations'
!
    endsubroutine get_gravx_heatcond
!***********************************************************************
    subroutine chit_profile
!
!  Calculate chit and its derivative where chit is the turbulent heat diffusivity.
!  Scaling with chi_t or chhi_t0 is not performed!
!
!  23-jan-2002/wolf: coded
!
      use Gravity, only: z1, z2
      use Sub, only: step,der_step
!
      real :: zbot, ztop
!
      if (lgravz) then
!
!  If zz1 and/or zz2 are not set, use z1 and z2 instead.
!
        if (.not.lmultilayer) call fatal_error('chit_profile', &
             "don't call for lgravz=T if you have only one layer")

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

        if (.not.allocated(chit_prof_stored)) &
          allocate(chit_prof_stored(nz),dchit_prof_stored(nz))

        chit_prof_stored = 1. + (chit_prof1-1.)*step(z(n1:n2),zbot,-widthss) &
                              + (chit_prof2-1.)*step(z(n1:n2),ztop, widthss)
        dchit_prof_stored =  (chit_prof1-1.)*der_step(z(n1:n2),zbot,-widthss) &
                           + (chit_prof2-1.)*der_step(z(n1:n2),ztop, widthss)
!
      elseif (lspherical_coords.or.lconvection_gravx) then

        if (.not.allocated(chit_prof_stored)) &
          allocate(chit_prof_stored(nx),dchit_prof_stored(nx))

        select case (ichit)
          case ('nothing')
            chit_prof_stored = 1. + (chit_prof1-1.)*step(x(l1:l2),xbot,-widthss) &
                                  + (chit_prof2-1.)*step(x(l1:l2),xtop, widthss)
            dchit_prof_stored =  (chit_prof1-1.)*der_step(x(l1:l2),xbot,-widthss) &
                               + (chit_prof2-1.)*der_step(x(l1:l2),xtop, widthss)
          case ('powerlaw','power-law')
            chit_prof_stored = (x(l1:l2)/xchit)**(-pclaw)
            !dchit_prof_stored = -chi_t*pclaw*(x(l1:l2)/xchit)**(-pclaw-1.)/xchit
            dchit_prof_stored = -chit_prof_stored*pclaw/x(l1:l2)

          case default
            call warning('chit_profile','no profiles calculated')
        endselect

      elseif (.not.lsphere_in_a_box) then
        call warning('chit_profile','no profiles calculated')
      endif
!
    endsubroutine chit_profile
!***********************************************************************
    subroutine get_prof_pencil(prof,dprof,l2D3D,amp,amp1,amp2,pos1,pos2,p,f,stored_prof,stored_dprof,llog)
!
!  Provides pencils for an either x or z dependent profile and its (optionally logarithmic) gradient:
!  prof and dprof. If 
!  lmultilayer=T there are real profiles, otherwise prof=amp is constant and dprof=0.
!  If the pencil case is provided (present(p)) and lhcond_global=T, these data are fetched from
!  the f-array, otherwise, for l2D3D=F fetched from stored_prof and stored_dprof.
!  For l2D3D=T they are calculated from amp,amp1,amp2,pos1,pos2.
!  In the latter case, the radial coordinate is calculated from x,y,z or taken from the pencil case,
!  if present.
!  If llog=T the logarithmic gradient is provided, otherwise the normal one.
! 
!  20-jan-20/MR: coded
!
      use Sub, only: step, der_step
      use General, only: loptest
!
      real, dimension(nx)  , intent(out) :: prof
      real, dimension(nx,3), intent(out) :: dprof
      logical,               intent(in ) :: l2D3D
      real, optional,        intent(in ) :: amp,amp1,amp2,pos1,pos2 
      logical, optional,     intent(in ) :: llog
      type(pencil_case),                 optional, intent(in) :: p
      real, dimension(mx,my,mz,mfarray), optional, intent(in) :: f
      real, dimension(:),                optional, intent(in) :: stored_prof, stored_dprof
!
      real, dimension(nx) :: r_mn, r_mn1
      integer :: j

      if (.not.lmultilayer) then
        prof=amp; dprof=0.
        return
      endif

      if (lgravz) then
        prof=stored_prof(n-nghost)
        dprof(:,3)=stored_dprof(n-nghost); dprof(:,1:2)=0.
      elseif (l2D3D) then

        if (present(p).and.lhcond_global) then
          prof = f(l1:l2,m,n,iglobal_hcond)
          dprof= f(l1:l2,m,n,iglobal_glhc:iglobal_glhc+2)
        else

          if (present(p)) then
            r_mn=p%r_mn
            r_mn1=p%r_mn1
          else
            r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
            r_mn1=1./r_mn
          endif

          prof =  (amp1-1.)*der_step(r_mn,pos1,-widthss) &
                 +(amp2-1.)*der_step(r_mn,pos2, widthss)
          dprof(:,1) = prof*x(l1:l2)*r_mn1
          dprof(:,2) = prof*y(  m  )*r_mn1

          if (lcylinder_in_a_box) then
            dprof(:,3) = 0.0
          else
            dprof(:,3) = prof*z(n)*r_mn1
          endif

          prof = 1.+(amp1-1.)*step(r_mn,pos1,-widthss) &
                   +(amp2-1.)*step(r_mn,pos2, widthss)

          if (loptest(llog)) then
            do j=1,3; dprof(:,j)=dprof(:,j)/prof; enddo
          else
            dprof = amp*dprof
          endif

          prof = amp*prof

        endif

      else  ! covers also lgravr=T 
        prof=stored_prof
        dprof(:,1)=stored_dprof; dprof(:,2:3)=0.
      endif

    endsubroutine get_prof_pencil
!***********************************************************************
    subroutine chit_profile_fluct
!
!  21-apr-2017/pete: aped from chit_profile
!
!  Calculates the chit_fluct profile and its derivative.
!  Scaled with chi_t1.
!
      use Gravity, only: z1, z2
      use Sub, only: step,der_step
      use Messages, only: warning
!
      real :: zbot, ztop
!
      if (lgravz) then
!
!  If zz1_fluct and/or zz2_fluct are not set, use z1 and z2 instead.
!
        if (zz1_fluct == impossible) then
          zbot=z1
        else
          zbot=zz1_fluct
        endif
!
        if (zz2_fluct == impossible) then
          ztop=z2
        else
          ztop=zz2_fluct
        endif
!
        if (.not.allocated(chit_prof_fluct_stored)) &
          allocate(chit_prof_fluct_stored(nz),dchit_prof_fluct_stored(nz))
        chit_prof_fluct_stored = chi_t1*(1. + (chit_fluct_prof1-1.)*step(z(n1:n2),zbot,-widthss) &
                                            + (chit_fluct_prof2-1.)*step(z(n1:n2),ztop, widthss))
        dchit_prof_fluct_stored = chi_t1*(  (chit_fluct_prof1-1.)*der_step(z(n1:n2),zbot,-widthss) &
                                          + (chit_fluct_prof2-1.)*der_step(z(n1:n2),ztop, widthss))
      elseif (lgravx) then
        if (.not.allocated(chit_prof_fluct_stored)) &
          allocate(chit_prof_fluct_stored(nx),dchit_prof_fluct_stored(nx))
        chit_prof_fluct_stored = chi_t1*(1. + (chit_fluct_prof1-1.)*step(x(l1:l2),xbot_chit1,-widthss) &
                                            + (chit_fluct_prof2-1.)*step(x(l1:l2),xtop_chit1, widthss))
        dchit_prof_fluct_stored = chi_t1*(  (chit_fluct_prof1-1.)*der_step(x(l1:l2),xbot_chit1,-widthss) &
                                          + (chit_fluct_prof2-1.)*der_step(x(l1:l2),xtop_chit1, widthss))
      elseif (lgravr) then
        if (.not.allocated(chit_prof_fluct_stored)) &
          allocate(chit_prof_fluct_stored(nx),dchit_prof_fluct_stored(nx))
        chit_prof_fluct_stored = chi_t1
        dchit_prof_fluct_stored = 0.
      else
        call warning('chit_profile_fluct','no profile calculated')
      endif
!
    endsubroutine chit_profile_fluct
!***********************************************************************
    subroutine chit_aniso_profile
!
!  Calculates the chit_aniso profile and its derivative.
!  Scaled with chi_t*chit_aniso.
!  
!  27-apr-2011/pete: adapted from gradlogchit_profile
!
      use Sub, only: step, der_step
      use Messages, only: warning
!
      if (lspherical_coords) then

        chit_aniso_prof = chi_t*chit_aniso* &
            (1. + (chit_aniso_prof1-1.)*step(x(l1:l2),xbot_aniso,-widthss) &
                + (chit_aniso_prof2-1.)*step(x(l1:l2),xtop_aniso, widthss))

        dchit_aniso_prof = chi_t*chit_aniso* &
              (  (chit_aniso_prof1-1.)*der_step(x(l1:l2),xbot_aniso,-widthss) &
               + (chit_aniso_prof2-1.)*der_step(x(l1:l2),xtop_aniso, widthss))
      else
        call warning('chit_aniso_profile','no profile calculated')
      endif
!
    endsubroutine chit_aniso_profile
!***********************************************************************
    subroutine newton_cool(df,p)
!
!  Keeps the temperature in the lower chromosphere
!  at a constant level using newton cooling.
!
!  15-dec-2004/bing: coded
!  25-sep-2006/bing: updated, using external data
!
      use Debug_IO, only: output_pencil
      use EquationOfState, only: lnrho0, gamma
      use Mpicomm, only: mpibcast_real
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: newton
      integer, parameter :: prof_nz=150
      real, dimension (prof_nz), save :: prof_lnT,prof_z
      real :: lnTTor,dummy=1.
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
          inquire(IOLENGTH=lend) dummy
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
            call output_pencil('newton.dat',newton*exp(p%lnrho+p%lnTT),1)
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
    subroutine read_cooling_profile_x(cs2cool_x)
!
!  Read sound speed profile from an ascii-file
!
!  11-dec-2014/pete: aped from read_hcond
!
      use Mpicomm, only: mpibcast_real_arr, MPI_COMM_WORLD
!
      real, dimension(nx), intent(out) :: cs2cool_x
      real, dimension(nxgrid) :: tmp1
      real :: var1
      logical :: exist
      integer :: stat, nn
!
!  Read cs2_cool_x and write into an array.
!  If file is not found in run directory, search under trim(directory).
!
      if (lroot ) then
        inquire(file='cooling_profile.dat',exist=exist)
        if (exist) then
          open(36,file='cooling_profile.dat')
        else
          inquire(file=trim(directory)//'/cooling_profile.ascii',exist=exist)
          if (exist) then
            open(36,file=trim(directory)//'/cooling_profile.ascii')
          else
            call fatal_error('read_cooling_profile_x','no input file "cooling_profile.*"')
          endif
        endif
!
!  Read profiles.
!
        do nn=1,nxgrid
          read(36,*,iostat=stat) var1
          if (stat<0) exit
          if (ip<5) print*,'cs2cool_x: ',var1
          tmp1(nn)=var1
        enddo
        close(36)
!
      endif
!
      call mpibcast_real_arr(tmp1, nxgrid, comm=MPI_COMM_WORLD)
!
!  Assuming no ghost zones in cooling_profile.dat
!
      do nn=l1,l2
        cs2cool_x(nn-nghost)=tmp1(ipx*nx+nn-nghost)
      enddo
!
    endsubroutine read_cooling_profile_x
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
!        if (zm<=z1) then
        if (zm<z1) then
!
!  radiative zone=polytropic stratification
!
          del=1./(mpoly1+1.)
        else
          if (zm<=z2) then
!
!  convective zone=mixing-length stratification
!
            del=delad+alpha_MLT*(fc/ &
                (exp(lnrhom(iz-1))*(gamma_m1*tempm(iz-1))**1.5))**twothird
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
      use EquationOfState, only: eoscalc, ilnrho_lnTT, lnrho0, get_cp1
      use Gravity, only: g0
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: lnrho,lnTT,TT,ss,r_mn
      real :: beta0,beta1,TT_bcz,cp1
      real :: lnrho_int,lnrho_ext,lnrho_bcz
!
      if (headtt) print*,'r_bcz in entropy.f90=',r_bcz
!
!  The temperature gradient is dT/dr=beta/r with
!  beta = (g/cp) /[(1-1/gamma)*(m+1)]
!
      call get_cp1(cp1)
      beta0=cp1*g0/(mpoly0+1)*gamma/gamma_m1
      beta1=cp1*g0/(mpoly1+1)*gamma/gamma_m1
      TT_bcz=TT_ext+beta0*(1./r_bcz-1./r_ext)
      lnrho_ext=lnrho0
      lnrho_bcz=lnrho0+ &
            mpoly0*log(TT_ext+beta0*(1./r_bcz-1./r_ext))-mpoly0*log(TT_ext)
      lnrho_int=lnrho_bcz + &
            mpoly1*log(TT_bcz+beta1*(1./r_int-1./r_bcz))-mpoly1*log(TT_bcz)
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
          TT = TT_ext+beta0*(1./r_mn-1./r_ext)
          f(l1:l2,m,n,ilnrho)=lnrho0+mpoly0*log(TT/TT_ext)
        endwhere
!
!  Radiative layer.
!
        where (r_mn <= r_bcz .AND. r_mn > r_int)
          TT = TT_bcz+beta1*(1./r_mn-1./r_bcz)
          f(l1:l2,m,n,ilnrho)=lnrho_bcz+mpoly1*log(TT/TT_bcz)
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
      if (ldensity_nolog) then 
        f(:,:,:,irho) = exp(f(:,:,:,irho))
        if (lreference_state) &
          f(l1:l2,:,:,irho) = f(l1:l2,:,:,irho)-spread(spread(reference_state(:,iref_rho),2,my),3,mz)
      endif
!
    endsubroutine shell_ss_layers
!***********************************************************************
    subroutine piecew_poly_cylind(f)
!
!  Initialise ss in a cylindrical ring using 2 superposed polytropic layers.
!
!  17-mar-07/dintrans: coded
!
      use EquationOfState, only: lnrho0, cs20, gamma, gamma_m1, cs2top, &
                                 cs2bot, get_cp1, eoscalc, ilnrho_TT
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: lnrho, TT, ss
      real :: beta0, beta1, TT_bcz, TT_ext, TT_int
      real :: cp1, lnrho_bcz
      real, pointer :: gravx
!
      if (headtt) print*,'r_bcz in piecew_poly_cylind=', r_bcz
!
!  beta is the (negative) temperature gradient
!  beta = (g/cp) 1./[(1-1/gamma)*(m+1)]
!
      call get_shared_variable('gravx', gravx, caller='piecew_poly_cylind')
      call get_cp1(cp1)
      beta0=cp1*gravx/(mpoly0+1.)*gamma/gamma_m1
      beta1=cp1*gravx/(mpoly1+1.)*gamma/gamma_m1
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
        call putlnrho(f(:,m,n,ilnrho),lnrho)
        call eoscalc(ilnrho_TT,lnrho,TT,ss=ss)
        f(l1:l2,m,n,iss)=ss
      enddo
      enddo
!
    endsubroutine piecew_poly_cylind 
!***********************************************************************
    subroutine single_polytrope(f)
!
!  Note: both entropy and density are initialized there (compared to layer_ss)
!
!  06-sep-07/dintrans: coded a single polytrope of index mpoly0
!
      use EquationOfState, only: eoscalc, ilnrho_TT, get_cp1, &
                                 gamma_m1, lnrho0
      use Gravity, only: gravz
      use SharedVariables, only: get_shared_variable
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

        call putlnrho(f(:,m,n,ilnrho),lnrho)
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
    subroutine read_hcond(nloc,ngrid,ipc,hcondtop,hcondbot)
!
!  Read radial profiles of hcond and glhc from an ascii-file.
!
!  11-jun-09/pjk: coded
!  19-mar-14/pjk: added reading z-dependent profiles
!  14-jan-20/MR: simplified
!
      integer,            intent(in) :: nloc,ngrid,ipc
      real,               intent(out):: hcondbot, hcondtop
!
      real, dimension(ngrid) :: tmp1,tmp2
      logical :: exist
      integer :: stat,offset

      if (.not.allocated(hcond_prof)) allocate(hcond_prof(nloc),dlnhcond_prof(nloc))
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
          call fatal_error('read_hcond','no input file "hcond_glhc.*"')
        endif
      endif
!
!  Read profiles.
!
      do n=1,ngrid
        read(31,*,iostat=stat) tmp1(n),tmp2(n)
        if (stat<0) exit
        if (ip<5) print*,'hcond, glhc: ',tmp1(n),tmp2(n)
      enddo
!
      close(31)
!
!  Assuming no ghost zones in hcond_glhc.dat.
!
      offset=nloc*ipc
      hcond_prof=tmp1(offset+1:offset+nloc)
      dlnhcond_prof=tmp2(offset+1:offset+nloc)/hcond_prof
!
      hcondbot=tmp1(1)
      hcondtop=tmp1(ngrid)

      FbotKbot=Fbot/hcondbot
      FtopKtop=Ftop/hcondtop
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
      if (.not.lreference_state .and. entropy_floor /= impossible) &
        where(f(:,:,:,iss) < entropy_floor) f(:,:,:,iss) = entropy_floor
!
    endsubroutine impose_energy_floor
!***********************************************************************
    subroutine dynamical_thermal_diffusion(uc)
!
!  Dynamically set thermal diffusion coefficient given fixed mesh Reynolds number.
!
!  02-aug-11/ccyang: coded
!
!  Input Argument
!      uc
!          Characteristic velocity of the system.
!
      real, intent(in) :: uc
!
      if (chi_hyper3 /= 0.0) chi_hyper3 = pi5_1 * uc * dxmax**5 / re_mesh
      if (chi_hyper3_mesh /= 0.0) chi_hyper3_mesh = pi5_1 * uc / re_mesh / sqrt(real(dimensionality))
!
    endsubroutine dynamical_thermal_diffusion
!***********************************************************************
    subroutine get_lnQ(lnTT,lnQ,delta_lnTT)
!
!  17-may-2015/piyali.chatterjee : copied from special/solar_corona.f90
!  input: lnTT in SI units
!  output: lnP  [p]=W/s * m^3
!
      real, dimension (nx), intent(in) :: lnTT
      real, dimension (nx), intent(out) :: lnQ, delta_lnTT
!
      real, parameter, dimension (37) :: intlnT = (/ &
          8.74982, 8.86495, 8.98008, 9.09521, 9.21034, 9.44060, 9.67086, &
          9.90112, 10.1314, 10.2465, 10.3616, 10.5919, 10.8221, 11.0524, &
          11.2827, 11.5129, 11.7432, 11.9734, 12.2037, 12.4340, 12.6642, &
          12.8945, 13.1247, 13.3550, 13.5853, 13.8155, 14.0458, 14.2760, &
          14.5063, 14.6214, 14.7365, 14.8517, 14.9668, 15.1971, 15.4273, &
          15.6576,  69.0776 /)
      real, parameter, dimension (37) :: intlnQ = (/ &
          -93.9455, -91.1824, -88.5728, -86.1167, -83.8141, -81.6650, &
          -80.5905, -80.0532, -80.1837, -80.2067, -80.1837, -79.9765, &
          -79.6694, -79.2857, -79.0938, -79.1322, -79.4776, -79.4776, &
          -79.3471, -79.2934, -79.5159, -79.6618, -79.4776, -79.3778, &
          -79.4008, -79.5159, -79.7462, -80.1990, -80.9052, -81.3196, &
          -81.9874, -82.2023, -82.5093, -82.5477, -82.4172, -82.2637, &
          -0.66650 /)
      real, parameter, dimension(9) :: pars = (/ &
          2.12040e+00, 3.88284e-01, 2.02889e+00, 3.35665e-01, 6.34343e-01, &
          1.94052e-01, 2.54536e+00, 7.28306e-01, -2.40088e+01 /)
!
      integer ::  px, z_ref
      real :: pos, frac
!
        do px = 1, nx
          pos = interpol_tabulated (lnTT(px), intlnT)
          z_ref = floor (pos)
          if (z_ref < 1) then
            lnQ(px) = -max_real
            delta_lnTT(px) = intlnT(2) - intlnT(1)
            cycle
          endif
          if (z_ref > 36) z_ref = 36
          frac = pos - z_ref
          lnQ(px) = intlnQ(z_ref) * (1.0-frac) + intlnQ(z_ref+1) * frac
          delta_lnTT(px) = intlnT(z_ref+1) - intlnT(z_ref)
        enddo
!
    endsubroutine get_lnQ
!***********************************************************************
    function interpol_tabulated (needle, haystack)
!
! Find the interpolated position of a given value in a tabulated values array.
! Bisection search algorithm with preset range guessing by previous value.
! Returns the interpolated position of the needle in the haystack.
! If needle is not inside the haystack, an extrapolated position is returned.
!
! 09-feb-2011/Bourdin.KIS: coded
! 17-may-2015/piyali.chatterjee : copied from special/solar_corona.f90

      real :: interpol_tabulated
      real, intent(in) :: needle
      real, dimension (:), intent(in) :: haystack
!
      integer, save :: lower=1, upper=1
      integer :: mid, num, inc
!
      num = size (haystack, 1)
      if (num < 2) call fatal_error ('interpol_tabulated', "Too few tabulated values!", .true.)
      if (lower >= num) lower = num - 1
      if ((upper <= lower) .or. (upper > num)) upper = num
!
      if (haystack(lower) > haystack(upper)) then
!
!  Descending array:
!
        ! Search for lower limit, starting from last known position
        inc = 2
        do while ((lower > 1) .and. (needle > haystack(lower)))
          upper = lower
          lower = lower - inc
          if (lower < 1) lower = 1
          inc = inc * 2
        enddo
!
        ! Search for upper limit, starting from last known position!
        inc = 2
        do while ((upper < num) .and. (needle < haystack(upper)))
          lower = upper
          upper = upper + inc
          if (upper > num) upper = num
          inc = inc * 2
        enddo
!
        if (needle < haystack(upper)) then
          ! Extrapolate needle value below range
          lower = num - 1
        elseif (needle > haystack(lower)) then
          ! Extrapolate needle value above range
          lower = 1
        else
          ! Interpolate needle value
          do while (lower+1 < upper)
            mid = lower + (upper - lower) / 2
            if (needle >= haystack(mid)) then
              upper = mid
            else
              lower = mid
            endif
          enddo
        endif
        upper = lower + 1
        interpol_tabulated = lower + (haystack(lower) - needle) / (haystack(lower) - haystack(upper))
!
      elseif (haystack(lower) < haystack(upper)) then
!
!  Ascending array:
!
        ! Search for lower limit, starting from last known position
        inc = 2
        do while ((lower > 1) .and. (needle < haystack(lower)))
          upper = lower
          lower = lower - inc
          if (lower < 1) lower = 1
          inc = inc * 2
        enddo
!
        ! Search for upper limit, starting from last known position
        inc = 2
        do while ((upper < num) .and. (needle > haystack(upper)))
          lower = upper
          upper = upper + inc
          if (upper > num) upper = num
          inc = inc * 2
        enddo
!
        if (needle > haystack(upper)) then
          ! Extrapolate needle value above range
          lower = num - 1
        elseif (needle < haystack(lower)) then
          ! Extrapolate needle value below range
          lower = 1
        else
          ! Interpolate needle value
          do while (lower+1 < upper)
            mid = lower + (upper - lower) / 2
            if (needle < haystack(mid)) then
              upper = mid
            else
              lower = mid
            endif
          enddo
        endif
        upper = lower + 1
        interpol_tabulated = lower + (needle - haystack(lower)) / (haystack(upper) - haystack(lower))
      else
        interpol_tabulated = -1.0
        call fatal_error ('interpol_tabulated', "Tabulated values are invalid!", .true.)
      endif
!
    endfunction interpol_tabulated
!***********************************************************************
    subroutine split_update_energy(f)
!
!  Dummy subroutine
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine
!***********************************************************************
   subroutine energy_after_timestep(f,df,dtsub)
!
     real, dimension(mx,my,mz,mfarray) :: f
     real, dimension(mx,my,mz,mvar) :: df
     real :: dtsub
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dtsub)
!
   endsubroutine energy_after_timestep
!***********************************************************************
    subroutine expand_shands_energy
!
!  Presently dummy, for possible use
!
    endsubroutine expand_shands_energy
!***********************************************************************
    subroutine pushdiags2c(p_diag)

    integer, parameter :: n_diags=0
    integer(KIND=ikind8), dimension(:) :: p_diag

    call keep_compiler_quiet(p_diag)

    endsubroutine pushdiags2c
!***********************************************************************
    subroutine pushpars2c(p_par)

    integer, parameter :: n_pars=1
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr_c(chi,p_par(1))

    endsubroutine pushpars2c
!***********************************************************************
endmodule Energy
