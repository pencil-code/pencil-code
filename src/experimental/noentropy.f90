! $Id: entropy.f90 11900 2009-10-13 23:02:31Z boris.dintrans $
! 
!  This module takes care of entropy (initial condition
!  and time advance)
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lentropy = .false.
! CPARAM logical, parameter :: ltemperature = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED ugss; Ma2; fpres(3); uglnTT
!
!***************************************************************
module Entropy
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include '../entropy.h'
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
  !parameters for Sedov type initial condition
  real :: center1_x=0., center1_y=0., center1_z=0.
  real :: center2_x=0., center2_y=0., center2_z=0.
  real :: kx_ss=1.,ky_ss=1.,kz_ss=1.
  real :: thermal_background=0., thermal_peak=0., thermal_scaling=1.
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
  real :: tau_cool=0.0, TTref_cool=0.0
  real :: cs0hs=0.0,H0hs=0.0,rho0hs=0.0
  real :: chit_aniso=0.0,xbot=0.0
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
  logical, pointer :: lpressuregradient_gas ! Shared with Hydro.
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
  character (len=5) :: iinit_str
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
  integer :: idiag_thermalpressure=0    ! DIAG_DOC:
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
  integer :: idiag_ssmxy=0      ! DIAG_DOC: $\left< s \right>_{z}$
  integer :: idiag_ssmxz=0      ! DIAG_DOC: $\left< s \right>_{y}$


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
      use SharedVariables
      use Sub
!
      integer :: ierr
!
!
!  identify version number
!
      if (lroot) call svn_id( &
          "$Id: entropy_anelastic.f90 11900 2009-10-13 23:02:31Z dhruba.mitra $")
!
!  Get the shared variable lpressuregradient_gas from Hydro module.
!
      call get_shared_variable('lpressuregradient_gas',lpressuregradient_gas,ierr)     
      if (ierr/=0) call fatal_error('register_entropy','there was a problem getting lpressuregradient_gas')
!
    endsubroutine register_entropy
!***********************************************************************
    subroutine initialize_entropy(f,lstarting)
!
!  Called by run.f90 after reading parameters, but before the time loop.
!
!  21-jul-02/wolf: coded
!
      use EquationOfState, only: cs0,  &
                                 beta_glnrho_global, beta_glnrho_scaled, &
                                 select_eos_variable,gamma,gamma_m1
      use FArrayManager
      use Gravity, only: gravz,g0
      use Mpicomm, only: stop_it
      use SharedVariables, only: put_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      real, dimension (nx,3) :: glhc
      real, dimension (nx) :: hcond
      real :: beta1, cp1, beta0, TT_crit
      integer :: i, ierr, q
      logical :: lnothing,lcompute_grav
      type (pencil_case) :: p
!
! Check any module dependencies
!
      if (.not. leos) then
        call fatal_error('initialize_entropy', &
            'EOS=noeos but entropy requires an EQUATION OF STATE for the fluid')
      endif
!
! Tell the equation of state that we're here and what f variable we use
!
        if (pretend_lnTT) then
          call select_eos_variable('lnTT',iss)
        else
          if (gamma_m1==0.) then
            call select_eos_variable('lnTT',-1)
          else
            call select_eos_variable('ss',iss)
          endif
        endif
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
      call put_shared_variable('lviscosity_heat',lviscosity_heat,ierr)
      if (ierr/=0) call stop_it("initialize_entropy: "//&
           "there was a problem when putting lviscosity_heat")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
      endsubroutine initialize_entropy
!***********************************************************************
    subroutine read_entropy_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_entropy_init_pars
!***********************************************************************
    subroutine read_entropy_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_entropy_run_pars
!***********************************************************************
    subroutine write_entropy_init_pars(unit)
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_entropy_init_pars
!***********************************************************************
    subroutine write_entropy_run_pars(unit)
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_entropy_run_pars
!***********************************************************************
    subroutine init_ss(f)
!
!  initialise entropy; called from start.f90
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_ss
!***********************************************************************
    subroutine pencil_criteria_entropy()
!
!  All pencils that the Entropy module depends on are specified here.
!
!  20-11-04/anders: coded
!
      use EquationOfState, only: beta_glnrho_scaled
!
      if (ldt) lpenc_requested(i_cs2)=.true.
      if (lpressuregradient_gas) lpenc_requested(i_fpres)=.true.
      if (ladvection_entropy) lpenc_requested(i_ugss)=.true.
!
      if (maxval(abs(beta_glnrho_scaled))/=0.0) lpenc_requested(i_cs2)=.true.
!
      if (idiag_ugradpm/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_uglnrho)=.true.
      endif
!
      if (idiag_thermalpressure/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_cs2)=.true.
        lpenc_diagnos(i_rcyl_mn)=.true.
      endif
!
      if (idiag_ethm/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_ee)=.true.
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
    endsubroutine pencil_interdep_entropy
!***********************************************************************
    subroutine calc_pencils_entropy(f,p)
!
! Do nothing 
! DM+PC

      use EquationOfState, only: gamma,gamma_m1,cs20,lnrho0,profz_eos
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!    
    endsubroutine calc_pencils_entropy
!**********************************************************************
    subroutine dss_dt(f,df,p)
! Do nothing
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      integer :: j,ju
      integer :: ix
!
      intent(in) :: f,p,df
!      intent(out) :: df
!

      call keep_compiler_quiet(f)
!

    endsubroutine dss_dt
!***********************************************************************
    subroutine calc_pencils_entropy_after_mn(f,p)
!
! Copied cal_pencils_entropy from noentropy to entropy_anelastic
! DM+PC
!
!  Calculate Entropy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-nov-04/anders: coded
!
      use EquationOfState, only: gamma,gamma_m1,cs20,lnrho0,profz_eos
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
! fpres (=pressure gradient force)
      if (lpencil(i_fpres)) then
        do j=1,3
          if (llocal_iso) then
            p%fpres(:,j)=-p%cs2*(p%glnrho(:,j)+p%glnTT(:,j))
          elseif (ldensity_anelastic) then
            p%fpres(:,j)=-p%gpp(:,j)/exp(p%lnrho)
          endif
          if (profz_eos(n)/=1.0) p%fpres(:,j)=profz_eos(n)*p%fpres(:,j)
        enddo
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    
    endsubroutine calc_pencils_entropy_after_mn
!**********************************************************************
    subroutine dss_dt_after_mn(f,df,p)
! Added dss_dt from noentropy to dss_dt_after_mn
! DM+PC
!
!  Calculate pressure gradient term for isothermal/polytropic equation
!  of state in the anelastic case.
!
      use EquationOfState, only: beta_glnrho_global, beta_glnrho_scaled
      use Diagnostics
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      integer :: j,ju
!
      intent(in) :: f,p
      intent(out) :: df
!
!  ``cs2/dx^2'' for timestep
!
      if (leos.and.ldensity) then ! no sound waves without equation of state
        if (lfirst.and.ldt) advec_cs2=p%cs2*dxyz_2
        if (headtt.or.ldebug) print*,'dss_dt: max(advec_cs2) =',maxval(advec_cs2)
      endif
!
!  Add isothermal/polytropic pressure term in momentum equation
!
      if (lhydro .and. lpressuregradient_gas) then
        do j=1,3
          ju=j+iuu-1
          df(l1:l2,m,n,ju)=df(l1:l2,m,n,ju)+p%fpres(:,j)
        enddo
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
     endif
!
!  Calculate entropy related diagnostics
!
      if (ldiagnos) then
        if (idiag_dtc/=0) &
            call max_mn_name(sqrt(advec_cs2)/cdt,idiag_dtc,l_dt=.true.)
        if (idiag_ugradpm/=0) &
            call sum_mn_name(p%rho*p%cs2*p%uglnrho,idiag_ugradpm)
        if (idiag_ethm/=0) call sum_mn_name(p%rho*p%ee,idiag_ethm)
      endif
!

      call keep_compiler_quiet(f)
!
    endsubroutine dss_dt_after_mn
!***********************************************************************
    subroutine rprint_entropy(lreset,lwrite)
!
!  reads and registers print parameters relevant to entropy
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Diagnostics
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
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_dtc=',idiag_dtc
        write(3,*) 'i_dtchi=',idiag_dtchi
        write(3,*) 'i_ethtot=',idiag_ethtot
        write(3,*) 'i_ethdivum=',idiag_ethdivum
        write(3,*) 'i_ethm=',idiag_ethm
        write(3,*) 'i_ssm=',idiag_ssm
        write(3,*) 'i_ss2m=',idiag_ss2m
        write(3,*) 'i_eem=',idiag_eem
        write(3,*) 'i_ppm=',idiag_ppm
        write(3,*) 'i_pdivum=',idiag_pdivum
        write(3,*) 'i_heatm=',idiag_heatm
        write(3,*) 'i_csm=',idiag_csm
        write(3,*) 'i_fconvm=',idiag_fconvm
        write(3,*) 'i_ugradpm=',idiag_ugradpm
        write(3,*) 'i_fradbot=',idiag_fradbot
        write(3,*) 'i_fradtop=',idiag_fradtop
        write(3,*) 'i_TTtop=',idiag_TTtop
        write(3,*) 'i_ssmphi=',idiag_ssmphi
        write(3,*) 'i_cs2mphi=',idiag_cs2mphi
        write(3,*) 'i_fturbz=',idiag_fturbz
        write(3,*) 'i_fconvz=',idiag_fconvz
        write(3,*) 'i_dcoolz=',idiag_dcoolz
        write(3,*) 'i_fradz=',idiag_fradz
        write(3,*) 'i_ssmz=',idiag_ssmz
        write(3,*) 'i_TTmz=',idiag_TTmz
        write(3,*) 'i_uxTTmz=',idiag_uxTTmz
        write(3,*) 'i_uyTTmz=',idiag_uyTTmz
        write(3,*) 'i_uzTTmz=',idiag_uzTTmz
        write(3,*) 'i_ssmr=',idiag_ssmr
        write(3,*) 'i_TTmr=',idiag_TTmr
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
        write(3,*) 'i_TTmxy=',idiag_TTmxy
        write(3,*) 'i_TTmxz=',idiag_TTmxz
        write(3,*) 'i_ssmxy=',idiag_ssmxy
        write(3,*) 'i_ssmxz=',idiag_ssmxz
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
      use EquationOfState, only: eoscalc, ilnrho_ss
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
        case ('ss')
          slices%yz =f(ix_loc,m1:m2,n1:n2,iss)
          slices%xz =f(l1:l2,iy_loc,n1:n2,iss)
          slices%xy =f(l1:l2,m1:m2,iz_loc,iss)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,iss)
          if (lwrite_slice_xy3) slices%xy3=f(l1:l2,m1:m2,iz3_loc,iss)
          if (lwrite_slice_xy4) slices%xy4=f(l1:l2,m1:m2,iz4_loc,iss)
          slices%ready=.true.
!
!
      endselect
!   
    endsubroutine get_slices_entropy
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
         if (lfirst .and. ip == 13) &
              call output_pencil(trim(directory)//'/newton.dat',newton*exp(p%lnrho+p%lnTT),1)
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
endmodule Entropy
