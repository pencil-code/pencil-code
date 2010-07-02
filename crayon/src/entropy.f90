! $Id: entropy.f90 13610 2010-04-09 14:34:59Z sven.bingert $
!
!  This module takes care of evolving the entropy.
!
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
! PENCILS PROVIDED ugss; Ma2; fpres(3); uglnTT
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
  use Viscosity
!
  implicit none
!
  include 'entropy.h'
!
  real :: radius_ss=0.1, ampl_ss=0.0, widthss=2*epsi, epsilon_ss=0.0
  real :: luminosity=0.0, wheat=0.1, cool=0.0, zcool=0.0, rcool=0.0, wcool=0.1
  real :: TT_int, TT_ext, cs2_int, cs2_ext
  real :: cool_int=0.0, cool_ext=0.0, ampl_TT=0.0
  real, target :: chi=0.0
  real :: chi_t=0.0, chi_shock=0.0, chi_th=0.0
  real :: Kgperp=0.0, Kgpara=0.0, tdown=0.0, allp=2.0
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
  real, target :: hcond0=impossible, hcond1=impossible
  real, target :: Fbot=impossible, FbotKbot=impossible
  real, target :: Ftop=impossible, FtopKtop=impossible
  real :: Kbot=impossible, Ktop=impossible
  real :: hcond2=impossible
  real :: chit_prof1=1.0, chit_prof2=1.0
  real :: tau_cor=0.0, TT_cor=0.0, z_cor=0.0
  real :: tauheat_buffer=0.0, TTheat_buffer=0.0
  real :: zheat_buffer=0.0, dheat_buffer1=0.0
  real :: heat_uniform=0.0, cool_RTV=0.0
  real :: deltaT_poleq=0.0, beta_hand=1.0, r_bcz=0.0
  real :: tau_cool=0.0, tau_diff=0.0, TTref_cool=0.0, tau_cool2=0.0
  real :: cs0hs=0.0, H0hs=0.0, rho0hs=0.0
  real :: xbot=0.0, xtop=0.0
  integer, parameter :: nheatc_max=4
  integer :: iglobal_hcond=0
  integer :: iglobal_glhc=0
  integer :: ippaux=0
  logical :: lturbulent_heat=.false.
  logical :: lheatc_Kprof=.false., lheatc_Kconst=.false.
  logical, target :: lheatc_chiconst=.false.
  logical :: lheatc_tensordiffusion=.false.
  logical :: lheatc_chitherm=.false.
  logical :: lheatc_shock=.false.
  logical :: lcooling_general=.false.
  logical, target :: lmultilayer=.true.
  logical :: ladvection_entropy=.true.
  logical, pointer :: lpressuregradient_gas
  logical :: lviscosity_heat=.true.
  logical :: lfreeze_sint=.false.,lfreeze_sext=.false.
  logical :: lhcond_global=.false.
  logical :: lfpres_from_pressure=.false.
  character (len=labellen), dimension(ninit) :: initss='nothing'
  character (len=labellen) :: pertss='zero'
  character (len=labellen) :: cooltype='Temp',cooling_profile='gaussian'
  character (len=labellen), dimension(nheatc_max) :: iheatcond='nothing'
  character (len=5) :: iinit_str
!
!  Input parameters.
!
  namelist /entropy_init_pars/ &
      initss, pertss, grads0, radius_ss, ampl_ss, widthss, epsilon_ss, &
      mixinglength_flux, chi_t, chi_th, pp_const, ss_left, ss_right, ss_const, mpoly0, &
      mpoly1, mpoly2, isothtop, khor_ss, thermal_background, thermal_peak, &
      thermal_scaling, cs2cool, center1_x, center1_y, center1_z, center2_x, &
      center2_y, center2_z, T0, ampl_TT, kx_ss, ky_ss, kz_ss, &
      beta_glnrho_global, ladvection_entropy, lviscosity_heat, r_bcz, &
      luminosity, wheat, hcond0, tau_cool, TTref_cool, lhcond_global, &
      cool_fac, cs0hs, H0hs, rho0hs, tau_cool2
!
!  Run parameters.
!
  namelist /entropy_run_pars/ &
      hcond0, hcond1, hcond2, widthss, mpoly0, mpoly1, mpoly2, &
      luminosity, wheat, cooling_profile, cooltype, cool, cs2cool, rcool, &
      wcool, Fbot, lcooling_general, chi_t, chi_th, chit_prof1, chit_prof2, chi_shock, &
      chi, iheatcond, Kgperp, Kgpara, cool_RTV, tau_ss_exterior, lmultilayer, &
      Kbot, tau_cor, TT_cor, z_cor, tauheat_buffer, TTheat_buffer, &
      zheat_buffer, dheat_buffer1, heat_uniform, cool_int, cool_ext, &
      lturbulent_heat, deltaT_poleq, tdown, allp, &
      beta_glnrho_global, ladvection_entropy, lviscosity_heat, r_bcz, &
      lfreeze_sint, lfreeze_sext, lhcond_global, tau_cool, TTref_cool, &
      mixinglength_flux, chiB, Ftop, xbot, xtop, tau_cool2, &
      tau_diff, lfpres_from_pressure
!
!  Diagnostic variables (need to be consistent with reset list below).
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
  integer :: idiag_yHm=0        ! DIAG_DOC:
  integer :: idiag_yHmax=0      ! DIAG_DOC:
  integer :: idiag_TTm=0        ! DIAG_DOC:
  integer :: idiag_TTmax=0      ! DIAG_DOC:
  integer :: idiag_TTmin=0      ! DIAG_DOC:
  integer :: idiag_fconvm=0     ! DIAG_DOC:
  integer :: idiag_fconvz=0     ! DIAG_DOC:
  integer :: idiag_dcoolz=0     ! DIAG_DOC:
  integer :: idiag_fradz=0      ! DIAG_DOC:
  integer :: idiag_fradz_Kprof=0 ! DIAG_DOC:
  integer :: idiag_fradxy_Kprof=0 ! DIAG_DOC:
  integer :: idiag_fturbz=0     ! DIAG_DOC:
  integer :: idiag_fturbxy=0    ! DIAG_DOC:
  integer :: idiag_ssmx=0       ! DIAG_DOC:
  integer :: idiag_ssmy=0       ! DIAG_DOC:
  integer :: idiag_ssmz=0       ! DIAG_DOC:
  integer :: idiag_ppmx=0       ! DIAG_DOC:
  integer :: idiag_ppmy=0       ! DIAG_DOC:
  integer :: idiag_ppmz=0       ! DIAG_DOC:
  integer :: idiag_TTp=0        ! DIAG_DOC:
  integer :: idiag_TTmx=0       ! DIAG_DOC:
  integer :: idiag_TTmy=0       ! DIAG_DOC:
  integer :: idiag_TTmz=0       ! DIAG_DOC:
  integer :: idiag_TTmxy=0      ! DIAG_DOC:
  integer :: idiag_TTmxz=0      ! DIAG_DOC:
  integer :: idiag_uxTTmz=0     ! DIAG_DOC:
  integer :: idiag_uyTTmz=0     ! DIAG_DOC:
  integer :: idiag_uzTTmz=0     ! DIAG_DOC:
  integer :: idiag_uxTTmxy=0    ! DIAG_DOC:
  integer :: idiag_ssmxy=0      ! DIAG_DOC: $\left< s \right>_{z}$
  integer :: idiag_ssmxz=0      ! DIAG_DOC: $\left< s \right>_{y}$
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
      use FArrayManager
      use SharedVariables
      use Sub
!
      integer :: ierr
!
      call farray_register_pde('ss',iss)
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id: entropy.f90 13610 2010-04-09 14:34:59Z sven.bingert $")
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
      use EquationOfState, only: cs0, get_soundspeed, get_cp1, &
                                 beta_glnrho_global, beta_glnrho_scaled, &
                                 mpoly, mpoly0, mpoly1, mpoly2, &
                                 select_eos_variable,gamma,gamma_m1
      use FArrayManager
      use Initcond
      use Gravity, only: gravz, g0
      use Mpicomm, only: stop_it
      use SharedVariables, only: put_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      real, dimension (nx) :: hcond
      integer :: i, ierr
      logical :: lnothing
!
!  Check any module dependencies.
!
      if (.not. leos) then
        call fatal_error('initialize_entropy', &
            'EOS=noeos but entropy requires an EQUATION OF STATE for the fluid')
      endif
!
!  Tell the equation of state that we're here and what f variable we use.
!
      call select_eos_variable('ss',iss)
!
!  Radiative diffusion: initialize flux etc.
!
      hcond=0.0
!
!  Kbot and hcond0 are used interchangibly, so if one is
!  =impossible, set it to the other's value
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
      lheatc_chiconst=.false.
      lheatc_shock=.false.
!
      lnothing=.false.
!
!  select which radiative heating we are using
!
      if (lroot) print*,'initialize_entropy: nheatc_max,iheatcond=',nheatc_max,iheatcond(1:nheatc_max)
      do i=1,nheatc_max
        select case (iheatcond(i))
        case ('K-const')
          lheatc_Kconst=.true.
          if (lroot) print*, 'heat conduction: K=cte'
        case ('chi-const')
          lheatc_chiconst=.true.
          if (lroot) print*, 'heat conduction: constant chi'
        case ('shock')
          lheatc_shock=.true.
          if (lroot) print*, 'heat conduction: shock'
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
      if (lheatc_chiconst .and. (chi==0.0 .and. chi_t==0.0)) then
        call warning('initialize_entropy','chi and chi_t are zero!')
      endif
      if (all(iheatcond=='nothing') .and. hcond0/=0.0) then
        call warning('initialize_entropy', &
            'No heat conduction, but hcond0 /= 0')
      endif
      if (lheatc_Kconst .and. Kbot==0.0) then
        call warning('initialize_entropy','Kbot is zero!')
      endif
      if (lheatc_shock .and. chi_shock==0.0) then
        call warning('initialize_entropy','chi_shock is zero!')
      endif
!
! Shared variables
!
      call put_shared_variable('hcond0',hcond0,ierr)
      if (ierr/=0) call stop_it("initialize_entropy: "//&
           "there was a problem when putting hcond0")
      call put_shared_variable('hcond1',hcond1,ierr)
      if (ierr/=0) call stop_it("initialize_entropy: "//&
           "there was a problem when putting hcond1")
      call put_shared_variable('Fbot',Fbot,ierr)
      if (ierr/=0) call stop_it("initialize_entropy: "//&
           "there was a problem when putting Fbot")
      call put_shared_variable('Ftop',Ftop,ierr)
      if (ierr/=0) call stop_it("initialize_entropy: "//&
           "there was a problem when putting Ftop")
      call put_shared_variable('FbotKbot',FbotKbot,ierr)
      if (ierr/=0) call stop_it("initialize_entropy: "//&
           "there was a problem when putting FbotKbot")
      call put_shared_variable('FtopKtop',FtopKtop,ierr)
      if (ierr/=0) call stop_it("initialize_entropy: "//&
           "there was a problem when putting FtopKtop")
      call put_shared_variable('chi',chi,ierr)
      if (ierr/=0) call stop_it("initialize_entropy: "//&
           "there was a problem when putting chi")
      call put_shared_variable('chi_t',chi_t,ierr)
      if (ierr/=0) call stop_it("initialize_entropy: "//&
           "there was a problem when putting chi_t")
      call put_shared_variable('chi_th',chi_th,ierr)
      if (ierr/=0) call stop_it("initialize_entropy: "//&
           "there was a problem when putting chi_th")
      call put_shared_variable('lmultilayer',lmultilayer,ierr)
      if (ierr/=0) call stop_it("initialize_entropy: "//&
           "there was a problem when putting lmultilayer")
      call put_shared_variable('lheatc_chiconst',lheatc_chiconst,ierr)
      if (ierr/=0) call stop_it("initialize_entropy: "//&
           "there was a problem when putting lcalc_heatcond_constchi")
      call put_shared_variable('lheatc_chitherm',lheatc_chitherm,ierr)
      if (ierr/=0) call stop_it("initialize_entropy: "//&
           "there was a problem when putting lcalc_heatcond_chitherm")
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
!  initialise entropy; called from start.f90
!  07-nov-2001/wolf: coded
!  24-nov-2002/tony: renamed for consistancy (i.e. init_[variable name])
!
      use Sub
      use Gravity
      use Initcond
      use General, only: chn
      use InitialCondition, only: initial_condition_ss
      use EquationOfState,  only: isothtop, &
                                mpoly0, mpoly1, mpoly2, cs2cool, cs0, &
                                rho0, lnrho0, &
                                eoscalc, ilnrho_pp, &
                                eosperturb
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      integer :: j
      logical :: lnothing=.true.
!
      intent(inout) :: f
!
      do j=1,ninit
!
        if (initss(j)=='nothing') cycle
!
        lnothing=.false.
        call chn(j,iinit_str)
!
!  select different initial conditions
!
        select case (initss(j))

          case ('zero', '0'); f(:,:,:,iss) = 0.
          case ('const_ss'); f(:,:,:,iss)=f(:,:,:,iss)+ss_const
          case ('gaussian-noise'); call gaunoise(ampl_ss,f,iss,iss)
          case ('blob')
            call blob(ampl_ss,f,iss,radius_ss,center1_x,center1_y,center1_z)
          case ('wave')
            do n=n1,n2; do m=m1,m2
              f(l1:l2,m,n,iss)=f(l1:l2,m,n,iss)+ss_const+ampl_ss*sin(kx_ss*x(l1:l2)+pi)
            enddo; enddo
          case('xjump'); call jump(f,iss,ss_left,ss_right,widthss,'x')
          case('yjump'); call jump(f,iss,ss_left,ss_right,widthss,'y')
          case('zjump'); call jump(f,iss,ss_left,ss_right,widthss,'z')
          case('sinxsinz'); call sinxsinz(ampl_ss,f,iss,kx_ss,ky_ss,kz_ss)
          case('hor-fluxtube')
            call htube(ampl_ss,f,iss,iss,radius_ss,epsilon_ss,center1_x,center1_z)
          case ('hor-tube')
            call htube2(ampl_ss,f,iss,iss,radius_ss,epsilon_ss)
          case ('sedov')
            if (lroot) print*,'init_ss: sedov - thermal background with gaussian energy burst'
          call blob(thermal_peak,f,iss,radius_ss,center1_x,center1_y,center1_z)
          case ('sedov-dual')
            if (lroot) print*,'init_ss: sedov - thermal background with gaussian energy burst'
          call blob(thermal_peak,f,iss,radius_ss,center1_x,center1_y,center1_z)
          call blob(thermal_peak,f,iss,radius_ss,center2_x,center2_y,center2_z)
          case ('blob_hs')
            print*,'init_ss: put blob in hydrostatic equilibrium: radius_ss,ampl_ss=',radius_ss,ampl_ss
            call blob(ampl_ss,f,iss,radius_ss,center1_x,center1_y,center1_z)
            call blob(-ampl_ss,f,ilnrho,radius_ss,center1_x,center1_y,center1_z)
          case default
!
!  Catch unknown values
!
            write(unit=errormsg,fmt=*) 'No such value for initss(' &
                             //trim(iinit_str)//'): ',trim(initss(j))
            call fatal_error('init_ss',errormsg)
        endselect
!
        if (lroot) print*,'init_ss: initss('//trim(iinit_str)//') = ', &
            trim(initss(j))
!
      enddo
!
      if (lnothing.and.lroot) print*,'init_ss: nothing'
!
!  Interface fow user's own initial condition
!
      if (linitial_condition) call initial_condition_ss(f)
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
      if (lheatc_Kconst .or. lheatc_chiconst .or. lheatc_Kprof .or. &
          tau_cor>0 .or. lheatc_chitherm) lpenc_requested(i_cp1)=.true.
      if (ldt) lpenc_requested(i_cs2)=.true.
      if (lpressuregradient_gas) lpenc_requested(i_fpres)=.true.
      if (ladvection_entropy) lpenc_requested(i_ugss)=.true.
      if (lviscosity.and.lviscosity_heat) then
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_visc_heat)=.true.
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
      if (lheatc_Kprof) then
        if (hcond0/=0) then
          lpenc_requested(i_rho1)=.true.
          lpenc_requested(i_glnTT)=.true.
          lpenc_requested(i_del2lnTT)=.true.
          lpenc_requested(i_cp1)=.true.
          if (lmultilayer) then
            if (lgravz) then
              lpenc_requested(i_z_mn)=.true.
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
        if (chi_th/=0) then
           lpenc_requested(i_del2ss)=.true.
           lpenc_requested(i_glnrho)=.true.
           lpenc_requested(i_gss)=.true.
        endif
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
      if (cooltype=='shell' .and. deltaT_poleq/=0.) then
        lpenc_requested(i_z_mn)=.true.
        lpenc_requested(i_rcyl_mn)=.true.
      endif
      if (tau_cool/=0.0) then
        lpenc_requested(i_cp)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_rho)=.true.
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
      if (idiag_pdivum/=0) then
        lpenc_diagnos(i_pp)=.true.
        lpenc_diagnos(i_divu)=.true.
      endif
      if (idiag_ssm/=0 .or. idiag_ss2m/=0 .or. idiag_ssmz/=0 .or. &
          idiag_ssmy/=0 .or. idiag_ssmx/=0) &
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
      if (idiag_fturbxy/=0) then
          lpenc_diagnos2d(i_rho)=.true.
          lpenc_diagnos2d(i_TT)=.true.
      endif
      if (idiag_TTm/=0 .or. idiag_TTmx/=0 .or. idiag_TTmy/=0 .or. &
          idiag_TTmz/=0 .or. idiag_TTmax/=0 .or. &
          idiag_TTmin/=0 .or. idiag_uxTTmz/=0 .or.idiag_uyTTmz/=0 .or. &
          idiag_uzTTmz/=0) &
          lpenc_diagnos(i_TT)=.true.
      if (idiag_TTmxy/=0 .or. idiag_TTmxz/=0 .or. idiag_uxTTmxy/=0) &
          lpenc_diagnos2d(i_TT)=.true.
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
          call u_dot_grad(iss,p%gss,p%uu,p%ugss)
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
!
    endsubroutine calc_pencils_entropy
!***********************************************************************
    subroutine dss_dt(df,p)
!
!  Calculate right hand side of entropy equation,
!  ds/dt = -u.grads + [H-C + div(K*gradT) + mu0*eta*J^2 + ...]
!
!  17-sep-01/axel: coded
!   9-jun-02/axel: pressure gradient added to du/dt already here
!   2-feb-03/axel: added possibility of ionization
!
      use Diagnostics
      use EquationOfState, only: beta_glnrho_global, beta_glnrho_scaled, gamma_inv, cs0
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: Hmax
      real :: uT
      integer :: j,ju
!
      intent(in)  :: p
      intent(out) :: df
!
      Hmax = 0.
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
      endif
!
!  Advection of entropy.
!
      if (ladvection_entropy) &
           df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) - p%ugss
!
!  Calculate viscous contribution to entropy
!
      if (lviscosity .and. lviscosity_heat) call calc_viscous_heat(df,p,Hmax)
!
!  Thermal conduction delegated to different subroutines.
!
      if (lheatc_Kconst)   call calc_heatcond_constK(df,p)
      if (lheatc_chiconst) call calc_heatcond_constchi(df,p)
      if (lheatc_shock)    call calc_heatcond_shock(df,p)
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
      endif
!
!  1-D averages.
!
      if (l1davgfirst) then
        if (idiag_fradz/=0) call xysum_mn_name_z(-hcond0*p%TT*p%glnTT(:,3),idiag_fradz)
        if (idiag_fconvz/=0) &
            call xysum_mn_name_z(p%rho*p%uu(:,3)*p%TT,idiag_fconvz)
        if (idiag_ssmx/=0)  call yzsum_mn_name_x(p%ss,idiag_ssmx)
        if (idiag_ssmy/=0)  call xzsum_mn_name_y(p%ss,idiag_ssmy)
        if (idiag_ssmz/=0)  call xysum_mn_name_z(p%ss,idiag_ssmz)
        if (idiag_ppmx/=0)  call yzsum_mn_name_x(p%pp,idiag_ppmx)
        if (idiag_ppmy/=0)  call xzsum_mn_name_y(p%pp,idiag_ppmy)
        if (idiag_ppmz/=0)  call xysum_mn_name_z(p%pp,idiag_ppmz)
        if (idiag_TTmx/=0)  call yzsum_mn_name_x(p%TT,idiag_TTmx)
        if (idiag_TTmy/=0)  call xzsum_mn_name_y(p%TT,idiag_TTmy)
        if (idiag_TTmz/=0)  call xysum_mn_name_z(p%TT,idiag_TTmz)
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
    endsubroutine dss_dt
!**********************************************************************
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
!
!  add heat conduction to entropy equation
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
      thdiff=chi_shock*(p%shock*(p%del2lnTT+g2)+gshockglnTT)
      if (chi_t/=0.) then
        call dot(p%glnrho+p%glnTT,p%gss,g2)
        thdiff=thdiff+chi_t*(p%del2ss+g2)
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
      use Diagnostics, only: max_mn_name
      use Sub, only: dot
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

      !  diffusion of the form:
      !  rho*T*Ds/Dt = ... + nab.(K*gradT)
      !        Ds/Dt = ... + K/rho*[del2lnTT+(glnTT)^2]
      !
      ! NB: chix = K/(cp rho) is needed for diffus_chi calculation
      !
      !  put empirical heat transport suppression by the B-field
      !
      if (chiB==0.) then
        chix = p%rho1*hcond*p%cp1
      else
        chix = p%rho1*hcond*p%cp1/(1.+chiB*p%b2)
      endif
      call dot(p%glnTT,p%glnTT,g2)
      !
      thdiff = p%rho1*hcond * (p%del2lnTT + g2)
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
    subroutine calc_heat_cool(df,p,Hmax)
!
!  Add combined heating and cooling to entropy equation.
!
!  02-jul-02/wolf: coded
!
      use Diagnostics, only: sum_mn_name, xysum_mn_name_z
      use Gravity, only: z2
      use IO, only: output_pencil
      use Sub, only: step, cubic_step, write_zprof
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: Hmax
!
      real, dimension (nx) :: heat,prof
      real :: zbot,ztop,profile_buffer
!
      intent(in) :: p
      intent(out) :: df
!
!  Vertical gravity determines some heat/cool models.
!
      if (headtt) print*, 'calc_heat_cool: lgravz= lgravx=', lgravz, lgravx
!
!  Initialize heating/cooling term.
!
      heat=0.
!
!  General spatially distributed cooling profiles (independent of gravity)
!
      if (lcooling_general) then
        select case (cooling_profile)
        case ('gaussian-z')
          prof=spread(exp(-0.5*((zcool-z(n))/wcool)**2), 1, l2-l1+1)
        endselect
        !heat=heat-cool*prof*(p%cs2-cs2cool)/cs2cool
        if (headtt) print*,'calc_heat_cool: cs20,cs2cool=',cs20,cs2cool
        heat=heat-cool*(p%cs2-(cs20-prof*cs2cool))/cs2cool
      endif
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
          prof = spread(exp(-0.5*((ztop-z(n))/wcool)**2), 1, l2-l1+1)
        case ('step')
          prof = step(spread(z(n),1,nx),z2,wcool)
        case ('cubic_step')
          prof = cubic_step(spread(z(n),1,nx),z2,wcool)
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
          if (idiag_dcoolz/=0) call xysum_mn_name_z(heat,idiag_dcoolz)
        endif
      endif
!
!  Add spatially uniform heating.
!
      if (heat_uniform/=0.) heat=heat+heat_uniform
!
!  Add cooling with constant time-scale to TTref_cool.
!
      if (tau_cool/=0.0) &
          heat=heat-p%rho*p%cp*gamma_inv*(p%TT-TTref_cool)/tau_cool
      if (tau_cool2/=0.0) &
          heat = heat - (p%cs2-cs2cool)/cs2cool/p%rho1/tau_cool2
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
!  Heating/cooling related diagnostics.
!
      if (ldiagnos) then
        if (idiag_heatm/=0) call sum_mn_name(heat,idiag_heatm)
      endif
!
    endsubroutine calc_heat_cool
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
          slices%ready=.true.
!
!  Pressure.
!
        case ('pp')
          do m=m1,m2; do n=n1,n2
            call eoscalc(ilnrho_ss,f(ix_loc,m,n,ilnrho),f(ix_loc,m,n,iss),pp=tmpval)
            slices%yz(m-m1+1,n-n1+1)=tmpval
          enddo; enddo
          do l=l1,l2; do n=n1,n2
            call eoscalc(ilnrho_ss,f(l,iy_loc,n,ilnrho),f(l,iy_loc,n,iss),pp=tmpval)
            slices%xz(l-l1+1,n-n1+1)=tmpval
          enddo; enddo
          do l=l1,l2; do m=m1,m2
            call eoscalc(ilnrho_ss,f(l,m,iz_loc,ilnrho),f(l,m,iz_loc,iss),pp=tmpval)
            slices%xy(l-l1+1,m-m1+1)=tmpval
            call eoscalc(ilnrho_ss,f(l,m,iz2_loc,ilnrho),f(l,m,iz2_loc,iss),pp=tmpval)
            slices%xy2(l-l1+1,m-m1+1)=tmpval
          enddo; enddo
          slices%ready=.true.
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
        hcond = hcond0
      endif
!
    endsubroutine heatcond
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
      integer :: iname,inamez,inamey,inamex,inamexy,inamexz
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtc=0; idiag_ethm=0; idiag_ethdivum=0; idiag_ssm=0; idiag_ss2m=0
        idiag_eem=0; idiag_ppm=0; idiag_csm=0; idiag_pdivum=0; idiag_heatm=0
        idiag_ugradpm=0; idiag_ethtot=0; idiag_dtchi=0 
        idiag_fradbot=0; idiag_fradtop=0; idiag_TTtop=0
        idiag_yHmax=0; idiag_yHm=0; idiag_TTmax=0; idiag_TTmin=0; idiag_TTm=0
        idiag_fconvm=0; idiag_fconvz=0; idiag_dcoolz=0; idiag_fradz=0
        idiag_fturbz=0; idiag_ppmx=0; idiag_ppmy=0; idiag_ppmz=0
        idiag_ssmx=0; idiag_ssmy=0; idiag_ssmz=0 
        idiag_TTmx=0; idiag_TTmy=0; idiag_TTmz=0; idiag_TTmxy=0; idiag_TTmxz=0
        idiag_uxTTmz=0; idiag_uyTTmz=0; idiag_uzTTmz=0
        idiag_ssmxy=0; idiag_ssmxz=0; idiag_fradz_Kprof=0; idiag_uxTTmxy=0
        idiag_fturbxy=0; idiag_fradxy_Kprof=0
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
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ssmx',idiag_ssmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ppmx',idiag_ppmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'TTmx',idiag_TTmx)
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
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ssmz',idiag_ssmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TTmz',idiag_TTmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ppmz',idiag_ppmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uxTTmz',idiag_uxTTmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uyTTmz',idiag_uyTTmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzTTmz',idiag_uzTTmz)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'TTmxy',idiag_TTmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'ssmxy',idiag_ssmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'uxTTmxy',idiag_uxTTmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'fturbxy',idiag_fturbxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'fradxy_Kprof',idiag_fradxy_Kprof)
      enddo
!
!  Check for those quantities for which we want y-averages.
!
      do inamexz=1,nnamexz
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'TTmxz',idiag_TTmxz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'ssmxz',idiag_ssmxz)
      enddo
!
!  Write column where which entropy variable is stored.
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
        write(3,*) 'i_fturbz=',idiag_fturbz
        write(3,*) 'i_fturbxy=',idiag_fturbxy
        write(3,*) 'i_fconvz=',idiag_fconvz
        write(3,*) 'i_dcoolz=',idiag_dcoolz
        write(3,*) 'i_fradz=',idiag_fradz
        write(3,*) 'i_fradz_Kprof=',idiag_fradz_Kprof
        write(3,*) 'i_fradxy_Kprof=',idiag_fradxy_Kprof
        write(3,*) 'i_ssmz=',idiag_ssmz
        write(3,*) 'i_TTmz=',idiag_TTmz
        write(3,*) 'i_uxTTmz=',idiag_uxTTmz
        write(3,*) 'i_uyTTmz=',idiag_uyTTmz
        write(3,*) 'i_uzTTmz=',idiag_uzTTmz
        write(3,*) 'nname=',nname
        write(3,*) 'iss=',iss
        write(3,*) 'i_yHmax=',idiag_yHmax
        write(3,*) 'i_yHm=',idiag_yHm
        write(3,*) 'i_TTmax=',idiag_TTmax
        write(3,*) 'i_TTmin=',idiag_TTmin
        write(3,*) 'i_TTm=',idiag_TTm
        write(3,*) 'i_TTp=',idiag_TTp
        write(3,*) 'i_TTmxy=',idiag_TTmxy
        write(3,*) 'i_TTmxz=',idiag_TTmxz
        write(3,*) 'i_uxTTmxy=',idiag_uxTTmxy
        write(3,*) 'i_ssmxy=',idiag_ssmxy
        write(3,*) 'i_ssmxz=',idiag_ssmxz
      endif
!
    endsubroutine rprint_entropy
!***********************************************************************
endmodule Entropy
