! $Id$
!
!  This module can replace the entropy module by using lnT or T (with
!  ltemperature_nolog=.true.) as dependent variable. For a perfect gas
!  with constant coefficients (no ionization) we have:
!  (1-1/gamma) * cp*T = cs02 * exp( (gamma-1)*ln(rho/rho0)-gamma*s/cp )
!
!  Note that to use lnTT as thermal variable, you may rather want to use
!  entropy.f90 with pretend_lnTT=.true. As of March 2007, entropy.f90
!  has way more options and features than temperature_idealgas.f90.
!
!  At a later point we may want to rename the module Entropy into Energy
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lentropy = .false.
! CPARAM logical, parameter :: ltemperature = .true.
! CPARAM logical, parameter :: lthermal_energy = .false.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED Ma2; uglnTT; ugTT; fpres(3); tcond;
!
!***************************************************************
module Entropy
!
  use Cdata
  use Cparam
  use EquationOfState, only: mpoly0, mpoly1, mpoly2
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'entropy.h'
!
  real :: radius_lnTT=0.1, ampl_lnTT=0.0, widthlnTT=2*epsi
  real :: lnTT_const=0.0, TT_const=1.0
  real :: Kgperp=0.0, Kgpara=0.0
  real :: chi=impossible
  real :: zbot=0.0, ztop=0.0
  real :: center1_x=0.0, center1_y=0.0, center1_z=0.0
  real :: r_bcz=0.0, chi_shock=0.0, chi_hyper3=0.0, chi_hyper3_mesh=5.0
  real :: Tbump=0.0, Kmin=0.0, Kmax=0.0, hole_slope=0.0, hole_width=0.0
  real :: hcond0=impossible, hcond1=1.0, hcond2=1.0, Fbot=impossible
  real :: luminosity=0.0, wheat=0.1, rcool=0.0, wcool=0.1, cool=0.0
  integer, parameter :: nheatc_max=3
  logical, pointer :: lpressuregradient_gas
  logical :: ladvection_temperature=.true.
  logical :: lupw_lnTT=.false., lcalc_heat_cool=.false., lheatc_hyper3=.false.
  logical :: lheatc_Kconst=.false., lheatc_Kprof=.false., lheatc_Karctan=.false.
  logical :: lheatc_tensordiffusion=.false., lheatc_hyper3_mesh=.false.
  logical :: lheatc_chiconst=.false., lheatc_chiconst_accurate=.false.
  logical :: lfreeze_lnTTint=.false., lfreeze_lnTText=.false.
  logical :: lhcond_global=.false.
  logical :: lheatc_shock=.false., lheatc_hyper3_polar=.false.
  logical :: lviscosity_heat=.true.
  integer :: iglobal_hcond=0
  integer :: iglobal_glhc=0
  logical :: linitial_log=.false.
  character (len=labellen), dimension(nheatc_max) :: iheatcond='nothing'
  character (len=labellen), dimension(ninit) :: initlnTT='nothing'
  character (len=5) :: iinit_str
  logical :: lADI_mixed=.false.
!
!  Input parameters.
!
  namelist /entropy_init_pars/ &
      initlnTT, radius_lnTT, ampl_lnTT, widthlnTT, &
      lnTT_const, TT_const, center1_x, center1_y, &
      center1_z, mpoly0, mpoly1, mpoly2, r_bcz, Fbot, Tbump, Kmin, Kmax, &
      hole_slope, hole_width, ltemperature_nolog, linitial_log, hcond0, &
      luminosity, wheat
!
!  Run parameters.
!
  namelist /entropy_run_pars/ &
      lupw_lnTT, ladvection_temperature, &
      chi, iheatcond, chi_hyper3_mesh, &
      lheatc_chiconst_accurate, hcond0, lcalc_heat_cool, lfreeze_lnTTint, &
      lfreeze_lnTText, widthlnTT, mpoly0, mpoly1, mpoly2, lhcond_global, &
      lviscosity_heat, chi_hyper3, chi_shock, Fbot, Tbump, Kmin, Kmax, &
      hole_slope, hole_width, Kgpara, Kgperp, lADI_mixed, &
      rcool, wcool, cool
!
!  Diagnostic variables for print.in
! (needs to be consistent with reset list below)
!
  integer :: idiag_TTmax=0    ! DIAG_DOC: $\max (T)$
  integer :: idiag_gTmax=0    ! DIAG_DOC: $\max (|\nabla T|)$
  integer :: idiag_TTmin=0    ! DIAG_DOC: $\min (T)$
  integer :: idiag_TTm=0      ! DIAG_DOC: $\left< T \right>$
  integer :: idiag_fradtop=0  ! DIAG_DOC: $<-K{dT\over dz}>_{\text{top}}$
                              ! DIAG_DOC: \quad(radiative flux at the top)
  integer :: idiag_yHmax=0, idiag_yHmin=0, idiag_yHm=0
  integer :: idiag_ethm=0     ! DIAG_DOC: $\left< e_{\text{th}}\right> =
                              ! DIAG_DOC:  \left< c_v \rho T \right> $
                              ! DIAG_DOC: \quad(mean thermal energy)
  integer :: idiag_eem=0      ! DIAG_DOC: $\left< e \right> =
                              ! DIAG_DOC:  \left< c_v T \right>$
                              ! DIAG_DOC: \quad(mean internal energy)
  integer :: idiag_ssm=0, idiag_thcool=0
  integer :: idiag_ppm=0, idiag_csm=0
  integer :: idiag_dtc=0        ! DIAG_DOC: $\delta t/[c_{\delta t}\,\delta_x
                                ! DIAG_DOC:   /\max c_{\rm s}]$
                                ! DIAG_DOC:   \quad(time step relative to
                                ! DIAG_DOC:   acoustic time step;
                                ! DIAG_DOC:   see \S~\ref{time-step})
  integer :: idiag_dtchi=0      ! DIAG_DOC: $\delta t / [c_{\delta t,{\rm v}}\,
                                ! DIAG_DOC:   \delta x^2/\chi_{\rm max}]$
                                ! DIAG_DOC:   \quad(time step relative to time
                                ! DIAG_DOC:   step based on heat conductivity;
                                ! DIAG_DOC:   see \S~\ref{time-step})
!
! xy averaged diagnostics given in xyaver.in
!
  integer :: idiag_ppmz=0       ! XYAVG_DOC: $\left<p\right>_{xy}$
  integer :: idiag_ppuzmz=0     ! XYAVG_DOC:
  integer :: idiag_TTmz=0       ! XYAVG_DOC: $\left<T\right>_{xy}$
  integer :: idiag_ethmz=0      ! XYAVG_DOC: $\left< e_{\text{th}}
                                ! XYAVG_DOC: \right>_{xy}$
  integer :: idiag_ethuxmz=0    ! XYAVG_DOC:
  integer :: idiag_ethuymz=0    ! XYAVG_DOC:
  integer :: idiag_ethuzmz=0    ! XYAVG_DOC:
  integer :: idiag_fpresxmz=0   ! XYAVG_DOC: $\left<(\nabla p)_x\right>_{xy}$
  integer :: idiag_fpresymz=0   ! XYAVG_DOC: $\left<(\nabla p)_y\right>_{xy}$
  integer :: idiag_fpreszmz=0   ! XYAVG_DOC: $\left<(\nabla p)_z\right>_{xy}$
!
! xz averaged diagnostics given in xzaver.in
!
  integer :: idiag_ppmy=0       ! XZAVG_DOC: $\left<p\right>_{xz}$
  integer :: idiag_TTmy=0       ! XZAVG_DOC: $\left<T\right>_{xz}$
!
! yz averaged diagnostics given in yzaver.in
!
  integer :: idiag_ppmx=0       ! YZAVG_DOC: $\left<p\right>_{yz}$
  integer :: idiag_TTmx=0       ! YZAVG_DOC: $\left<T\right>_{yz}$
  integer :: idiag_ethuxmx=0    ! YZAVG_DOC:
!
! variables for slices given in video.in
!
  real, dimension(nx,nz) :: pp_xz
  real, dimension(ny,nz) :: pp_yz
  real, dimension(nx,ny) :: pp_xy,pp_xy2,pp_xy3,pp_xy4
!
! y averaged diagnostics given in yaver.in
!
  integer :: idiag_TTmxz=0      ! YAVG_DOC: $\left<T\right>_{y}$
!
! z averaged diagnostics given in zaver.in
!
  integer :: idiag_TTmxy=0      ! ZAVG_DOC: $\left<T\right>_{z}$
!
  contains
!***********************************************************************
    subroutine register_entropy()
!
!  Initialise variables which should know that we solve an entropy
!  equation: ilnTT, etc; increase nvar accordingly.
!
!  6-nov-01/wolf: coded
!
      use FArrayManager, only: farray_register_pde
      use SharedVariables, only: get_shared_variable
!
      integer :: ierr
!
      call farray_register_pde('lnTT',ilnTT)
!
!  logical variable lpressuregradient_gas shared with hydro modules
!
      call get_shared_variable('lpressuregradient_gas',lpressuregradient_gas,ierr)
      if (ierr/=0) call fatal_error('register_entropy','lpressuregradient_gas')
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Writing files for use with IDL.
!
      if (lroot) then
        if (maux == 0) then
           if (nvar < mvar) write(4,*) ',lnTT $'
           if (nvar == mvar) write(4,*) ',lnTT'
        else
           write(4,*) ',lnTT $'
        endif
        write(15,*) 'lnTT = fltarr(mx,my,mz)*one'
      endif
!
    endsubroutine register_entropy
!***********************************************************************
    subroutine initialize_entropy(f,lstarting)
!
!  Called by run.f90 after reading parameters, but before the time loop.
!
!  21-jul-2002/wolf: coded
!
      use FArrayManager, only: farray_register_global
      use Gravity, only: g0, gravz, compute_gravity_star
      use EquationOfState, only : cs2bot, cs2top, gamma, gamma_m1, &
                                  select_eos_variable
      use Sub, only: step,der_step
      use SharedVariables, only: put_shared_variable
      use Mpicomm, only: stop_it
!
      logical :: lstarting
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: hcond, dhcond
      logical :: lnothing
      integer :: i, ierr
      real, dimension(5) :: hole_params
      real :: star_cte
!
!  Set iTT equal to ilnTT if we are considering non-logarithmic temperature.
!
      if (ltemperature_nolog) iTT=ilnTT
!
      if (.not. leos) then
         call fatal_error('initialize_entropy', &
             'EOS=noeos but temperature_idealgas requires an EQUATION OF STATE')
      endif
!
      if (ltemperature_nolog) then
        call select_eos_variable('TT',iTT)
      else
        call select_eos_variable('lnTT',ilnTT)
      endif
!
!  Freeze temperature.
!
      if (lfreeze_lnTTint) lfreeze_varint(ilnTT)=.true.
      if (lfreeze_lnTText) lfreeze_varext(ilnTT)=.true.
!
!  Check whether we want heat conduction.
!
      lheatc_Kconst= .false.
      lheatc_Kprof= .false.
      lheatc_Karctan= .false.
      lheatc_tensordiffusion=.false.
      lheatc_chiconst = .false.
!
!  Initialize thermal diffusion.
!
      lheatc_shock=.false.
      lheatc_hyper3=.false.
      lheatc_hyper3_mesh=.false.
      lheatc_hyper3_polar=.false.
!
!  initialize lnothing. It is needed to prevent multiple output.
!
      lnothing = .false.
!
!  Different choices of heat conduction (if any).
!
      do i=1,nheatc_max
      select case (iheatcond(i))
        case ('K-const')
          lheatc_Kconst=.true.
          if (lroot) call information('initialize_entropy', &
          ' heat conduction: K=cst --> gamma*K/rho/TT/cp*div(T*grad lnTT)')
        case ('K-profile')
          lheatc_Kprof=.true.
!  11-Aug-2008/dintrans: better somewhere else?
          hcond1=(mpoly1+1.)/(mpoly0+1.)
          hcond2=(mpoly2+1.)/(mpoly0+1.)
          Fbot=-gamma/(gamma-1.)*hcond0*g0/(mpoly0+1.)
          if (lroot) &
              call information('initialize_entropy',' heat conduction: K=K(r)')
        case ('K-arctan')
          lheatc_Karctan=.true.
          if (.not. ltemperature_nolog) &
            call fatal_error('initialize_entropy', &
              'K-arctan only valid for TT')
          if (lADI_mixed .and. .not. lADI) &
            call fatal_error('initialize_entropy', &
              'K-arctan with lADI_mixed=T while lADI=F?')
          if (lroot) call information('initialize_entropy', &
              'heat conduction: arctan profile')
        case ('chi-const')
          lheatc_chiconst=.true.
          if (lroot) call information('initialize_entropy', &
              ' heat conduction: constant chi')
        case ('chi-hyper3')
          lheatc_hyper3=.true.
          if (lroot) call information('initialize_entropy','hyper conductivity')
        case ('hyper3_mesh','hyper3-mesh')
          lheatc_hyper3_mesh=.true.
          if (lroot) call information('initialize_entropy','hyper mesh conductivity')
        case ('hyper3_cyl','hyper3-cyl','hyper3_sph','hyper3-sph')
          lheatc_hyper3_polar=.true.
          if (lroot) call information('initialize_entropy', &
              'hyper conductivity: polar coords')
        case ('shock','chi-shock')
          lheatc_shock=.true.
          if (lroot) call information('initialize_entropy','shock conductivity')
        case ('tensor-diffusion')
          lheatc_tensordiffusion=.true.
          if (lroot) print*, 'heat conduction: tensor diffusion'
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
!  Compute and store hcond and dhcond if hcond_global=.true.
!
      if (lhcond_global) then
        call farray_register_global("hcond",iglobal_hcond)
        call farray_register_global("glhc",iglobal_glhc)
        do n=n1,n2
        do m=m1,m2
          hcond = 1. + (hcond1-1.)*step(x(l1:l2),r_bcz,-widthlnTT)
          hcond = hcond0*hcond
          dhcond = hcond0*(hcond1-1.)*der_step(x(l1:l2),r_bcz,-widthlnTT)
          f(l1:l2,m,n,iglobal_hcond)=hcond
          f(l1:l2,m,n,iglobal_glhc)=dhcond
        enddo
        enddo
      endif
!
      if (initlnTT(1)=='gaussian') then
!
!  Needed when one only works with temperature_idealgas to check the radiative
!  diffusion term, i.e. one solves d(TT)/dt=gamma*chi*del2(TT) with bcz='cT'
!  (all other modules are down).
!
        cs2bot=gamma_m1*f(l1,4,n1,ilnTT)
        cs2top=gamma_m1*f(l1,4,n2,ilnTT)
      endif
!
!  Some tricks regarding Fbot and hcond0 when bcz1='c1' (constant flux).
!
      if (bcz1(ilnTT)=='c1' .and. lrun) then
        if (Fbot==impossible .and. hcond0 /= impossible) then
          Fbot=-gamma/gamma_m1*hcond0*gravz/(mpoly0+1.0)
          if (lroot) print*, &
              'initialize_entropy: Calculated Fbot = ', Fbot
        endif
        if (hcond0==impossible .and. Fbot /= impossible) then
          hcond0=-Fbot*gamma_m1/gamma*(mpoly0+1.0)/gravz
          if (lroot) print*, &
              'initialize_entropy: Calculated hcond0 = ', hcond0
        endif
        if (Fbot==impossible .and. hcond0==impossible) &
          call fatal_error('temperature_idealgas',  &
              'Both Fbot and hcond0 are unknown')
      endif
!
      if (initlnTT(1)=='star_heat') then
        if (lroot) print*,'star_heat: compute the gravity profile'
        ! compute the gravity profile
        star_cte=(mpoly0+1.)/hcond0*gamma_m1/gamma
        call compute_gravity_star(f, wheat, luminosity, star_cte)
        if (rcool==0.) rcool=r_ext
      endif
!
!  Now we share several variables.
!
      call put_shared_variable('hcond0', hcond0, ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting hcond0')
      call put_shared_variable('Fbot', Fbot, ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting Fbot')
      call put_shared_variable('lADI_mixed', lADI_mixed, ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting lADI_mixed')
      call put_shared_variable('lviscosity_heat',lviscosity_heat,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting lviscosity_heat')
!
!  Share the 4 parameters of the radiative conductivity hole (kappa-mechanism
!  problem).
!
      hole_params=(/Tbump,Kmax,Kmin,hole_slope,hole_width/)
      call put_shared_variable('hole_params',hole_params,ierr)
      if (ierr/=0) call fatal_error('initialize_entropy', &
          'there was a problem when putting the hole_params array')
!
!  A word of warning...
!
      if (lheatc_Kconst .and. hcond0==0.0) then
        call warning('initialize_entropy', 'hcond0 is zero!')
      endif
      if (lheatc_Kprof .and. hcond0==0.0) then
        call warning('initialize_entropy', 'hcond0 is zero!')
      endif
      if (lheatc_chiconst .and. chi==0.0) then
        call warning('initialize_entropy','chi is zero!')
      endif
      if (lrun) then
        if (lheatc_hyper3 .and. chi_hyper3==0.0) &
            call fatal_error('initialize_entropy', &
            'Conductivity coefficient chi_hyper3 is zero!')
        if (lheatc_shock .and. chi_shock==0.0) &
            call fatal_error('initialize_entropy', &
            'Conductivity coefficient chi_shock is zero!')
      endif
      if (iheatcond(1)=='nothing') then
        if (hcond0/=impossible) call warning('initialize_entropy', &
            'No heat conduction, but hcond0/=0')
        if (chi/=impossible) call warning('initialize_entropy', &
            'No heat conduction, but chi/=0')
      endif
      if (lADI_mixed .and. iheatcond(1) /= 'K-arctan') then
        call stop_it("temperature_idealgas: "//&
          "lADI_mixed=T while iheatcond /= K-arctan?")
      endif
!
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_entropy
!***********************************************************************
    subroutine init_ss(f)
!
!  Initialise lnTT or TT; called from start.f90.
!
!  13-dec-2002/axel+tobi: adapted from init_ss
!
!  initialise entropy; called from start.f90
!  07-nov-2001/wolf: coded
!  24-nov-2002/tony: renamed for consistancy (i.e. init_[variable name])
!
      use General,  only: chn
      use Sub,      only: blob
      use InitialCondition, only: initial_condition_ss
      use EquationOfState, only: gamma, gamma_m1, cs2bot, cs2top, cs20, &
                                 lnrho0, get_cp1
      use Gravity, only: gravz
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
!
      integer :: j
      logical :: lnothing=.true.
      real :: haut, Rgas, cp1, Ttop, alpha, beta, expo
!
      do j=1,ninit
!
        if (initlnTT(j)/='nothing') then
!
          lnothing=.false.
!
          call chn(j,iinit_str)
!
!  Select between various initial conditions.
!
          select case (initlnTT(j))
          case ('zero', '0'); f(:,:,:,ilnTT) = 0.
!
          case ('const_lnTT'); f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+lnTT_const
!
          case ('const_TT')
            if (ltemperature_nolog) then
              f(:,:,:,iTT)=f(:,:,:,iTT)+TT_const
            else
              f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+log(TT_const)
            endif
!
          case ('single_polytrope'); call single_polytrope(f)
!
          case ('piecew-poly'); call piecew_poly(f)
!
          case ('gaussian')
            do n=n1,n2
              f(l1:l2,4,n,ilnTT)=exp(-(x(l1:l2)/radius_lnTT)**2)* &
                  exp(-((z(n)-0.5)/radius_lnTT)**2)
            enddo
            cs2bot=gamma_m1*f(l1,4,n1,ilnTT)
            cs2top=gamma_m1*f(l1,4,n2,ilnTT)
!
          case ('rad_equil')
            call rad_equil(f)
!
          case ('blob_hs')
            if (lroot) print*, 'init_lnTT: hydrostatic blob with ', &
                radius_lnTT, ampl_lnTT, center1_x, center1_y, center1_z
            call blob(ampl_lnTT,f,ilnTT,radius_lnTT, &
                center1_x, center1_y,center1_z)
            call blob(-ampl_lnTT,f,ilnrho,radius_lnTT, &
                center1_x,center1_y,center1_z)
!
          case ('blob')
            if (lroot) print*, 'init_lnTT: blob ', &
                radius_lnTT, ampl_lnTT, center1_x, center1_y, center1_z
            call blob(ampl_lnTT,f,ilnTT,radius_lnTT, &
                center1_x,center1_y,center1_z)
!
          case ('isothermal')
            if (lroot) print*, 'init_lnTT: isothermal atmosphere'
            if (ltemperature_nolog) then
              f(:,:,:,iTT)  =cs20/gamma_m1
            else
              f(:,:,:,ilnTT)=log(cs20/gamma_m1)
            endif
            haut=-cs20/gamma/gravz
            do n=n1,n2
              f(:,:,n,ilnrho)=lnrho0+(1.-z(n))/haut
            enddo
!
          case ('hydro_rad')
            if (lroot) print*, 'init_lnTT: hydrostatic+radiative equilibria'
            if (Fbot==impossible .or. hcond0==impossible) &
                call stop_it("initialize_lnTT: Fbot or hcond0 not initialized")
            call get_cp1(cp1)
            Rgas=(1.-1./gamma)/cp1
            Ttop=cs20/gamma_m1
            beta=-Fbot/hcond0
            alpha=Ttop-beta
            expo=-gravz/beta/Rgas
            do n=n1,n2
              if (ltemperature_nolog) then
                f(:,:,n,iTT)  =beta*z(n)+alpha
              else
                f(:,:,n,ilnTT)=log(beta*z(n)+alpha)
              endif
              f(:,:,n,ilnrho)=lnrho0+ &
                  (1.+expo)*log((1.+alpha/beta)/(z(n)+alpha/beta))
            enddo
!
          case ('star_heat')
            call star_heat(f)
!
          case default
!
!  Catch unknown values.
!
            write(unit=errormsg,fmt=*) 'No such value for initss(' &
                           //trim(iinit_str)//'): ',trim(initlnTT(j))
            call fatal_error('init_ss',errormsg)
!
          endselect
!
          if (lroot) print*,'init_ss: initss(' &
              //trim(iinit_str)//') = ',trim(initlnTT(j))
        endif
      enddo
!
!  Interface for user's own initial condition.
!
      if (linitial_condition) call initial_condition_ss(f)
!
      if (lnothing.and.lroot) print*,'init_ss: nothing'
!
      if (ltemperature_nolog.and.linitial_log) f(:,:,:,iTT)=exp(f(:,:,:,ilnTT))
!
    endsubroutine init_ss
!***********************************************************************
    subroutine pencil_criteria_entropy()
!
!  All pencils that the Entropy module depends on are specified here.
!
!  20-11-04/anders: coded
!
      if (dvid/=0.0) lpenc_video(i_pp)=.true.
!
      if (ldt) lpenc_requested(i_cs2)=.true.
!
      if (lpressuregradient_gas) lpenc_requested(i_fpres)=.true.
!
      if (lviscosity.and.lviscosity_heat) then
        lpenc_requested(i_cv1)=.true.
        lpenc_requested(i_visc_heat)=.true.
        if (.not.ltemperature_nolog) &
            lpenc_requested(i_TT1)=.true.
      endif
!
      if (ldensity) then
        lpenc_requested(i_divu)=.true.
        if (ltemperature_nolog) lpenc_requested(i_TT)=.true.
      endif
!
      if (lcalc_heat_cool) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_cv1)=.true.
        if (lgravr) lpenc_requested(i_r_mn)=.true.
      endif
!
      if (lheatc_chiconst) then
        if (ltemperature_nolog) then
          lpenc_requested(i_del2TT)=.true.
          lpenc_requested(i_gTT)=.true.
        else
          lpenc_requested(i_del2lnTT)=.true.
          lpenc_requested(i_glnTT)=.true.
        endif
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_cp1)=.true.
      endif
!
      if (lheatc_Kconst) then
        if (ldensity) lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_cp1)=.true.
        if (ltemperature_nolog) then
          lpenc_requested(i_del2TT)=.true.
        else
          lpenc_requested(i_glnTT)=.true.
          lpenc_requested(i_del2lnTT)=.true.
        endif
      endif
!
      if (lheatc_Kprof) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_cp1)=.true.
        if (ltemperature_nolog) then
          lpenc_requested(i_gTT)=.true.
          lpenc_requested(i_del2TT)=.true.
        else
          lpenc_requested(i_glnTT)=.true.
          lpenc_requested(i_del2lnTT)=.true.
        endif
        if (lgravz) lpenc_requested(i_z_mn)=.true.
        if (lgravr) lpenc_requested(i_r_mn)=.true.
      endif
!
      if (lheatc_Karctan) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_gTT)=.true.
        if (.not. lADI_mixed) lpenc_requested(i_del2TT)=.true.
      endif
!
      if (lheatc_shock) then
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_del2lnrho)=.true.
        lpenc_requested(i_gshock)=.true.
        if (ltemperature_nolog) then
          lpenc_requested(i_gTT)=.true.
          lpenc_requested(i_del2TT)=.true.
        else
          lpenc_requested(i_glnTT)=.true.
          lpenc_requested(i_del2lnTT)=.true.
        endif
      endif
!
      if (lheatc_tensordiffusion) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_bij)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_cp1)=.true.
      endif
!
      if (lheatc_hyper3) then
        if (ltemperature_nolog) then
          lpenc_requested(i_del6TT)=.true.
        else
          lpenc_requested(i_del6lnTT)=.true.
        endif
      endif
!
      if (lheatc_hyper3_mesh) lpenc_requested(i_TT1)=.true.
!
      if (ladvection_temperature) then
        if (ltemperature_nolog) then
          lpenc_requested(i_ugTT)=.true.
        else
          lpenc_requested(i_uglnTT)=.true.
        endif
      endif
!
!  Diagnostic pencils.
!
      if (idiag_TTmax/=0) lpenc_diagnos(i_TT)  =.true.
      if (idiag_gTmax/=0) then
         lpenc_diagnos(i_glnTT) =.true.
         lpenc_diagnos(i_TT) =.true.
      endif
      if (idiag_TTmin/=0) lpenc_diagnos(i_TT)  =.true.
      if (idiag_TTm/=0)   lpenc_diagnos(i_TT)  =.true.
      if (idiag_fradtop/=0) then
        lpenc_diagnos(i_TT) =.true.  ! for hcond computation
        lpenc_diagnos(i_glnTT) =.true.
      endif
      if (idiag_yHmax/=0) lpenc_diagnos(i_yH)  =.true.
      if (idiag_yHmin/=0) lpenc_diagnos(i_yH)  =.true.
      if (idiag_yHm/=0)   lpenc_diagnos(i_yH)  =.true.
      if (idiag_ethm/=0.or.idiag_ethmz/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_ee)  =.true.
      endif
      if (idiag_ethuxmz/=0.or.idiag_ethuymz/=0.or.idiag_ethuzmz/=0.or.&
          idiag_ethuxmx/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_ee) =.true.
        lpenc_diagnos(i_uu) =.true.
      endif
      if (idiag_ssm/=0)    lpenc_diagnos(i_ss)  =.true.
      if (idiag_dtchi/=0)  lpenc_diagnos(i_cs2)=.true.
      if (idiag_csm/=0)    lpenc_diagnos(i_cs2)=.true.
      if (idiag_eem/=0)    lpenc_diagnos(i_ee) =.true.
      if (idiag_ppm/=0 .or. idiag_ppmx/=0 .or. idiag_ppmy/=0 .or. &
          idiag_ppmz/=0 .or. idiag_ppuzmz/=0) lpenc_diagnos(i_pp) =.true.
      if (idiag_thcool/=0) lpenc_diagnos(i_rho)=.true.
      if (idiag_TTmx/=0 .or. idiag_TTmy/=0 .or. idiag_TTmz/=0) &
          lpenc_diagnos(i_TT)=.true.
      if (idiag_dtchi/=0) then
        lpenc_diagnos(i_rho1)=.true.
        lpenc_diagnos(i_cv1) =.true.
      endif
      if (idiag_fpresxmz/=0 .or. idiag_fpresymz/=0 .or. &
          idiag_fpreszmz/=0) lpenc_requested(i_fpres)=.true.
!
      if (idiag_TTmxy/=0 .or. idiag_TTmxz/=0) lpenc_diagnos2d(i_TT)=.true.
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
      if (lpencil_in(i_Ma2)) then
        lpencil_in(i_u2)=.true.
        lpencil_in(i_cs2)=.true.
      endif
      if (lpencil_in(i_uglnTT)) then
        lpencil_in(i_glnTT)=.true.
        lpencil_in(i_uu)=.true.
      endif
      if (lpencil_in(i_ugTT)) then
        lpencil_in(i_gTT)  =.true.
        lpencil_in(i_uu)=.true.
      endif
!      
      if (lpencil_in(i_fpres)) then
        lpencil_in(i_cs2)=.true.
        lpencil_in(i_glnrho)=.true.
        lpencil_in(i_glnTT)=.true.
        lpencil_in(i_glnmumol)=.true.
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
      use Sub, only: u_dot_grad
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      type (pencil_case), intent (inout) :: p
      integer :: j
! Ma2
      if (lpencil(i_Ma2)) p%Ma2=p%u2/p%cs2
! uglnTT
      if (lpencil(i_uglnTT)) &
          call u_dot_grad(f,ilnTT,p%glnTT,p%uu,p%uglnTT,UPWIND=lupw_lnTT)
! ugTT
      if (lpencil(i_ugTT)) &
          call u_dot_grad(f,ilnTT,p%gTT,p%uu,p%ugTT,UPWIND=lupw_lnTT)
! fpres
      if (lpencil(i_fpres)) then
        do j=1,3
          p%fpres(:,j)=-gamma_inv*p%cs2* &
              (p%glnrho(:,j)+p%glnTT(:,j)-p%glnmumol(:,j))
        enddo
      endif
! tcond
      if (lpencil(i_tcond)) then
        if (lheatc_chiconst) then
          p%tcond=chi*p%rho/p%cp1
        elseif (lheatc_Kconst) then
          p%tcond=hcond0
        else
          call fatal_error('calc_pencils_entropy',  &
              'This heatcond is not implemented to work with lpencil(i_cond)!')
        endif
      endif
!
    endsubroutine calc_pencils_entropy
!***********************************************************************
    subroutine dss_dt(f,df,p)
!
!  Calculate right hand side of temperature equation.
!
!  lnTT version: DlnTT/Dt = -gamma_m1*divu + gamma*cp1*rho1*TT1*RHS
!    TT version:   DTT/Dt = -gamma_m1*TT*divu + gamma*cp1*rho1*RHS
!
!  13-dec-02/axel+tobi: adapted from entropy
!
      use Deriv, only: der6
      use Diagnostics
      use EquationOfState, only: gamma_m1
      use ImplicitPhysics, only: heatcond_TT
      use Special, only: special_calc_entropy
      use Sub, only: dot2,identify_bcs
      use Viscosity, only: calc_viscous_heat
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: Hmax=0.0, hcond, thdiff, tmp
      real :: fradtop
      integer :: j
!
      intent(inout) :: f,p
      intent(out) :: df
!
! Initialization of thdiff in the declaration the
! variable gets the SAVE attribute,
! so in the next call thdiff is not initialized anymore.
!
      thdiff = 0.0
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*, 'SOLVE dlnTT_dt'
      if (headtt) call identify_bcs('lnTT',ilnTT)
      if (headtt) print*, 'dss_dt: lnTT,cs2=', p%lnTT(1), p%cs2(1)
!
!  Sound speed squared.
!
      if (headtt) print*, 'dss_dt: cs20=', p%cs2(1)
!
!  ``cs2/dx^2'' for timestep
!
      if (lfirst.and.ldt) advec_cs2=p%cs2*dxyz_2
      if (headtt.or.ldebug) print*,'dss_dt: max(advec_cs2) =',maxval(advec_cs2)
!
!  Add pressure gradient term in momentum equation.
!
      if (lpressuregradient_gas) &
          df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + p%fpres
!
!  Advection term and PdV-work.
!
      if (ladvection_temperature) then
        if (ltemperature_nolog) then
          df(l1:l2,m,n,iTT)   = df(l1:l2,m,n,iTT)   - p%ugTT
        else
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%uglnTT
        endif
      endif
!
!  Need to add left-hand-side of the continuity equation (see manual).
!  Check this:
!
      if (ldensity) then
        if (ltemperature_nolog) then
          df(l1:l2,m,n,iTT)   = df(l1:l2,m,n,iTT)   - gamma_m1*p%TT*p%divu
        else
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - gamma_m1*p%divu
        endif
      endif
!
!  Calculate viscous contribution to temperature.
!
      if (lviscosity.and.lviscosity_heat) call calc_viscous_heat(f,df,p,Hmax)
!
!  Various heating conduction contributions.
!
      if (lcalc_heat_cool)  call calc_heat_cool(df,p)
!
!  Thermal conduction
!
      if (lheatc_chiconst) call calc_heatcond_constchi(df,p)
      if (lheatc_Kconst)   call calc_heatcond_constK(df,p)
      if (lheatc_Kprof)    call calc_heatcond(f,df,p)
      if (lheatc_Karctan)  call calc_heatcond_arctan(df,p)
      if (lheatc_tensordiffusion) call calc_heatcond_tensor(df,p)
!
!  Hyper diffusion.
!
      if (lheatc_hyper3) then
        if (ltemperature_nolog) then
          thdiff=thdiff+chi_hyper3*p%del6TT
        else
          thdiff=thdiff+chi_hyper3*p%del6lnTT
        endif
        if (lfirst.and.ldt) diffus_chi3=diffus_chi3+chi_hyper3*dxyz_6
        if (headtt) print*,'dss_dt: chi_hyper3=', chi_hyper3
      endif
!
      if (lheatc_hyper3_mesh) then
        do j=1,3
          call der6(f,ilnTT,tmp,j,IGNOREDX=.true.)
          if (.not.ltemperature_nolog) tmp=tmp*p%TT1
          thdiff = thdiff + chi_hyper3_mesh*pi5_1/60.*tmp*dline_1(:,j)
        enddo
        if (lfirst.and.ldt) &
            advec_hypermesh_ss=chi_hyper3_mesh*pi5_1*sqrt(dxyz_2)
        if (headtt) print*,'dss_dt: chi_hyper3_mesh=', chi_hyper3_mesh
      endif
!
      if (lheatc_hyper3_polar) then
        do j=1,3
          call der6(f,ilnTT,tmp,j,IGNOREDX=.true.)
          if (.not.ltemperature_nolog) tmp=tmp*p%TT1
          thdiff = thdiff + chi_hyper3*pi4_1/60.*tmp*dline_1(:,j)**2
        enddo
        if (lfirst.and.ldt) &
             diffus_chi3=diffus_chi3+chi_hyper3*pi4_1*dxyz_2
        if (headtt) print*,'dss_dt: chi_hyper3=', chi_hyper3
      endif
!
!  Shock diffusion.
!
      if (lheatc_shock) call calc_heatcond_shock(df,p)
!
!  Entry possibility for "personal" entries.
!  In that case you'd need to provide your own "special" routine.
!
      if (lspecial) call special_calc_entropy(f,df,p)
!
!  Add thermal diffusion to temperature equation.
!
      if (ltemperature_nolog) then
        df(l1:l2,m,n,iTT)   = df(l1:l2,m,n,iTT)   + thdiff
      else
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + thdiff
      endif
!
!  Information on the timescales.
!
      if (lfirst.and.ldt.and.(headtt.or.ldebug)) then
        print*, 'dss_dt: max(diffus_chi ) =', maxval(diffus_chi)
        print*, 'dss_dt: max(diffus_chi3) =', maxval(diffus_chi3)
      endif
!
!  Calculate temperature related diagnostics.
!
      if (ldiagnos) then
        if (idiag_TTm/=0)   call sum_mn_name(p%TT,idiag_TTm)
        if (idiag_TTmax/=0) call max_mn_name(p%TT,idiag_TTmax)
        if (idiag_TTmin/=0) call max_mn_name(-p%TT,idiag_TTmin,lneg=.true.)
        if (idiag_ssm/=0)   call sum_mn_name(p%ss,idiag_ssm)
        if (idiag_eem/=0)   call sum_mn_name(p%ee,idiag_eem)
        if (idiag_ppm/=0)   call sum_mn_name(p%pp,idiag_ppm)
        if (idiag_ethm/=0)  call sum_mn_name(p%rho*p%ee,idiag_ethm)
        if (idiag_csm/=0)   call sum_mn_name(p%cs2,idiag_csm,lsqrt=.true.)
        if (idiag_dtc/=0) then
          call max_mn_name(sqrt(advec_cs2)/cdt,idiag_dtc,l_dt=.true.)
        endif
        if (idiag_gTmax/=0) then
          call dot2(p%glnTT,tmp)
          call max_mn_name(p%TT*sqrt(tmp),idiag_gTmax)
        endif
        if (idiag_fradtop/=0.and.n==n2) then
          call heatcond_TT(p%TT,hcond)
          fradtop=sum(-hcond*p%glnTT(:,3))/nx
          call save_name(fradtop, idiag_fradtop)
        endif
      endif
!
!  1-D averages.
!
      if (l1davgfirst) then
        if (idiag_ppmx/=0)  call yzsum_mn_name_x(p%pp,idiag_ppmx)
        if (idiag_ppmy/=0)  call xzsum_mn_name_y(p%pp,idiag_ppmy)
        if (idiag_ppmz/=0)  call xysum_mn_name_z(p%pp,idiag_ppmz)
        if (idiag_TTmx/=0)  call yzsum_mn_name_x(p%TT,idiag_TTmx)
        if (idiag_TTmy/=0)  call xzsum_mn_name_y(p%TT,idiag_TTmy)
        if (idiag_TTmz/=0)  call xysum_mn_name_z(p%TT,idiag_TTmz)
        if (idiag_ppuzmz/=0)  call xysum_mn_name_z(p%pp*p%uu(:,3),idiag_ppuzmz)
        if (idiag_ethmz/=0)   call xysum_mn_name_z(p%rho*p%ee,idiag_ethmz)
        if (idiag_ethuxmx/=0) call yzsum_mn_name_x(p%rho*p%ee*p%uu(:,1), &
            idiag_ethuxmx)
        if (idiag_ethuxmz/=0) call xysum_mn_name_z(p%rho*p%ee*p%uu(:,1), &
            idiag_ethuxmz)
        if (idiag_ethuymz/=0) call xysum_mn_name_z(p%rho*p%ee*p%uu(:,2), &
            idiag_ethuymz)
        if (idiag_ethuzmz/=0) call xysum_mn_name_z(p%rho*p%ee*p%uu(:,3), &
            idiag_ethuzmz)
        if (idiag_fpresxmz/=0) call xysum_mn_name_z(p%fpres(:,1),idiag_fpresxmz)
        if (idiag_fpresymz/=0) call xysum_mn_name_z(p%fpres(:,2),idiag_fpresymz)
        if (idiag_fpreszmz/=0) call xysum_mn_name_z(p%fpres(:,3),idiag_fpreszmz)
      endif
!
!  2-D averages.
!
      if (l2davgfirst) then
        if (idiag_TTmxy/=0) call zsum_mn_name_xy(p%TT,idiag_TTmxy)
        if (idiag_TTmxz/=0) call ysum_mn_name_xz(p%TT,idiag_TTmxz)
      endif
!
      if (lvideo.and.lfirst) then
        pp_yz(m-m1+1,n-n1+1)=p%pp(ix_loc-l1+1)
        if (m==iy_loc)  pp_xz(:,n-n1+1)=p%pp
        if (n==iz_loc)  pp_xy(:,m-m1+1)=p%pp
        if (n==iz2_loc) pp_xy2(:,m-m1+1)=p%pp
        if (n==iz3_loc) pp_xy3(:,m-m1+1)=p%pp
        if (n==iz4_loc) pp_xy4(:,m-m1+1)=p%pp
      endif
!
    endsubroutine dss_dt
!***********************************************************************
    subroutine calc_lentropy_pars(f)
!
!  Dummy routine.
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lentropy_pars
!***********************************************************************
    subroutine calc_heatcond_shock(df,p)
!
!  Add shock diffusion to the energy equation,
!
!    De/Dt = ... + div(K*grad(T)) = ... + div(cv*T*Xi*grad(T)) ,
!
!  where e=rho*cv*T and Xi is a regular shock diffusion coefficient.
!
!  01-aug-08/wlad: adapted from entropy
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: thdiff, g2, gshockgTT, gshockglnTT
!
      intent(in) :: p
      intent(out) :: df
!
      if (headtt) print*, 'calc_heatcond_shock: chi_shock=', chi_shock
!
!  Shock energy diffusivity.
!
      if (ltemperature_nolog) then
        call dot(p%gshock,p%gTT,gshockgTT)
        call dot(p%glnrho,p%gTT,g2)
        thdiff=chi_shock*(p%shock*(p%del2TT+g2)+gshockgTT)
      else
        call dot(p%gshock,p%glnTT,gshockglnTT)
        call dot(p%glnrho+p%glnTT,p%glnTT,g2)
        thdiff=chi_shock*(p%shock*(p%del2lnTT+g2)+gshockglnTT)
      endif
!
      df(l1:l2,m,n,ilntt) = df(l1:l2,m,n,ilntt) + thdiff
!
      if (headtt) print*,'calc_heatcond_shock: added thdiff'
!
      if (lfirst.and.ldt) then
        if (leos_idealgas) then
          diffus_chi=diffus_chi+(gamma*chi_shock*p%shock)*dxyz_2
        else
          diffus_chi=diffus_chi+(chi_shock*p%shock)*dxyz_2
        endif
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_shock
!***********************************************************************
    subroutine rad_equil(f)
!
!  Compute the radiative and hydrostatic equilibria for a given radiative
!  profile defined in heatcond_TT.
!
!  16-may-07/gastine+dintrans: coded
!
      use Gravity, only: gravz
      use EquationOfState, only: lnrho0,cs20,cs2top,cs2bot,gamma, &
                                 gamma_m1,eoscalc,ilnrho_TT
      use ImplicitPhysics, only: heatcond_TT
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nzgrid) :: temp,lnrho
      real :: hcond, dtemp, dlnrho, ss
      integer :: i,n,iz
!
      if (.not. ltemperature_nolog) &
        call fatal_error('temperature_idealgas',  &
                         'rad_equil not implemented for lnTT')
      if (lroot) print*,'init_ss: rad_equil for kappa-mechanism pb'
!
!  Integrate from top to bottom: z(n2) --> z(n1).
!
      temp(nzgrid)=cs20/gamma_m1
      lnrho(nzgrid)=lnrho0
!
!  Calculate the n2-1 gridpoint thanks to a 1st order forward Euler scheme.
!
      call heatcond_TT(temp(nzgrid), hcond)
      dtemp=Fbot/hcond
      temp(nzgrid-1)=temp(nzgrid)+dz*dtemp
      dlnrho=(-gamma/gamma_m1*gravz-dtemp)/temp(nzgrid)
      lnrho(nzgrid-1)=lnrho(nzgrid)+dz*dlnrho
!
!  Now we use a 2nd order centered scheme for the other gridpoints.
!
      do i=nzgrid-1,2,-1
        call heatcond_TT(temp(i), hcond)
        dtemp=Fbot/hcond
        temp(i-1)=temp(i+1)+2.*dz*dtemp
        dlnrho=(-gamma/gamma_m1*gravz-dtemp)/temp(i)
        lnrho(i-1)=lnrho(i+1)+2.*dz*dlnrho
      enddo
!
      do n=1,nz
        iz=ipz*nz+n
        f(:,:,nghost+n,ilnTT)=temp(iz)
        f(:,:,nghost+n,ilnrho)=lnrho(iz)
      enddo
!
!  Initialize cs2bot by taking into account the new bottom value of temperature
!  Note: cs2top=cs20 already defined in eos_idealgas.
!
      cs2bot=gamma_m1*temp(1)
      print*,'cs2top, cs2bot=', cs2top, cs2bot
!
      if (lroot) then
        print*,'--> write the initial setup in data/proc0/setup.dat'
        open(unit=11,file=trim(directory)//'/setup.dat')
        write(11,'(5a14)') 'z','rho','temp','ss','hcond'
        do i=n2,n1,-1
          call eoscalc(ilnrho_TT,lnrho(i),temp(i),ss=ss)
          call heatcond_TT(temp(i), hcond)
          write(11,'(5e14.5)') z(i),exp(lnrho(i)),temp(i),ss,hcond
        enddo
        close(11)
      endif
!
    endsubroutine rad_equil
!***********************************************************************
    subroutine calc_heat_cool(df,p)
!
      use Diagnostics, only: sum_lim_mn_name
      use EquationOfState, only: cs20
      use Sub, only: step
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: tau, cooling, kappa, a1, a3, prof, heat
      real :: a2, kappa0, kappa0_cgs
!
!  Initialize
!
      intent(in) :: p
      intent(out) :: df
!
      if (headtt) print*,'enter calc_heat_cool', rcool, wcool, cool, cs20
!
      if (lgravr) then
        ! 2-D heating/cooling profiles
        prof = exp(-0.5*(p%r_mn/wheat)**2) * (2*pi*wheat**2)**(-1.)
        heat = luminosity*prof
        prof = step(p%r_mn,rcool,wcool)
        heat = heat - cool*prof*(p%cs2-cs20)/cs20
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + p%cv1*p%TT1*heat
      else
        kappa0_cgs=2e-6  !cm2/g
        kappa0=kappa0_cgs*unit_density*unit_length
        kappa=kappa0*p%TT**2
!
!  Optical Depth tau=kappa*rho*H.
!  If we are using 2D, the pencil value p%rho is actually sigma, the column
!  density, sigma=rho*2*H
!
        if (nzgrid==1) then
          tau = .5*kappa*p%rho
        else
          call fatal_error("calc_heat_cool", &
              "opacity not yet implemented for 3D")
        endif
!
!  Analytical gray description of Hubeny (1990)
!  a1 is the optically thick contribution,
!  a3 the optically thin one.
!
        a1=0.375*tau ; a2=0.433013 ; a3=0.25/tau
!
!  Cooling for energy: 2*sigmaSB*p%TT**4/(a1+a2+a3)
!
        cooling = 2*sigmaSB*p%rho1*p%TT**4/(a1+a2+a3)
!
!  This cooling has dimension of energy over time.
!
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%cv1*p%TT1*cooling
      endif
!
      if (ldiagnos) then
         !cooling power - energy radiated away (luminosity)
         if (idiag_thcool/=0) call sum_lim_mn_name(cooling*p%rho,idiag_thcool,p)
      endif
!
    endsubroutine calc_heat_cool
!***********************************************************************
    subroutine calc_heatcond_constchi(df,p)
!
!  Calculate the radiative diffusion term for chi=cte:
!  lnTT version: cp*chi*Div(rho*T*glnTT)/(rho*cv*TT)
!           = gamma*chi*(g2.glnTT+g2lnTT) where g2=glnrho+glnTT
!    TT version: cp*chi*Div(rho*gTT)/(rho*cv)
!           = gamma*chi*(g2.gTT+g2TT) where g2=glnrho
!
!  01-mar-07/dintrans: adapted from temperature_ionization
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: g2
!
      intent(in) :: p
      intent(inout) :: df
!
      if (ltemperature_nolog) then
        call dot(p%glnrho,p%gTT,g2)
        g2=g2+p%del2TT
      else
        call dot(p%glnTT+p%glnrho,p%glnTT,g2)
        g2=g2+p%del2lnTT
      endif
!
!  Add heat conduction to RHS of temperature equation.
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + gamma*chi*g2
!
!  Check maximum diffusion from thermal diffusion.
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chi*dxyz_2
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_constchi
!***********************************************************************
    subroutine calc_heatcond_constK(df,p)
!
!  Calculate the radiative diffusion term for K=cte:
!
!  lnTT version: gamma*K/rho/TT/cp*div(T*grad lnTT)
!                =gamma*K/rho/cp*(gradlnTT.gradlnTT + del2ln TT)
!    TT version: gamma*K/rho/cp*del2(TT)=gamma*chi*del2(TT)
!
!  Note: if ldensity=.false. then rho=1 and chi=K/cp
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot
!
      real, dimension(mx,my,mz,mvar) :: df
      type (pencil_case)  :: p
      real, dimension(nx) :: g2, chix
!
      intent(in) :: p
      intent(inout) :: df
!
!  Add heat conduction to RHS of temperature equation.
!
      if (ldensity) then
        chix=p%rho1*hcond0*p%cp1
      else
        chix=hcond0*p%cp1
      endif
!
      if (ltemperature_nolog) then
        df(l1:l2,m,n,iTT)   = df(l1:l2,m,n,iTT)   + gamma*chix*p%del2TT
      else
        call dot(p%glnTT,p%glnTT,g2)
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + gamma*chix*(g2 + p%del2lnTT)
      endif
!
!  Check maximum diffusion from thermal diffusion.
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
    subroutine calc_heatcond_arctan(df,p)
!
!  Radiative diffusion with an arctan profile for the conductivity
!
!  Calculate gamma/(rho*cp)*div(K * grad TT)=
!      gamma*K/(rho*cp)*(grad LnK.grad TT + del2 TT)
!
!  16-may-07/gastine+dintrans: coded
!  01-mar-10/dintrans: introduced a mixed version with the ADI scheme that only
!  computes *during the explicit step* the term
!  gamma/(rho*cp)*grad(K).grad(TT) with grad(K)=dK/dT*grad(TT),
!  this term being less restrictive for the explicit timestep
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot, multsv
      use ImplicitPhysics, only: heatcond_TT
!
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension (nx)   :: hcond, dhcond, g1, chix
      real, dimension (nx,3) :: gLnhcond=0.
      type (pencil_case)     :: p
!
      intent(in) :: p
      intent(inout) :: df
!
      call heatcond_TT(p%TT, hcond, dhcond)
!  must specify the new bottom value of hcond for the 'c1' BC
!     if (n == n1) hcond0=hcond(1)
      if (lADI_mixed) then
        call dot(p%gTT, p%gTT, g1)
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + gamma*p%rho1*p%cp1*dhcond*g1
        chix=0.
      else
!  grad LnK=grad_T Ln K.grad(TT)
        dhcond=dhcond/hcond
        call multsv(dhcond, p%gTT, gLnhcond)
        call dot(gLnhcond, p%gTT, g1)
        chix=p%rho1*p%cp1*hcond
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + gamma*chix*(g1+p%del2TT)
      endif
!
!
!  Check maximum diffusion from thermal diffusion.
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chix*dxyz_2
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_arctan
!***********************************************************************
    subroutine calc_heatcond(f,df,p)
!
!  Calculate the radiative diffusion term for a variable K:
!    ivar=lnTT --> 1/(rho*cv*T)*div(K*grad TT)
!    ivar=TT   --> 1/(rho*cv)*div(K*grad TT)
!
!  12-Mar-07/dintrans: coded
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot, step, der_step
      use Gravity, only: z1, z2
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension(nx) :: g2, hcond, dhcond, chix
      real, dimension (nx,3) :: glhc=0.,glnThcond
!
      intent(in) :: f,p
      intent(inout) :: df
!
      if (lhcond_global) then
        hcond=f(l1:l2,m,n,iglobal_hcond)
        glhc(:,1)=f(l1:l2,m,n,iglobal_glhc)
      else
        if (lgravz) then
          hcond = 1. + (hcond1-1.)*step(p%z_mn,z1,-widthlnTT) &
                     + (hcond2-1.)*step(p%z_mn,z2,widthlnTT)
          hcond = hcond0*hcond
          glhc(:,3) = (hcond1-1.)*der_step(p%z_mn,z1,-widthlnTT) &
                    + (hcond2-1.)*der_step(p%z_mn,z2,widthlnTT)
          glhc(:,3) = hcond0*glhc(:,3)
        elseif (lcylindrical_coords) then
          hcond = 1. + (hcond1-1.)*step(rcyl_mn,r_bcz,-widthlnTT)
          hcond = hcond0*hcond
          glhc(:,1) = hcond0*(hcond1-1.)*der_step(rcyl_mn,r_bcz,-widthlnTT)
        elseif (lgravr) then
          hcond = 1. + (hcond1-1.)*step(p%r_mn,r_bcz,-widthlnTT) &
                     + (hcond2-1.)*step(p%r_mn,r_ext,widthlnTT)
          hcond = hcond0*hcond
          dhcond=(hcond1-1.)*der_step(p%r_mn,r_bcz,-widthlnTT) &
                 + (hcond2-1.)*der_step(p%r_mn,r_ext,widthlnTT)
          dhcond=hcond0*dhcond
          glhc(:,1) = x(l1:l2)/p%r_mn*dhcond
          glhc(:,2) = y(m)/p%r_mn*dhcond
          glhc(:,3) = z(n)/p%r_mn*dhcond
        endif
      endif
!
      if (ltemperature_nolog) then
        glnThcond = glhc/spread(hcond,2,3)              ! grad ln(hcond)
        call dot(p%gTT,glnThcond,g2)
        g2 = g2 + p%del2TT
      else
        glnThcond = p%glnTT + glhc/spread(hcond,2,3)    ! grad ln(T*hcond)
        call dot(p%glnTT,glnThcond,g2)
        g2 = g2 + p%del2lnTT
      endif
!
!  Add heat conduction to RHS of temperature equation.
!
      chix=p%rho1*hcond*p%cp1
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + gamma*chix*g2
!
!  Check maximum diffusion from thermal diffusion.
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chix*dxyz_2
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond
!***********************************************************************
    subroutine calc_heatcond_tensor(df,p)
!
!  Calculates heat conduction parallel and perpendicular (isotropic)
!  to magnetic field lines.
!
!  25-aug-09/bing: moved from dss_dt to here
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot,dot2,tensor_diffusion_coef
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: cosbgT,gT2,b2
      real, dimension (nx) :: vKpara,vKperp,rhs
!
      vKpara(:) = Kgpara
      vKperp(:) = Kgperp
!
      call tensor_diffusion_coef(p%glnTT,p%hlnTT,p%bij,p%bb, &
          vKperp,vKpara,rhs,llog=.true.)
!
      df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT)+rhs*p%rho1*gamma*p%cp1
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
        diffus_chi=diffus_chi+cosbgT*gamma*Kgpara*p%rho1*p%cp1*dxyz_2
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_tensor
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
    subroutine rprint_entropy(lreset,lwrite)
!
!  Reads and registers print parameters relevant to entropy.
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Diagnostics, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname, inamex, inamey, inamez, inamexy, inamexz
      logical :: lwr
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_TTmax=0; idiag_TTmin=0; idiag_TTm=0; idiag_fradtop=0
        idiag_yHmax=0; idiag_yHmin=0; idiag_yHm=0; idiag_gTmax=0
        idiag_ethm=0; idiag_ssm=0; idiag_thcool=0
        idiag_dtchi=0; idiag_dtc=0
        idiag_eem=0; idiag_ppm=0; idiag_csm=0
        idiag_ppmx=0; idiag_ppmy=0; idiag_ppmz=0; idiag_ppuzmz=0
        idiag_TTmx=0; idiag_TTmy=0; idiag_TTmz=0; idiag_ethuxmx=0
        idiag_ethmz=0; idiag_ethuxmz=0; idiag_ethuymz=0; idiag_ethuzmz=0
        idiag_TTmxy=0; idiag_TTmxz=0
        idiag_fpresxmz=0; idiag_fpresymz=0; idiag_fpreszmz=0;
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'TTmax',idiag_TTmax)
        call parse_name(iname,cname(iname),cform(iname),'gTmax',idiag_gTmax)
        call parse_name(iname,cname(iname),cform(iname),'TTmin',idiag_TTmin)
        call parse_name(iname,cname(iname),cform(iname),'TTm',idiag_TTm)
        call parse_name(iname,cname(iname),cform(iname),'fradtop',idiag_fradtop)
        call parse_name(iname,cname(iname),cform(iname),'ethm',idiag_ethm)
        call parse_name(iname,cname(iname),cform(iname),'ssm',idiag_ssm)
        call parse_name(iname,cname(iname),cform(iname),'dtchi',idiag_dtchi)
        call parse_name(iname,cname(iname),cform(iname),'dtc',idiag_dtc)
        call parse_name(iname,cname(iname),cform(iname),'eem',idiag_eem)
        call parse_name(iname,cname(iname),cform(iname),'ppm',idiag_ppm)
        call parse_name(iname,cname(iname),cform(iname),'csm',idiag_csm)
        call parse_name(iname,cname(iname),cform(iname),'thcool',idiag_thcool)
      enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ppmx',idiag_ppmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'TTmx',idiag_TTmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ethuxmx', &
            idiag_ethuxmx)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'ppmy',idiag_ppmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'TTmy',idiag_TTmy)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ppmz',idiag_ppmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TTmz',idiag_TTmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ppuzmz', &
            idiag_ppuzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ethmz', &
            idiag_ethmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ethuxmz', &
            idiag_ethuxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ethuymz', &
            idiag_ethuymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ethuzmz', &
            idiag_ethuzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fpresxmz', &
            idiag_fpresxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fpresymz', &
            idiag_fpresymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fpreszmz', &
            idiag_fpreszmz)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'TTmxy', &
            idiag_TTmxy)
      enddo
!
!  Check for those quantities for which we want y-averages.
!
      do inamexz=1,nnamexz
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'TTmxz', &
            idiag_TTmxz)
      enddo
!
!  Write column where which variable is stored.
!
      if (lwr) then
        write(3,*) 'nname=',nname
        write(3,*) 'ilnTT=',ilnTT
        write(3,*) 'iyH=',iyH
        write(3,*) 'iss=',iss
      endif
!
    endsubroutine rprint_entropy
!***********************************************************************
    subroutine get_slices_entropy(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Temperature.
!
        case ('TT')
          if (iTT>0) then
            slices%yz =f(ix_loc,m1:m2,n1:n2,iTT)
            slices%xz =f(l1:l2,iy_loc,n1:n2,iTT)
            slices%xy =f(l1:l2,m1:m2,iz_loc,iTT)
            slices%xy2=f(l1:l2,m1:m2,iz2_loc,iTT)
            if (lwrite_slice_xy3) slices%xy3=f(l1:l2,m1:m2,iz3_loc,iTT)
            if (lwrite_slice_xy4) slices%xy4=f(l1:l2,m1:m2,iz4_loc,iTT)
          else
            slices%yz =exp(f(ix_loc,m1:m2,n1:n2,ilnTT))
            slices%xz =exp(f(l1:l2,iy_loc,n1:n2,ilnTT))
            slices%xy =exp(f(l1:l2,m1:m2,iz_loc,ilnTT))
            slices%xy2=exp(f(l1:l2,m1:m2,iz2_loc,ilnTT))
            if (lwrite_slice_xy3) slices%xy3=exp(f(l1:l2,m1:m2,iz3_loc,ilnTT))
            if (lwrite_slice_xy4) slices%xy4=exp(f(l1:l2,m1:m2,iz4_loc,ilnTT))
          endif
          slices%ready=.true.
!  lnTT
        case ('lnTT')
          slices%yz =f(ix_loc,m1:m2,n1:n2,ilnTT)
          slices%xz =f(l1:l2,iy_loc,n1:n2,ilnTT)
          slices%xy =f(l1:l2,m1:m2,iz_loc,ilnTT)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,ilnTT)
          if (lwrite_slice_xy3) slices%xy3=f(l1:l2,m1:m2,iz3_loc,ilnTT)
          if (lwrite_slice_xy4) slices%xy4=f(l1:l2,m1:m2,iz4_loc,ilnTT)
          slices%ready=.true.
!  Pressure
        case ('pp')
          slices%yz =pp_yz
          slices%xz =pp_xz
          slices%xy =pp_xy
          slices%xy2=pp_xy2
          if (lwrite_slice_xy3) slices%xy3=pp_xy3
          if (lwrite_slice_xy4) slices%xy4=pp_xy4
          slices%ready=.true.
!
      endselect
!
    endsubroutine get_slices_entropy
!***********************************************************************
    subroutine single_polytrope(f)
!
!  04-aug-07/dintrans: a single polytrope with index mpoly0
!
      use Gravity, only: gravz
      use EquationOfState, only: cs20, lnrho0, gamma, gamma_m1, get_cp1, &
                                 cs2bot, cs2top
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real :: beta, zbot, ztop, cp1, T0, temp
!
!  beta is the (negative) temperature gradient
!  beta = -(g/cp) /[(1-1/gamma)*(m+1)]
!  gamma*(Rgas/mu)T0 = cs2(ad) = cp*T0*gamma_m1,
!  so T0 = cs20*cp1/gamma_m1
!
      call get_cp1(cp1)
      beta=-cp1*gravz/(mpoly0+1.)*gamma/gamma_m1
      ztop=xyz0(3)+Lxyz(3)
      zbot=xyz0(3)
      T0=cs20*cp1/gamma_m1
      print*, 'polytrope: mpoly0, beta, T0=', mpoly0, beta, T0
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        temp=T0+beta*(ztop-z(n))
        if (ltemperature_nolog) then
          f(:,m,n,iTT)  =temp
        else
          f(:,m,n,ilnTT)=log(temp)
        endif
        f(:,m,n,ilnrho)=lnrho0+mpoly0*log(temp/T0)
      enddo
      cs2bot=gamma_m1*(T0+beta*(ztop-zbot))
      cs2top=cs20
!
    endsubroutine single_polytrope
!***********************************************************************
    subroutine piecew_poly(f)
!
!  Computes piecewise polytropic and hydrostatic atmosphere.
!  Adapted from single_polytrope.
!
!  19-jan-10/bing: coded
!  The layout is the same than in entropy.f90:
!  ------------ ztop Ttop
!     mpoly2
!  ------------ z2   T2, lnrho2
!     mpoly0
!  ------------ z1   T1, lnrho1
!     mpoly1
!  ------------ zbot
!
      use Gravity, only: gravz, z1, z2
      use EquationOfState, only: cs2top, cs2bot, gamma, gamma_m1, lnrho0, &
                                 get_cp1
!
      real, dimension(mx,my,mz,mfarray) :: f
      real :: Ttop, T1, T2, beta0, beta1, beta2, cp1, temp
      real :: lnrhotop, lnrho1, lnrho2, ztop
      integer :: i
!
      call get_cp1(cp1)
!
!  Top boundary values.
!
      Ttop=cs2top*cp1/gamma_m1
      lnrhotop = lnrho0
      ztop=xyz0(3)+Lxyz(3)
!
!  Temperature gradients.
!
      beta0 =-cp1*gravz/(mpoly0+1.)*gamma/gamma_m1
      beta1 =-cp1*gravz/(mpoly1+1.)*gamma/gamma_m1
      beta2 =-cp1*gravz/(mpoly2+1.)*gamma/gamma_m1
!
      T2 = Ttop + beta2*(ztop-z2)
      T1 = T2   + beta0*(z2-z1)
!
      lnrho2 = lnrhotop+mpoly2*log(T2/Ttop)
      lnrho1 = lnrho2  +mpoly0*log(T1/T2)
!
      do  i=n2,n1,-1
        if (z(i) >= z2)                 temp = Ttop + beta2*(ztop-z(i))
        if (z(i) < z2 .and. z(i) >= z1) temp = T2   + beta0*(z2-z(i))
        if (z(i) < z1)                  temp = T1   + beta1*(z1-z(i))
!
        if (ltemperature_nolog) then
          f(:,:,i,iTT)  =temp
        else
          f(:,:,i,ilnTT)=log(temp)
        endif
!
        if (z(i) >= z2) f(:,:,i,ilnrho)=lnrhotop+mpoly2*log(temp/Ttop)
        if (z(i) < z2 .and. z(i) >= z1 ) &
            f(:,:,i,ilnrho)=lnrho2+mpoly0*log(temp/T2)
        if (z(i) < z1) f(:,:,i,ilnrho)=lnrho1+mpoly1*log(temp/T1)
      enddo
!
! one also needs to refresh cs2bot in case of a 'cT' BC for the temperature
!
      cs2bot=gamma_m1*(T1 + beta1*(z1-xyz0(3)))
!
    endsubroutine piecew_poly
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
    subroutine star_heat(f)
!
!  Initialize entropy for two superposed polytropes with a central heating
!
!  04-fev-2011/dintrans: coded
!
      use EquationOfState, only: rho0, lnrho0, get_soundspeed, eoscalc, &
                                 ilnrho_TT, gamma, gamma_m1
      use Sub, only: step, interp1, erfunc
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer, parameter   :: nr=100
      integer              :: i,l,iter
      real, dimension (nr) :: r, lnrho, temp, lumi, g, hcond
      real                 :: u,r_mn,lnrho_r,temp_r,cs2,ss
      real                 :: rhotop, rbot,rt_old,rt_new,rhobot
      real                 :: rb_old,rb_new,crit,r_max
!
!  Define the radial grid r=[0,r_max], luminosity and gravity
!
      r_max=sqrt(xyz1(1)**2+xyz1(2)**2+xyz1(3)**2)
      r(1)=0. ; lumi(1)=0. ; g(1)=0.
      do i=2,nr
        r(i)=r_max*float(i-1)/(nr-1)
        u=r(i)/sqrt(2.)/wheat
        lumi(i)=luminosity*(1.-exp(-u**2))
        g(i)=-lumi(i)/(2.*pi*r(i))*(mpoly0+1.)/hcond0*gamma_m1/gamma
      enddo
!
      hcond1=(mpoly1+1.)/(mpoly0+1.)
      hcond2=(mpoly2+1.)/(mpoly0+1.)
      hcond = 1. + (hcond1-1.)*step(r,r_bcz,-widthlnTT) &
                 + (hcond2-1.)*step(r,r_ext,widthlnTT)
      hcond = hcond0*hcond
!
      rbot=rho0
      rt_old=0.01*rbot
      rt_new=0.012*rbot
      rhotop=rt_old
      call strat_heat(nr, r, lumi, g, hcond, temp, lnrho, rhotop, rhobot)
      print*, 'find rhobot=', rhobot
      rb_old=rhobot
!
      rhotop=rt_new
      call strat_heat(nr, r, lumi, g, hcond, temp, lnrho, rhotop, rhobot)
      print*, 'find rhobot=', rhobot
      rb_new=rhobot
!
      do 10 iter=1,10
        rhotop=rt_old+(rt_new-rt_old)/(rb_new-rb_old)*(rbot-rb_old)
!
        crit=abs(rhotop/rt_new-1.)
        if (crit<=1e-4) goto 20
        call strat_heat(nr, r, lumi, g, hcond, temp, lnrho, rhotop, rhobot)
!
!  Update new estimates.
!
        rt_old=rt_new
        rb_old=rb_new
        rt_new=rhotop
        rb_new=rhobot
 10 continue
 20 print*,'- iteration completed: rhotop,crit=',rhotop,crit
!
!  One needs to refresh rho0 and lnrho0 because the density top value
!  has changed --> important for the future EOS calculations (ss, ...)
!
      lnrho0=lnrho(nr)
      rho0=exp(lnrho0)
      print*,'new values for lnrho0 and rho0:', lnrho0, rho0
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        do l=l1,l2
          r_mn=sqrt(x(l)**2+y(m)**2+z(n)**2)
          lnrho_r=interp1(r,lnrho,nr,r_mn)
          temp_r=interp1(r,temp,nr,r_mn)
          f(l,m,n,ilnrho)=lnrho_r
          f(l,m,n,ilnTT)=temp_r
        enddo
      enddo
!
      if (lroot) then
        print*,'--> writing initial setup to data/proc0/setup.dat'
        open(unit=11,file=trim(directory)//'/setup.dat')
        write(11,'(a1,a5,6a12)') '#','r','rho','ss','cs2','grav', &
          'lumi','hcond'
        do i=nr,1,-1
          u=r(i)/sqrt(2.)/wheat
          call get_soundspeed(log(temp(i)),cs2)
          call eoscalc(ilnrho_TT,lnrho(i),temp(i),ss=ss)
          write(11,'(f6.3,6e12.3)') r(i),exp(lnrho(i)),ss,cs2,g(i), &
            lumi(i), hcond(i)
        enddo
        close(11)
      endif
!
    endsubroutine star_heat
!***********************************************************************
    subroutine strat_heat(nr,r,lumi,g,hcond,temp,lnrho,rhotop,rhobot)
!
      use EquationOfState, only: gamma, gamma_m1, cs20
      use Sub, only: interp1
!
      integer              :: nr, i
      real, dimension (nr) :: r, lnrho, temp, lumi, g, hcond
      real                 :: dr,dtemp,dlnrho
      real                 :: rhotop,rhobot,lnrhobot
!
      temp(nr)=cs20/gamma_m1 ; lnrho(nr)=alog(rhotop)
      dr=r(2)
      do i=nr-1,1,-1
        if (r(i+1) > r_ext) then
          ! Isothermal exterior: mpoly2 but force T=cte
          dtemp=0.
          dlnrho=-gamma*g(i+1)/cs20
        elseif (r(i+1) > r_bcz) then
          ! Convection zone: mpoly0
! adiabatic stratification
!          dtemp=-g(i+1)
!          dlnrho=3./2.*dtemp/temp(i+1)
          dtemp=lumi(i+1)/(2.*pi*r(i+1))/hcond(i+1)
          dlnrho=mpoly0*dtemp/temp(i+1)
        else
          ! Radiative zone: mpoly1
          dtemp=lumi(i+1)/(2.*pi*r(i+1))/hcond(i+1)
          dlnrho=mpoly1*dtemp/temp(i+1)
        endif
        temp(i)=temp(i+1)+dtemp*dr
        lnrho(i)=lnrho(i+1)+dlnrho*dr
      enddo
!
      lnrhobot=interp1(r,lnrho,nr,r_ext)
      rhobot=exp(lnrhobot)
!
    endsubroutine strat_heat
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
