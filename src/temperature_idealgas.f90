! $Id$
!
!  This module can replace the energy module by using lnT or T (with
!  ltemperature_nolog=.true.) as dependent variable. For a perfect gas
!  with constant coefficients (no ionization) we have:
!  (1-1/gamma) * cp*T = cs20 * exp( (gamma-1)*ln(rho/rho0)-gamma*s/cp )
!
!  Note that to use lnTT as thermal variable, you may rather want to use
!  energy.f90 with pretend_lnTT=.true. As of March 2007, entropy.f90
!  has way more options and features than temperature_idealgas.f90.
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
! PENCILS PROVIDED Ma2; uglnTT; ugTT; fpres(3); tcond; sglnTT(3); d2xlnTT; d2zlnTT
!
!***************************************************************
module Energy
!
! 12-may-12/MR: made ampl_lnTT a vector; added parameters for initialization
!               by mode to input pars
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'energy.h'
!
  real :: radius_lnTT=0.1, widthlnTT=2*epsi
  real, dimension (ninit) :: ampl_lnTT=0.0
  real :: lnTT_const=0.0, TT_const=1.0
  real :: Kgperp=0.0, Kgpara=0.0
  real :: chi=impossible, chi_jump=1., chi_z0=0.0, chi_zwidth=0.0
  real :: zbot=0.0, ztop=0.0
  real :: center1_x=0.0, center1_y=0.0, center1_z=0.0
  real :: r_bcz=0.0, chi_shock=0.0, chi_hyper3=0.0, chi_hyper3_mesh=5.0
  real :: Tbump=0.0, Kmin=0.0, Kmax=0.0, hole_slope=0.0, hole_width=0.0
  real, dimension(5) :: hole_params
  real, dimension(nz) :: zmask_temp, zmask_emiss
  real, dimension(nzgrid) :: zmask_temp_global
  real, dimension(2) :: temp_zaver_range=(/-max_real,max_real/)
  real :: mu=1.0
  real :: hcond0=impossible, hcond1=1.0, hcond2=1.0, Fbot=impossible,Ftop=impossible
  real :: luminosity=0.0, wheat=0.1, rcool=0.0, wcool=0.1, cool=0.0
  real :: beta_bouss=-1.0
  integer :: temp_zmask_count=1
  integer, parameter :: nheatc_max=3
  logical, pointer :: lpressuregradient_gas
  logical :: ladvection_temperature=.true.
  logical :: lupw_lnTT=.false., lcalc_heat_cool=.false., lheatc_hyper3=.false.
  logical :: lheatc_Kconst=.false., lheatc_Kprof=.false., lheatc_Karctan=.false.
  logical :: lheatc_tensordiffusion=.false., lheatc_hyper3_mesh=.false.
  logical :: lheatc_chiconst=.false., lheatc_chiconst_accurate=.false.
  logical :: lfreeze_lnTTint=.false., lfreeze_lnTText=.false.
  logical :: lhcond_global=.false., lheatc_chicubicstep=.false.
  logical :: lheatc_shock=.false., lheatc_hyper3_polar=.false.
  logical :: lheatc_Ktherm=.false.
  logical, target :: lheatc_kramers=.false.
  real, target :: hcond0_kramers=0.0, nkramers=0.0
  logical :: lviscosity_heat=.true.
  logical :: lcalc_TTmean=.false.
  integer :: iglobal_hcond=0
  integer :: iglobal_glhc=0
  logical :: lenergy_slope_limited=.false.
  logical :: linitial_log=.false.
  logical, pointer :: lreduced_sound_speed   !=>.false.
!  logical, pointer :: lscale_to_cs2top
  character (len=labellen), dimension(nheatc_max) :: iheatcond='nothing'
  character (len=labellen) :: borderss='nothing'
  character (len=labellen), dimension(ninit) :: initlnTT='nothing'
  complex :: coef_lnTT=0.
  character (len=intlen) :: iinit_str
  real    :: kx_lnTT=1.,ky_lnTT=1.,kz_lnTT=1.
  logical :: lADI_mixed=.false., lmultilayer=.false.
  real, pointer :: PrRa   ! preliminary
  real, pointer :: reduce_cs2
  real, target :: mpoly0=1.5, mpoly1=1.5, mpoly2=1.5
!
  real, dimension(nz) :: TTmz, gTTmz 
!
  real, dimension(3) :: gradTT0=(/0.0,0.0,0.0/)
!
!
!  Init parameters.
!
  namelist /entropy_init_pars/ &
      initlnTT, radius_lnTT, ampl_lnTT, widthlnTT, lnTT_const, TT_const, &
      center1_x, center1_y, center1_z, mpoly0, mpoly1, mpoly2, r_bcz, Fbot, &
      Tbump, Kmin, Kmax, hole_slope, hole_width, ltemperature_nolog, &
      linitial_log, hcond0, luminosity, wheat, coef_lnTT, kx_lnTT, ky_lnTT, kz_lnTT, &
      temp_zaver_range
!
!  Run parameters.
!
  namelist /entropy_run_pars/ &
      lupw_lnTT, ladvection_temperature, chi, iheatcond, chi_hyper3_mesh, &
      lheatc_chiconst_accurate, hcond0, lcalc_heat_cool, lfreeze_lnTTint, &
      lfreeze_lnTText, widthlnTT, mpoly0, mpoly1, mpoly2, lhcond_global, &
      lviscosity_heat, chi_hyper3, chi_shock, Fbot, Tbump, Kmin, Kmax, &
      hole_slope, hole_width, Kgpara, Kgperp, lADI_mixed, rcool, wcool, &
      cool, beta_bouss, borderss, lmultilayer, lcalc_TTmean, &
      temp_zaver_range, &
      gradTT0, w_sldchar_ene, chi_z0, chi_jump, chi_zwidth, &
      hcond0_kramers, nkramers
!
!  Diagnostic variables for print.in
! (needs to be consistent with reset list below)
!
  integer :: idiag_TTmax=0        ! DIAG_DOC: $\max (T)$
  integer :: idiag_gTmax=0        ! DIAG_DOC: $\max (|\nabla T|)$
  integer :: idiag_TTmin=0        ! DIAG_DOC: $\min (T)$
  integer :: idiag_TTm=0          ! DIAG_DOC: $\left< T \right>$
  integer :: idiag_TTzmask=0      ! DIAG_DOC: $\left< T \right>$ for
                                  ! DIAG_DOC: the temp_zaver_range
  integer :: idiag_TT2m=0         ! DIAG_DOC: $\left< T^2 \right>$
  integer :: idiag_TugTm=0        ! DIAG_DOC: $\left< T\uv\cdot\nabla T \right>$
  integer :: idiag_Trms=0         ! DIAG_DOC: $\sqrt{\left< T^2 \right>}$
  integer :: idiag_uxTm=0         ! DIAG_DOC: $\left< u_x T \right>$
  integer :: idiag_uyTm=0         ! DIAG_DOC: $\left< u_y T \right>$
  integer :: idiag_uzTm=0         ! DIAG_DOC: $\left< u_z T \right>$
  integer :: idiag_gT2m=0         ! DIAG_DOC: $\left< (\nabla T)^2 \right>$
  integer :: idiag_guxgTm=0       ! DIAG_DOC: $\left< \nabla u_x \cdot \nabla T \right>$
  integer :: idiag_guygTm=0       ! DIAG_DOC: $\left< \nabla u_y \cdot \nabla T \right>$
  integer :: idiag_guzgTm=0       ! DIAG_DOC: $\left< \nabla u_z \cdot \nabla T \right>$
  integer :: idiag_Tugux_uxugTm=0 ! DIAG_DOC: $\left< T \uv\cdot\nabla u_x + u_x \uv\cdot\nabla T \right>
                                  ! DIAG_DOC:  =\left< \uv\cdot\nabla(u_x T) \right>$
  integer :: idiag_Tuguy_uyugTm=0 ! DIAG_DOC: $\left< T \uv\cdot\nabla u_y + u_y \uv\cdot\nabla T \right>
                                  ! DIAG_DOC:  =\left< \uv\cdot\nabla(u_y T) \right>$
  integer :: idiag_Tuguz_uzugTm=0 ! DIAG_DOC: $\left< T \uv\cdot\nabla u_z + u_z \uv\cdot\nabla T \right>
                                  ! DIAG_DOC:  =\left< \uv\cdot\nabla(u_z T) \right>$
  integer :: idiag_Tdxpm=0        ! DIAG_DOC: $\left< T dp/dx \right>$
  integer :: idiag_Tdypm=0        ! DIAG_DOC: $\left< T dp/dy \right>$
  integer :: idiag_Tdzpm=0        ! DIAG_DOC: $\left< T dp/dz \right>$
!
  integer :: idiag_fradtop=0  ! DIAG_DOC: $<-K{dT\over dz}>_{\text{top}}$
                              ! DIAG_DOC: \quad(top radiative flux)
  integer :: idiag_fradbot=0  ! DIAG_DOC: $<-K{dT\over dz}>_{\text{bot}}$
                              ! DIAG_DOC: \quad(bottom radiative flux)
  integer :: idiag_yHmax=0    ! DIAG_DOC: DOCUMENT ME
  integer :: idiag_yHmin=0    ! DIAG_DOC: DOCUMENT ME
  integer :: idiag_yHm=0      ! DIAG_DOC: DOCUMENT ME
  integer :: idiag_ethm=0     ! DIAG_DOC: $\left< e_{\text{th}}\right> =
                              ! DIAG_DOC:  \left< c_v \rho T \right> $
                              ! DIAG_DOC: \quad(mean thermal energy)
  integer :: idiag_eem=0      ! DIAG_DOC: $\left< e \right> =
                              ! DIAG_DOC:  \left< c_v T \right>$
                              ! DIAG_DOC: \quad(mean internal energy)
  integer :: idiag_ethtot=0   ! DIAG_DOC: $\int_V\varrho e\,dV$
                              ! DIAG_DOC:   \quad(total thermal energy)
  integer :: idiag_ssm=0      ! DIAG_DOC: $\overline{S}$
  integer :: idiag_thcool=0   ! DIAG_DOC: $\tau_{\rm cool}$
  integer :: idiag_ppm=0      ! DIAG_DOC: $\overline{P}$
  integer :: idiag_csm=0      ! DIAG_DOC: $\overline{c}_{\rm s}$
  integer :: idiag_csmax=0    ! DIAG_DOC: $\max (c_{\rm s})$
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
!
! xy averaged diagnostics given in xyaver.in written every it1d timestep
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
  integer :: idiag_TT2mz=0      ! XYAVG_DOC: $\left<T^2\right>_{xy}$
  integer :: idiag_uxTmz=0      ! XYAVG_DOC: $\left<u_x T\right>_{xy}$
  integer :: idiag_uyTmz=0      ! XYAVG_DOC: $\left<u_y T\right>_{xy}$
  integer :: idiag_uzTmz=0      ! XYAVG_DOC: $\left<u_z T\right>_{xy}$
  integer :: idiag_fradmz=0     ! XYAVG_DOC: $\left<F_{\rm rad}\right>_{xy}$
  integer :: idiag_fconvmz=0    ! XYAVG_DOC: $\left<c_p \varrho u_z T \right>_{xy}$
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
  real, dimension(:,:), allocatable :: pp_xz,pp_yz,pp_xy,pp_xy2,pp_xy3,pp_xy4,pp_xz2
!
! y averaged diagnostics given in yaver.in
!
  integer :: idiag_TTmxz=0      ! YAVG_DOC: $\left<T\right>_{y}$
  integer :: idiag_EmAIA94mxz=0 ! YAVG_DOC: Emission off AIA 94 channel
                                ! YAVG_DOC: integrated over y direction
  integer :: idiag_EmAIA131mxz=0! YAVG_DOC: Emission off AIA 131 channel
                                ! YAVG_DOC: integrated over y direction
  integer :: idiag_EmAIA171mxz=0! YAVG_DOC: Emission off AIA 171 channel
                                ! YAVG_DOC: integrated over y direction
  integer :: idiag_EmAIA193mxz=0! YAVG_DOC: Emission off AIA 193 channel
                                ! YAVG_DOC: integrated over y direction
  integer :: idiag_EmAIA211mxz=0! YAVG_DOC: Emission off AIA 211 channel
                                ! YAVG_DOC: integrated over y direction
  integer :: idiag_EmAIA304mxz=0! YAVG_DOC: Emission off AIA 304 channel
                                ! YAVG_DOC: integrated over y direction
  integer :: idiag_EmAIA335mxz=0! YAVG_DOC: Emission off AIA 335 channel
                                ! YAVG_DOC: integrated over y direction
  integer :: idiag_EmXRTmxz=0   ! YAVG_DOC: Emission off XRT Al-poly channel
                                ! YAVG_DOC: integrated over y direction
!
! z averaged diagnostics given in zaver.in
!
  integer :: idiag_TTmxy=0      ! ZAVG_DOC: $\left<T\right>_{z}$
  integer :: idiag_EmAIA94mxy=0 ! ZAVG_DOC: Emission off AIA 94 channel
                                ! ZAVG_DOC: integrated over z direction
  integer :: idiag_EmAIA131mxy=0! ZAVG_DOC: Emission off AIA 131 channel
                                ! ZAVG_DOC: integrated over z direction
  integer :: idiag_EmAIA171mxy=0! ZAVG_DOC: Emission off AIA 171 channel
                                ! ZAVG_DOC: integrated over z direction
  integer :: idiag_EmAIA193mxy=0! ZAVG_DOC: Emission off AIA 193 channel
                                ! ZAVG_DOC: integrated over z direction
  integer :: idiag_EmAIA211mxy=0! ZAVG_DOC: Emission off AIA 211 channel
                                ! ZAVG_DOC: integrated over z direction
  integer :: idiag_EmAIA304mxy=0! ZAVG_DOC: Emission off AIA 304 channel
                                ! ZAVG_DOC: integrated over z direction
  integer :: idiag_EmAIA335mxy=0! ZAVG_DOC: Emission off AIA 335 channel
                                ! ZAVG_DOC: integrated over z direction
  integer :: idiag_EmXRTmxy=0   ! ZAVG_DOC: Emission off XRT Al-poly channel
                                ! ZAVG_DOC: integrated over z direction
!
  integer :: ivid_pp=0
!
  real, dimension(nx) :: diffus_chi,diffus_chi3,hcond
!
  contains
!***********************************************************************
    subroutine register_energy
!
!  Initialise variables which should know that we solve an energy
!  equation: ilnTT, etc; increase nvar accordingly.
!
!  6-nov-01/wolf: coded
! 18-may-12/MR: shared variable PrRa fetched from hydro
!
      use BorderProfiles, only: request_border_driving
      use FArrayManager, only: farray_register_pde, farray_index_append
      use SharedVariables, only: get_shared_variable
!
!  Register TT or lnTT, depending on whether or not ltemperature_nolog
!
      if (ltemperature_nolog) then
        call farray_register_pde('TT',iTT)
        ilnTT=iTT
        !call farray_index_append('ilnTT',ilnTT)
      else
        call farray_register_pde('lnTT',ilnTT)
      endif
!
!  logical variable lpressuregradient_gas shared with hydro modules
!
      call get_shared_variable('lpressuregradient_gas',lpressuregradient_gas, &
                               caller='register_energy')
!
!  real variable PrRa shared with hydro modules, used for Boussinesq
!
      if (lboussinesq.and.lviscosity_heat) &
        call get_shared_variable('PrRa',PrRa)
!
!  Tell the BorderProfiles module if we intend to use border driving, so
!  that the module can request the right pencils.
!
      if (borderss/='nothing') call request_border_driving(borderss)
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
    endsubroutine register_energy
!***********************************************************************
    subroutine initialize_energy(f)
!
!  Called by run.f90 after reading parameters, but before the time loop.
!
!  21-jul-2002/wolf: coded
!
      use FArrayManager, only: farray_register_global
      use Gravity, only: gravz, compute_gravity_star
      use EquationOfState, only : cs2bot, cs2top, gamma, gamma_m1, &
                                  select_eos_variable
      use Sub, only: step,der_step
      use SharedVariables, only: put_shared_variable, get_shared_variable
      use Mpicomm, only: stop_it
      !use Slices_methods, only: alloc_slice_buffers
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: dhcond
      logical :: lnothing
      integer :: i
      logical, pointer :: lrss
      real :: star_cte

      if (lstart .and. .not.leos) &
        call warning('initialize_energy', &
            'EOS=noeos, but possibly an EQUATION OF STATE is needed')
!
!  Set iTT equal to ilnTT if we are considering non-logarithmic temperature.
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
      lheatc_chicubicstep = .false.
      lheatc_kramers= .false.
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
          if (lroot) call information('initialize_energy', &
          ' heat conduction: K=cst --> gamma*K/rho/TT/cp*div(T*grad lnTT)')
        case ('K-profile')
          lheatc_Kprof=.true.
          lmultilayer=.true.
          if (lroot) call information('initialize_energy', &
                ' heat conduction: K=K(z) or K=K(r)')
        case ('K-arctan')
          lheatc_Karctan=.true.
          if (.not. ltemperature_nolog) &
            call fatal_error('initialize_energy', &
              'K-arctan only valid for TT')
          if (lADI_mixed .and. .not. lADI) &
            call fatal_error('initialize_energy', &
              'K-arctan with lADI_mixed=T while lADI=F?')
          if (lroot) call information('initialize_energy', &
              'heat conduction: arctan profile')
        case ('chi-const')
          lheatc_chiconst=.true.
          if (lroot) call information('initialize_energy', &
              ' heat conduction: constant chi')
        case ('chi-cubicstep')
          lheatc_chicubicstep=.true.
          if (lroot) call information('initialize_energy', &
              ' heat conduction: cubic step profile of chi')
        case('K-therm')
          lheatc_Ktherm=.true.
          if (lroot) call information('initialize_energy', &
              ' heat conduction: temperature dependent K')
        case ('chi-hyper3')
          lheatc_hyper3=.true.
          if (lroot) call information('initialize_energy','hyper conductivity')
        case ('hyper3_mesh','hyper3-mesh')
          lheatc_hyper3_mesh=.true.
          if (lroot) call information('initialize_energy','hyper mesh conductivity')
        case ('hyper3_cyl','hyper3-cyl','hyper3_sph','hyper3-sph')
          lheatc_hyper3_polar=.true.
          if (lroot) call information('initialize_energy', &
              'hyper conductivity: polar coords')
        case ('shock','chi-shock')
          lheatc_shock=.true.
          if (lroot) call information('initialize_energy','shock conductivity')
        case ('tensor-diffusion')
          lheatc_tensordiffusion=.true.
          if (lroot) print*, 'heat conduction: tensor diffusion'
        case ('kramers')
          lheatc_kramers=.true.
          if (lroot) print*, 'heat conduction: kramers'
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
      if (bcz12(ilnTT,1)=='c1' .and. lrun) then

        if (.not. leos) &
          call fatal_error('initialize_energy', &
              'EOS=noeos, but BC "c1" requires an EQUATION OF STATE')
!
        if (Fbot==impossible .and. hcond0 /= impossible) then
          Fbot=-gamma/gamma_m1*hcond0*gravz/(mpoly0+1.0)
          if (lroot) print*, &
              'initialize_energy: Calculated Fbot = ', Fbot
        endif
        if (hcond0==impossible .and. Fbot /= impossible) then
          hcond0=-Fbot*gamma_m1/gamma*(mpoly0+1.0)/gravz
          if (lroot) print*, &
              'initialize_energy: Calculated hcond0 = ', hcond0
        endif
        if (Fbot==impossible .and. hcond0==impossible) &
          call fatal_error('temperature_idealgas',  &
              'Both Fbot and hcond0 are unknown')
      endif
!
      if (initlnTT(1)=='star_heat') then

        if (.not. leos) &
          call fatal_error('initialize_energy', &
              'EOS=noeos, but star heating requires an EQUATION OF STATE')
!
        if (lroot) print*,'star_heat: compute the gravity profile'
        ! compute the gravity profile
        star_cte=(mpoly0+1.)/hcond0*gamma_m1/gamma
        call compute_gravity_star(f, wheat, luminosity, star_cte)
        if (rcool==0.) rcool=r_ext
      endif
!
      if (lmultilayer) then

        if (.not. leos) &
          call fatal_error('initialize_energy', &
              'EOS=noeos, but multilayer requires an EQUATION OF STATE')
!
        hcond1=(mpoly1+1.)/(mpoly0+1.)
        hcond2=(mpoly2+1.)/(mpoly0+1.)
      endif
!
!  Now we share several variables.
!
      call put_shared_variable('hcond0', hcond0, caller='initialize_energy')
      call put_shared_variable('hcond1', hcond1)
      call put_shared_variable('hcond2', hcond2)
      call put_shared_variable('lmultilayer', lmultilayer)
      call put_shared_variable('widthlnTT', widthlnTT)
      call put_shared_variable('Fbot', Fbot)
      call put_shared_variable('lADI_mixed', lADI_mixed)
      call put_shared_variable('lviscosity_heat',lviscosity_heat)
      call put_shared_variable('mpoly0', mpoly0)
      call put_shared_variable('mpoly1', mpoly1)
      call put_shared_variable('mpoly2', mpoly2)
      call put_shared_variable('hcond0_kramers', hcond0_kramers)
      call put_shared_variable('nkramers', nkramers)
      call put_shared_variable('lheatc_kramers', lheatc_kramers)
      if (lsolid_cells) then
        if (.not. lchemistry) call put_shared_variable('chi', chi)
        call put_shared_variable('ladvection_temperature',ladvection_temperature)
        call put_shared_variable('lheatc_chiconst',lheatc_chiconst)
        call put_shared_variable('lupw_lnTT',lupw_lnTT)
      endif
!
!  Share the 4 parameters of the radiative conductivity hole (kappa-mechanism
!  problem).
!
      hole_params=(/Tbump,Kmin,Kmax,hole_slope,hole_width/)
      call put_shared_variable('hole_params',hole_params)
!
!  A word of warning...
!
      if (lheatc_Kconst .and. hcond0==0.0) &
        call warning('initialize_energy', 'hcond0 is zero!')
      if (lheatc_Ktherm .and. hcond0==0.0) &
        call warning('initialize_energy','hcond0 is zero!')
      if (lheatc_Kprof .and. hcond0==0.0) &
        call warning('initialize_energy', 'hcond0 is zero!')
      if (lheatc_chiconst .and. chi==0.0) &
        call warning('initialize_energy','chi is zero!')
      if (lheatc_chicubicstep .and. chi==0.0) &
        call warning('initialize_energy','chi is zero!')
      if (lheatc_kramers .and. hcond0==0.0) &
        call warning('initialize_energy','hcond0 is zero!')
      if (lrun) then
        if (lheatc_hyper3 .and. chi_hyper3==0.0) &
            call fatal_error('initialize_energy', &
            'Conductivity coefficient chi_hyper3 is zero!')
        if (lheatc_shock .and. chi_shock==0.0) &
            call fatal_error('initialize_energy', &
            'Conductivity coefficient chi_shock is zero!')
      endif
      if (iheatcond(1)=='nothing') then
        if (hcond0/=impossible) call warning('initialize_energy', &
            'No heat conduction, but hcond0/=0')
        if (chi/=impossible) call warning('initialize_energy', &
            'No heat conduction, but chi/=0')
      endif
      if (lADI_mixed .and. iheatcond(1) /= 'K-arctan') &
        call stop_it("temperature_idealgas: "//&
          "lADI_mixed=T while iheatcond /= K-arctan?")
!
      if (llocal_iso) &
           call fatal_error('initialize_energy', &
           'llocal_iso switches on the local isothermal approximation. ' // &
           'Use ENERGY=noenergy in src/Makefile.local')
!
!  Compute mask for z-averaging where z is in temp_zaver_range.
!  Normalize such that the average over the full domain
!  gives still unity.
!
      if (n1 == n2) then
        zmask_temp = 1.
        zmask_temp_global=1.
      else
        where (z(n1:n2) >= temp_zaver_range(1) .and. z(n1:n2)  <= temp_zaver_range(2))
          zmask_temp = 1.
        elsewhere
          zmask_temp = 0.
        endwhere
        where (zglobal(nghost+1:mzgrid-nghost) >= temp_zaver_range(1) .and. &
               zglobal(nghost+1:mzgrid-nghost) <= temp_zaver_range(2))
          zmask_temp_global = 1.
        elsewhere
          zmask_temp_global = 0.
        endwhere
        temp_zmask_count =  max(count(zmask_temp_global ==  1.0),1)
        zmask_temp = zmask_temp*float(nzgrid)/float(temp_zmask_count)
      endif
!
!  debug output
!
      if (lroot.and.ip<14) then
        print*,'zmask_temp=' ,zmask_temp
      endif
!
!
!  Check if reduced sound speed is used
!
      if (ldensity) then
        call get_shared_variable('lreduced_sound_speed', lreduced_sound_speed)  !lrss)
        !!!lreduced_sound_speed=lrss
        if (lreduced_sound_speed) then
          call get_shared_variable('reduce_cs2',reduce_cs2)
!          call get_shared_variable('lscale_to_cs2top',lscale_to_cs2top)
        endif
      endif
!
      if (ivid_pp/=0) then
        !call alloc_slice_buffers(pp_xy,pp_xz,pp_yz,pp_xy2,pp_xy3,pp_xy4,pp_xz2)
        if (lwrite_slice_xy .and..not.allocated(pp_xy) ) allocate(pp_xy (nx,ny))
        if (lwrite_slice_xz .and..not.allocated(pp_xz) ) allocate(pp_xz (nx,nz))
        if (lwrite_slice_yz .and..not.allocated(pp_yz) ) allocate(pp_yz (ny,nz))
        if (lwrite_slice_xy2.and..not.allocated(pp_xy2)) allocate(pp_xy2(nx,ny))
        if (lwrite_slice_xy3.and..not.allocated(pp_xy3)) allocate(pp_xy3(nx,ny))
        if (lwrite_slice_xy4.and..not.allocated(pp_xy4)) allocate(pp_xy4(nx,ny))
        if (lwrite_slice_xz2.and..not.allocated(pp_xz2)) allocate(pp_xz2(nx,nz))
      endif
!
    endsubroutine initialize_energy
!***********************************************************************
    subroutine init_energy(f)
!
!  Initialise lnTT or TT; called from start.f90.
!
!  13-dec-2002/axel+tobi: adapted from init_energy
!
!  initialise energy; called from start.f90
!  07-nov-01/wolf: coded
!  24-nov-02/tony: renamed for consistency (i.e. init_[variable name])
!  12-may-12/MR: initialization with mode added
!
      use General,  only: itoa
      use Sub,      only: blob
      use InitialCondition, only: initial_condition_ss
      use EquationOfState, only: gamma, gamma_m1, cs2bot, cs2top, cs20, &
                                 lnrho0, get_cp1, rho0
      use Gravity, only: gravz
      use Mpicomm, only: stop_it
      use Initcond, only: modes
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      real, dimension (mz) :: TTz
!
      integer :: j
      logical :: lnothing=.true.
      real :: haut, Rgas, cp1, Ttop, alpha, beta, expo, ztop
!
      do j=1,ninit
!
        if (initlnTT(j)/='nothing') then
!
          lnothing=.false.
!
          iinit_str=itoa(j)
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
            cs2bot=gamma_m1*TT_const
            cs2top=gamma_m1*TT_const
!
          case ('const_dTTdz')
            TTz=TT_const+z
            if (ltemperature_nolog) then
              f(:,:,:,iTT)=f(:,:,:,iTT)+spread(spread(TTz,1,mx),2,my)
            else
              f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+spread(spread(log(TTz),1,mx),2,my)
            endif
            cs2bot=gamma_m1*TTz(n1)
            cs2top=gamma_m1*TTz(n2)
!
          case ('linz_fromBC')
!
!  Linear profile between temperature values at bottom and top supposed to be given
!  in fbcz.
!
            if (ltemperature_nolog) then
              TTz=fbcz(iTT,1)+(z-xyz0(3))/Lxyz(3)*(fbcz(iTT,2)-fbcz(iTT,1))
              f(:,:,:,iTT)=f(:,:,:,iTT)+spread(spread(TTz,1,mx),2,my)
              cs2bot=gamma_m1*fbcz(iTT,1)
              cs2top=gamma_m1*fbcz(iTT,2)
            else
              call fatal_error('init_ss','BC "linz_fromBC" not implemented for logarithmic temperature')
            endif
!
          case ('mode')
!
            if (ltemperature_nolog) then
              call modes(ampl_lnTT(j),coef_lnTT,f,iTT,kx_lnTT,ky_lnTT,kz_lnTT)
            else
              call modes(ampl_lnTT(j),coef_lnTT,f,ilnTT,kx_lnTT,ky_lnTT,kz_lnTT)
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
                radius_lnTT, ampl_lnTT(j), center1_x, center1_y, center1_z
            call blob(ampl_lnTT(j),f,ilnTT,radius_lnTT, &
                center1_x, center1_y,center1_z)
            call blob(-ampl_lnTT(j),f,ilnrho,radius_lnTT, &
                center1_x,center1_y,center1_z)
!
          case ('blob')
            if (lroot) print*, 'init_lnTT: blob ', &
                radius_lnTT, ampl_lnTT(j), center1_x, center1_y, center1_z
            call blob(ampl_lnTT(j),f,ilnTT,radius_lnTT, &
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
            ztop=xyz1(3)
            if (ldensity_nolog) then
              do n=n1,n2
                f(:,:,n,irho)=rho0*exp((ztop-z(n))/haut)
              enddo
            else
              do n=n1,n2
                f(:,:,n,ilnrho)=lnrho0+(ztop-z(n))/haut
              enddo
            endif
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
            write(unit=errormsg,fmt=*) 'No such value for initlnTT(' &
                           //trim(iinit_str)//'): ',trim(initlnTT(j))
            call fatal_error('init_energy',errormsg)
!
          endselect
!
          if (lroot) print*,'init_energy: initlnTT(' &
              //trim(iinit_str)//') = ',trim(initlnTT(j))
        endif
      enddo
!
!  Interface for user's own initial condition.
!
      if (linitial_condition) call initial_condition_ss(f)
!
      if (lnothing.and.lroot) print*,'init_energy: nothing'
!
      if (ltemperature_nolog.and.linitial_log) f(:,:,:,iTT)=exp(f(:,:,:,ilnTT))
!
    endsubroutine init_energy
!***********************************************************************
    subroutine pencil_criteria_energy
!
!  All pencils that the Energy module depends on are specified here.
!
!  20-11-04/anders: coded
!
      if (lwrite_slices.and.ivid_pp/=0) lpenc_video(i_pp)=.true.
!
!  cs2 affects time step only if the continuity equation is being solved,
!  i.e. not for boussinesq or anelastic.
!
      if (ldt.and.ldensity) lpenc_requested(i_cs2)=.true.
!
      if (lpressuregradient_gas) lpenc_requested(i_fpres)=.true.
      if (lparticles_temperature) lpenc_requested(i_tcond)=.true.
!
      if (lviscosity.and.lviscosity_heat) then
        lpenc_requested(i_cv1)=.true.
        lpenc_requested(i_visc_heat)=.true.
        if (.not.ltemperature_nolog) &
            lpenc_requested(i_TT1)=.true.
      endif
!
      if (ldensity.or.lanelastic) then
        lpenc_requested(i_divu)=.true.
        if (ltemperature_nolog) lpenc_requested(i_TT)=.true.
      endif
!
      if (lcalc_heat_cool) then
        if (ldensity.or.lboussinesq.or.lanelastic) then
          lpenc_requested(i_rho1)=.true.
          lpenc_requested(i_TT)=.true.
          lpenc_requested(i_TT1)=.true.
          lpenc_requested(i_cv1)=.true.
        endif
        if (lgravr) then
          lpenc_requested(i_r_mn)=.true.
          if (lboussinesq) lpenc_requested(i_evr)=.true.
        endif
      endif
!
      if (lheatc_chiconst.or.lheatc_chicubicstep) then
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
      if (lheatc_chicubicstep) lpenc_requested(i_z_mn)=.true.
!
      if (lheatc_Kconst) then
        if (ldensity.or.lboussinesq.or.lanelastic) lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_cp1)=.true.
        if (ltemperature_nolog) then
          lpenc_requested(i_del2TT)=.true.
        else
          lpenc_requested(i_glnTT)=.true.
          lpenc_requested(i_del2lnTT)=.true.
        endif
      endif
!
      if (lheatc_Ktherm) then
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
      if (lheatc_kramers) then

        !lpenc_requested(i_rho)=.true.
        !lpenc_requested(i_lnrho)=.true.
        !lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_glnTT)=.true.
        !lpenc_requested(i_cv1)=.true.
        lpenc_requested(i_del2lnTT)=.true.
      endif
!
      if (ladvection_temperature) then
        if (ltemperature_nolog) then
          lpenc_requested(i_ugTT)=.true.
          lpenc_requested(i_del2TT)=.true.
        else
          lpenc_requested(i_uglnTT)=.true.
        endif
      endif
!
!  Diagnostic pencils.
!
      if ( idiag_TTmax/=0.or.idiag_TTmin/=0 .or.idiag_TTm/=0  .or.idiag_TugTm/=0 .or. &
           idiag_Trms/=0 .or.idiag_uxTm/=0  .or.idiag_uyTm/=0 .or.idiag_uzTm/=0  .or. &
           idiag_TT2mz/=0 .or.idiag_uxTmz/=0.or.idiag_uyTmz/=0.or.idiag_uzTmz/=0 .or. &
           idiag_TT2m/=0 .or.idiag_TTzmask/=0 ) lpenc_diagnos(i_TT)=.true.
!
      if ( idiag_TugTm/=0 .or. idiag_gT2m/=0 .or. &
           idiag_guxgTm/=0 .or. idiag_guygTm/=0 .or. idiag_guzgTm/=0 ) lpenc_diagnos(i_gTT)=.true.
!
      if ( lhydro .or. lhydro_kinematic ) then
        if ( idiag_guxgTm/=0 .or. idiag_guygTm/=0 .or. idiag_guzgTm/=0 ) lpenc_diagnos(i_uij)=.true.
        if ( idiag_Tugux_uxugTm/=0 .or. idiag_Tuguy_uyugTm/=0 .or. idiag_Tuguz_uzugTm/=0 ) then
          lpenc_requested(i_ugu)=.true.; lpenc_requested(i_ugTT)=.true.
        endif
      endif
!
      if ( idiag_Tdxpm/=0 .or. idiag_Tdypm/=0 .or. idiag_Tdzpm/=0 ) then
        lpenc_diagnos(i_fpres)=.true.
        lpenc_diagnos(i_TT)=.true.
      endif
      if (idiag_gTmax/=0) then
         lpenc_diagnos(i_glnTT) =.true.
         lpenc_diagnos(i_TT) =.true.
      endif
      if (idiag_fradtop/=0.or.idiag_fradbot/=0.or.idiag_fradmz/=0) then
        lpenc_diagnos(i_TT) =.true.
        lpenc_diagnos(i_glnTT) =.true.
      endif
      if (idiag_fconvmz/=0) then
        lpenc_diagnos(i_TT) =.true.
        lpenc_diagnos(i_cp)=.true.
        lpenc_diagnos(i_uu)=.true.
        lpenc_diagnos(i_rho)=.true.
      endif
      if (idiag_yHmax/=0) lpenc_diagnos(i_yH)  =.true.
      if (idiag_yHmin/=0) lpenc_diagnos(i_yH)  =.true.
      if (idiag_yHm/=0)   lpenc_diagnos(i_yH)  =.true.
      if (idiag_ethm/=0.or.idiag_ethmz/=0.or.idiag_ethtot/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_ee)  =.true.
      endif
      if (idiag_ethuxmz/=0.or.idiag_ethuymz/=0.or.idiag_ethuzmz/=0.or.&
          idiag_ethuxmx/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_ee) =.true.
        lpenc_diagnos(i_uu) =.true.
      endif
      if (idiag_uxTmz/=0.or.idiag_uyTmz/=0.or.idiag_uzTmz/=0) &
        lpenc_diagnos(i_uu) =.true.
      if (idiag_ssm/=0)    lpenc_diagnos(i_ss)  =.true.
      if (idiag_dtchi/=0)  lpenc_diagnos(i_cs2)=.true.
      if (idiag_csm/=0 .or. idiag_csmax/=0)   lpenc_diagnos(i_cs2)=.true.
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
      if (idiag_EmAIA94mxy/=0 .or. idiag_EmAIA94mxz/=0 .or.  &
          idiag_EmAIA131mxy/=0 .or. idiag_EmAIA131mxz/=0 .or. &
          idiag_EmAIA171mxy/=0 .or. idiag_EmAIA171mxz/=0 .or. &
          idiag_EmAIA193mxy/=0 .or. idiag_EmAIA193mxz/=0 .or. &
          idiag_EmAIA211mxy/=0 .or. idiag_EmAIA211mxz/=0 .or. &
          idiag_EmAIA304mxy/=0 .or. idiag_EmAIA304mxz/=0 .or. &
          idiag_EmAIA335mxy/=0 .or. idiag_EmAIA335mxz/=0 .or. &
          idiag_EmXRTmxy/=0 .or. idiag_EmXRTmxz/=0 ) then
        lpenc_diagnos2d(i_TT)=.true.
        lpenc_diagnos2d(i_rho)=.true.
      endif
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
      if (lpencil_in(i_fpres).and. .not.lboussinesq) then
        lpencil_in(i_cs2)=.true.
        lpencil_in(i_glnrho)=.true.
        lpencil_in(i_glnTT)=.true.
        lpencil_in(i_glnmumol)=.true.
      endif

      if (lpencil_in(i_glnTT).and.ltemperature_nolog.and..not.leos) &
        lpencil_in(i_TT1) =.true.
!
    endsubroutine pencil_interdep_energy
!***********************************************************************
    subroutine calc_pencils_energy(f,p)
!
!  Calculate Energy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!  31-01-18/MR: made calculation of p%gTT corrrect also for log temperature
!
      use EquationOfState, only: gamma1
      use Sub, only: u_dot_grad,grad,multmv,del2
      use Deriv, only: der2
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      type (pencil_case), intent (inout) :: p
      integer :: j
      real, dimension(nx,3) :: gpp
      real, dimension(nx) :: temp
!
      if (.not.leos) then
        if (lpencil(i_lnTT)) then
          if (ltemperature_nolog) then
            p%lnTT=alog(f(l1:l2,m,n,iTT))
          else
            p%lnTT=f(l1:l2,m,n,ilnTT)
          endif
        endif
        if (lpencil(i_TT)) then
          if (ltemperature_nolog) then
            p%TT=f(l1:l2,m,n,iTT)
          else
            p%TT=exp(f(l1:l2,m,n,ilnTT))
          endif
        endif
        if (lpencil(i_TT1)) then
          if (ltemperature_nolog) then
            p%TT1=1./f(l1:l2,m,n,iTT)
          else
            p%TT1=exp(-f(l1:l2,m,n,ilnTT))
          endif
        endif
        if (lpencil(i_glnTT)) then
          call grad(f,ilnTT,p%glnTT)
          if (ltemperature_nolog) then
            do j=1,3
              p%glnTT(:,j) = p%glnTT(:,j)*p%TT1
            enddo
          endif
        endif

        if (ltemperature_nolog) then
          if (lpencil(i_del2TT)) call del2(f,iTT,p%del2TT)
        else
          call fatal_error('calc_pencils_energy','del2TT, not implmented for logaritthmic temp and noeos')
        endif

        if (lpencil(i_hlnTT).or.lpencil(i_del2lnTT).or.lpencil(i_del6TT).or.lpencil(i_del6lnTT)) &
          call fatal_error('calc_pencils_energy','hlnTT, del[26]TT, del[26]lnTT, not implmented for noeos')

      endif
! Ma2
      if (lpencil(i_Ma2)) p%Ma2=p%u2/p%cs2
! uglnTT
      if (lpencil(i_uglnTT)) &
          call u_dot_grad(f,ilnTT,p%glnTT,p%uu,p%uglnTT,UPWIND=lupw_lnTT)
! ugTT
      if (lpencil(i_ugTT)) &
          call u_dot_grad(f,ilnTT,p%gTT,p%uu,p%ugTT,UPWIND=lupw_lnTT)
! glnTT
      if (lpencil(i_gTT)) then
        call grad(f,ilnTT,p%gTT)
        if (.not.ltemperature_nolog) then
          temp=exp(f(l1:l2,m,n,ilnTT))
          do j=1,3
            p%gTT(:,j) = p%gTT(:,j)*temp
          enddo
        endif
        do j=1,3
          if (gradTT0(j)/=0.) p%gTT(:,j)=p%gTT(:,j)+gradTT0(j)
        enddo
      endif

      if (lpencil(i_d2xlnTT)) call der2(f,ilnTT,p%d2xlnTT,1)
      if (lpencil(i_d2zlnTT)) call der2(f,ilnTT,p%d2zlnTT,3)
!
! fpres
      if (lpencil(i_fpres)) then
        if (lboussinesq) then
!
! subroutine boussinesq has already been called at this instant,
! so f(:,:,:,ipp) contains div(uu), which is in fact the change of div(uu) due to the last
! timestep, hence p.dt
!
            call grad(f,ipp,gpp)
            if (dt==0.) then
              p%fpres=gpp
            else
              p%fpres=gpp/dt
            endif
        else
          do j=1,3
            p%fpres(:,j)=-gamma1*p%cs2* &
                (p%glnrho(:,j)+p%glnTT(:,j)-p%glnmumol(:,j))
          enddo
        endif
      endif
! tcond
      if (lpencil(i_tcond)) then
        if (lheatc_chiconst) then
          p%tcond=chi*p%rho/p%cp1
        elseif (lheatc_Kconst) then
          p%tcond=hcond0
        else
          call fatal_error('calc_pencils_energy',  &
              'This heatcond is not implemented to work with lpencil(i_tcond)!')
        endif
      endif
! sglnTT 
      if (lpencil(i_sglnTT)) then
        call multmv(p%sij,p%glnTT,p%sglnTT)
      endif
!
    endsubroutine calc_pencils_energy
!***********************************************************************
    subroutine denergy_dt(f,df,p)
!
!  Calculate right hand side of temperature equation.
!
!  lnTT version: DlnTT/Dt = -gamma_m1*divu + gamma*cp1*rho1*TT1*RHS
!    TT version:   DTT/Dt = -gamma_m1*TT*divu + gamma*cp1*rho1*RHS
!
!  13-dec-02/axel+tobi: adapted from energy
!  18-may-12/MR: compression work as heat sink added for boussinesq
!
      use Deriv, only: der6
      use EquationOfState, only: gamma_m1, lpres_grad
      use ImplicitPhysics, only: heatcond_TT
      use Special, only: special_calc_energy
      use Sub, only: identify_bcs
      use Viscosity, only: calc_viscous_heat
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: Hmax=0.0, thdiff, tmp, advec_hypermesh_ss
      integer :: j
!
      intent(inout) :: f,p,df
!
! When initializating  thdiff in the declaration the
! variable gets the SAVE attribute,
! so in the next call thdiff is not initialized anymore.
!
      thdiff = 0.0
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*, 'SOLVE dlnTT_dt'
      if (headtt) then
        if (ltemperature_nolog) then
          print*, 'denergy_dt: TT,cs2=', p%TT(1), p%cs2(1)
          call identify_bcs('TT',iTT)
        else
          print*, 'denergy_dt: lnTT,cs2=', p%lnTT(1), p%cs2(1)
          call identify_bcs('lnTT',ilnTT)
        endif
      endif
!
!  ``cs2/dx^2'' for timestep
!
      if (lfirst.and.ldt.and.ldensity.and.lhydro) then
        if (lreduced_sound_speed) then
!          if (lscale_to_cs2top) then
!            advec_cs2=reduce_cs2*cs2top*dxyz_2
!          else
            advec_cs2=reduce_cs2*p%cs2*dxyz_2
!          endif
        else
          advec_cs2=p%cs2*dxyz_2
        endif
        if (headtt.or.ldebug) print*,'denergy_dt: max(advec_cs2) =',maxval(advec_cs2)
      endif
!
!  Sound speed squared.
!
      if (headtt) print*, 'denergy_dt: cs20=', p%cs2(1)
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
!  Add divu term.
!
      if (ldensity.or.lanelastic) then
        if (ltemperature_nolog) then
          df(l1:l2,m,n,iTT)   = df(l1:l2,m,n,iTT)   - gamma_m1*p%TT*p%divu
        else
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - gamma_m1*p%divu
        endif
!
! if viscous heating is allowed, compression work must be taken into account, too
!
      elseif (lboussinesq.and.lviscosity.and.lviscosity_heat) then
        if (ltemperature_nolog) then
          df(l1:l2,m,n,iTT)   = df(l1:l2,m,n,iTT)   - p%cv1*PrRa*f(l1:l2,m,n,iTT)*p%uu(:,3)
        else
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%cv1*PrRa*p%uu(:,3)
        endif
      endif
!
!  Calculate viscous contribution to temperature.
!
      if (lviscosity.and.lviscosity_heat) call calc_viscous_heat(df,p,Hmax)
!
!  Various heating conduction contributions.
!
      if (lcalc_heat_cool)  call calc_heat_cool(f,df,p)
!
!  Thermal conduction
!
      diffus_chi=0.; diffus_chi3=0.
      if (lheatc_chiconst) call calc_heatcond_constchi(df,p)
      if (lheatc_chicubicstep) call calc_heatcond_cubicstepchi(df,p)
      if (lheatc_Kconst)   call calc_heatcond_constK(df,p)
      if (lheatc_Kprof)    call calc_heatcond(f,df,p)
      if (lheatc_Karctan)  call calc_heatcond_arctan(df,p)
      if (lheatc_tensordiffusion) call calc_heatcond_tensor(df,p)
      if (lheatc_Ktherm) call calc_heatcond_Ktherm (df,p)
      if (lheatc_kramers) call calc_heatcond_kramers (df,p)
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
        if (headtt) print*,'denergy_dt: chi_hyper3=', chi_hyper3
      endif
!
      if (lheatc_hyper3_mesh) then
        do j=1,3
          call der6(f,ilnTT,tmp,j,IGNOREDX=.true.)
          if (.not.ltemperature_nolog) tmp=tmp*p%TT1
          thdiff = thdiff + chi_hyper3_mesh*pi5_1/60.*tmp*dline_1(:,j)
        enddo
        if (lfirst.and.ldt) then
          advec_hypermesh_ss=chi_hyper3_mesh*pi5_1*sqrt(dxyz_2)
          advec2_hypermesh=advec2_hypermesh+advec_hypermesh_ss**2
        endif
        if (headtt) print*,'denergy_dt: chi_hyper3_mesh=', chi_hyper3_mesh
      endif
!
      if (lheatc_hyper3_polar) then
        do j=1,3
          call der6(f,ilnTT,tmp,j,IGNOREDX=.true.)
          if (.not.ltemperature_nolog) tmp=tmp*p%TT1
          thdiff = thdiff + chi_hyper3*pi4_1*tmp*dline_1(:,j)**2
        enddo
        if (lfirst.and.ldt) &
             diffus_chi3=diffus_chi3+chi_hyper3*pi4_1*dxmin_pencil**4
        if (headtt) print*,'denergy_dt: chi_hyper3=', chi_hyper3
      endif
!
!  Shock diffusion.
!
      if (lheatc_shock) call calc_heatcond_shock(df,p)
!
!  Entry possibility for "personal" entries.
!  In that case you'd need to provide your own "special" routine.
!
      if (lspecial) call special_calc_energy(f,df,p)
!
!  Add thermal diffusion to temperature equation.
!
      if (ltemperature_nolog) then
        df(l1:l2,m,n,iTT)   = df(l1:l2,m,n,iTT)   + thdiff
      else
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + thdiff
      endif
!
!  Boussinesq approximation: - u.grad T_0 added
!
      if (lboussinesq) then
        if (lsphere_in_a_box) then
          df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) - p%r_mn*beta_bouss*( &
          f(l1:l2,m,n,iux)*p%evr(:,1)+f(l1:l2,m,n,iuy)*p%evr(:,2)+    &
          f(l1:l2,m,n,iuz)*p%evr(:,3))
        else
!
! background temperature gradient in z direction
!
          df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) - beta_bouss*f(l1:l2,m,n,iuz)
        endif
      endif
!
! Store  pressure gradient in f-array if requested
!
      if (lpres_grad) then
        f(l1:l2,m,n,igpx) = -p%fpres(:,1)*f(l1:l2,m,n,irho)
        f(l1:l2,m,n,igpy) = -p%fpres(:,2)*f(l1:l2,m,n,irho)
      endif
!
!  Information on the timescales.
!
      if (lfirst.and.ldt) then

        maxdiffus=max(maxdiffus,diffus_chi)
        maxdiffus3=max(maxdiffus3,diffus_chi3)
 
        if (headtt.or.ldebug) then
          print*, 'denergy_dt: max(diffus_chi ) =', maxval(diffus_chi)
          print*, 'denergy_dt: max(diffus_chi3) =', maxval(diffus_chi3)
        endif
      endif
!
!  Apply border profile
!
      if (lborder_profiles) call set_border_entropy(f,df,p)

      call calc_diagnostics_energy(p)
!
    endsubroutine denergy_dt
!***********************************************************************
    subroutine calc_diagnostics_energy(p)

      use Slices_methods, only: store_slices

      type(pencil_case) :: p

      call calc_2d_diagnostics_energy(p)
      call calc_1d_diagnostics_energy(p)
      call calc_0d_diagnostics_energy(p)

      if (lvideo.and.lfirst) then
        if (ivid_pp/=0) call store_slices(p%pp,pp_xy,pp_xz,pp_yz,pp_xy2,pp_xy3,pp_xy4,pp_xz2)
      endif
!
    endsubroutine calc_diagnostics_energy
!***********************************************************************
    subroutine calc_2d_diagnostics_energy(p)
!
!  2-D averages.
!
      use Diagnostics

      type(pencil_case) :: p      
!
      real, dimension (nx) :: NN2_resp, logT_resp, resp, resp2

      if (l2davgfirst) then

        call zsum_mn_name_xy(p%TT,idiag_TTmxy)
        call ysum_mn_name_xz(p%TT,idiag_TTmxz)
!
!   calculating emission with tabulated response files
!
        if (idiag_EmAIA94mxy/=0  .or. idiag_EmAIA94mxz/=0  .or. &
            idiag_EmAIA131mxy/=0 .or. idiag_EmAIA131mxz/=0 .or. &
            idiag_EmAIA171mxy/=0 .or. idiag_EmAIA171mxz/=0 .or. &
            idiag_EmAIA193mxy/=0 .or. idiag_EmAIA193mxz/=0 .or. &
            idiag_EmAIA211mxy/=0 .or. idiag_EmAIA211mxz/=0 .or. &
            idiag_EmAIA304mxy/=0 .or. idiag_EmAIA304mxz/=0 .or. &
            idiag_EmAIA335mxy/=0 .or. idiag_EmAIA335mxz/=0 .or. &
            idiag_EmXRTmxy/=0    .or. idiag_EmXRTmxz/=0        ) then
!
! particle density squared in SI
!
          NN2_resp = (6./7.*p%rho/m_p*unit_density/unit_mass)**2.
!
! log10 temperature in SI (K)
!
          logT_resp =log10(p%TT*unit_temperature)
!
          if (idiag_EmAIA94mxy/=0 .or. idiag_EmAIA94mxz/=0) then
            call get_AIA_tab_resp('094',logT_resp,NN2_resp,resp)
            resp = resp*unit_length/1e6
            if (idiag_EmAIA94mxy/=0) then
              resp2 = resp*Lz
              call zsum_mn_name_xy(resp2,idiag_EmAIA94mxy)
            endif
            if (idiag_EmAIA94mxz/=0) then
              resp2 = resp*Ly
              call ysum_mn_name_xz(resp2,idiag_EmAIA94mxz)
            endif
          endif
!
          if (idiag_EmAIA131mxy/=0 .or. idiag_EmAIA131mxz/=0) then
            call get_AIA_tab_resp('131',logT_resp,NN2_resp,resp)
            resp = resp*unit_length/1e6
            if (idiag_EmAIA131mxy/=0) then
              resp2 = resp*Lz
              call zsum_mn_name_xy(resp2,idiag_EmAIA131mxy)
            endif
            if (idiag_EmAIA131mxz/=0) then
              resp2 = resp*Ly
              call ysum_mn_name_xz(resp2,idiag_EmAIA131mxz)
            endif
          endif
!
          if (idiag_EmAIA171mxy/=0 .or. idiag_EmAIA171mxz/=0) then
            call get_AIA_tab_resp('171',logT_resp,NN2_resp,resp)
            resp = resp*unit_length/1e6
            if (idiag_EmAIA171mxy/=0) then
              resp2 = resp*Lz
              call zsum_mn_name_xy(resp2,idiag_EmAIA171mxy)
            endif
            if (idiag_EmAIA171mxz/=0) then
              resp2 = resp*Ly
              call ysum_mn_name_xz(resp2,idiag_EmAIA171mxz)
            endif
          endif
!
          if (idiag_EmAIA193mxy/=0 .or. idiag_EmAIA193mxz/=0) then
            call get_AIA_tab_resp('193',logT_resp,NN2_resp,resp)
            resp = resp*unit_length/1e6
            if (idiag_EmAIA193mxy/=0) then
              resp2 = resp*Lz
              call zsum_mn_name_xy(resp2,idiag_EmAIA193mxy)
            endif
            if (idiag_EmAIA193mxz/=0) then
              resp2 = resp*Ly
              call ysum_mn_name_xz(resp2,idiag_EmAIA193mxz)
            endif
          endif
!
          if (idiag_EmAIA211mxy/=0 .or. idiag_EmAIA211mxz/=0) then
            call get_AIA_tab_resp('211',logT_resp,NN2_resp,resp)
              resp = resp*unit_length/1e6
            if (idiag_EmAIA211mxy/=0) then
              resp2 = resp*Lz
              call zsum_mn_name_xy(resp2,idiag_EmAIA211mxy)
            endif
            if (idiag_EmAIA211mxz/=0) then
              resp2 = resp*Ly
              call ysum_mn_name_xz(resp2,idiag_EmAIA211mxz)
            endif
          endif
!
          if (idiag_EmAIA304mxy/=0 .or. idiag_EmAIA304mxz/=0) then
            call get_AIA_tab_resp('304',logT_resp,NN2_resp,resp)
            resp = resp*unit_length/1e6
            if (idiag_EmAIA304mxy/=0) then
              resp2 = resp*Lz
              call zsum_mn_name_xy(resp2,idiag_EmAIA304mxy)
            endif
            if (idiag_EmAIA304mxz/=0) then
              resp2 = resp*Ly
              call ysum_mn_name_xz(resp2,idiag_EmAIA304mxz)
            endif
          endif
!
          if (idiag_EmAIA335mxy/=0 .or. idiag_EmAIA335mxz/=0) then
            call get_AIA_tab_resp('335',logT_resp,NN2_resp,resp)
            resp = resp*unit_length/1e6
            if (idiag_EmAIA335mxy/=0) then
              resp2 = resp*Lz
              call zsum_mn_name_xy(resp2,idiag_EmAIA335mxy)
            endif
            if (idiag_EmAIA335mxz/=0) then
              resp2 = resp*Ly
              call ysum_mn_name_xz(resp2,idiag_EmAIA335mxz)
            endif
          endif
!
          if (idiag_EmXRTmxy/=0 .or. idiag_EmXRTmxz/=0) then
            call get_AIA_tab_resp('XRT',logT_resp,NN2_resp,resp)
            resp = resp*unit_length/1e6
            if (idiag_EmXRTmxy/=0) then
              resp2 = resp*Lz
              call zsum_mn_name_xy(resp2,idiag_EmXRTmxy)
            endif
            if (idiag_EmXRTmxz/=0) then
              resp2 = resp*Ly
              call ysum_mn_name_xz(resp2,idiag_EmXRTmxz)
            endif
          endif
!
        endif

      endif
    endsubroutine calc_2d_diagnostics_energy
!***********************************************************************
    subroutine calc_0d_diagnostics_energy(p)
!
!  Calculate temperature related 0D-diagnostics.
!
      use Diagnostics
      use Sub, only: dot2, dot, dot_mn

      type(pencil_case) :: p

      real, dimension(nx) :: tmp
      real :: fradtop, fradbot

      if (ldiagnos) then
!
        call sum_mn_name(p%TT,idiag_TTm)
        if (idiag_TTzmask/=0) call sum_mn_name(p%TT*zmask_temp(n-n1+1),idiag_TTzmask)
        call max_mn_name(p%TT,idiag_TTmax)
        if (idiag_TTmin/=0) call max_mn_name(-p%TT,idiag_TTmin,lneg=.true.)
        call sum_mn_name(p%ss,idiag_ssm)
        call sum_mn_name(p%ee,idiag_eem)
        call sum_mn_name(p%pp,idiag_ppm)
        if (idiag_ethm/=0)  call sum_mn_name(p%rho*p%ee,idiag_ethm)
        if (idiag_ethtot/=0) call integrate_mn_name(p%rho*p%ee,idiag_ethtot)
        call sum_mn_name(p%cs2,idiag_csm,lsqrt=.true.)
        call max_mn_name(p%cs2,idiag_csmax,lsqrt=.true.)
        if (idiag_TugTm/=0) call sum_mn_name(p%TT*p%ugTT,idiag_TugTm)
        if (idiag_Trms/=0)  call sum_mn_name(p%TT**2,idiag_Trms,lsqrt=.true.)
        if (idiag_TT2m/=0)  call sum_mn_name(p%TT**2,idiag_TT2m)
        call sum_mn_name(p%uu(:,1)*p%TT,idiag_uxTm)
        call sum_mn_name(p%uu(:,2)*p%TT,idiag_uyTm)
        call sum_mn_name(p%uu(:,3)*p%TT,idiag_uzTm)
        if (idiag_Tugux_uxugTm/=0) call sum_mn_name(p%TT*p%ugu(:,1)+p%uu(:,1)*p%ugTT,idiag_Tugux_uxugTm)
        if (idiag_Tuguy_uyugTm/=0) call sum_mn_name(p%TT*p%ugu(:,2)+p%uu(:,2)*p%ugTT,idiag_Tuguy_uyugTm)
        if (idiag_Tuguz_uzugTm/=0) call sum_mn_name(p%TT*p%ugu(:,3)+p%uu(:,3)*p%ugTT,idiag_Tuguz_uzugTm)
        if (idiag_Tdxpm/=0)  call sum_mn_name(p%TT*p%fpres(:,1),idiag_Tdxpm)
        if (idiag_Tdypm/=0)  call sum_mn_name(p%TT*p%fpres(:,2),idiag_Tdypm)
        if (idiag_Tdzpm/=0)  call sum_mn_name(p%TT*p%fpres(:,3),idiag_Tdzpm)
!
        if (idiag_gT2m/=0) then
          call dot(p%gTT,p%gTT,tmp)
          call sum_mn_name(tmp,idiag_gT2m)
        endif
!
        if (idiag_guxgTm/=0) then
          call dot_mn(p%gTT(:,:),p%uij(:,1,:),tmp)
          call sum_mn_name(tmp,idiag_guxgTm)
        endif
!
        if (idiag_guygTm/=0) then
          call dot_mn(p%gTT(:,:),p%uij(:,2,:),tmp)
          call sum_mn_name(tmp,idiag_guygTm)
        endif
!
        if (idiag_guzgTm/=0) then
          call dot_mn(p%gTT(:,:),p%uij(:,3,:),tmp)
          call sum_mn_name(tmp,idiag_guzgTm)
        endif
!
        if (idiag_dtc/=0) &
          call max_mn_name(sqrt(advec_cs2)/cdt,idiag_dtc,l_dt=.true.)

        if (idiag_gTmax/=0) then
          call dot2(p%glnTT,tmp)
          call max_mn_name(p%TT*sqrt(tmp),idiag_gTmax)
        endif
!
        if (idiag_fradtop/=0.and.llast_proc_z.and.n==n2) then
          if (lADI) then
            call heatcond_TT(p%TT,hcond)
          else
            hcond=hcond0
          endif
          fradtop=sum(-hcond*p%TT*p%glnTT(:,3)*dsurfxy)
          call surf_mn_name(fradtop,idiag_fradtop)
        endif
!
        if (idiag_fradbot/=0.and.lfirst_proc_z.and.n==n1) then
          if (lADI) then
            call heatcond_TT(p%TT,hcond)
          else
            hcond=hcond0
          endif
          fradbot=sum(-hcond*p%TT*p%glnTT(:,3)*dsurfxy)
          call surf_mn_name(fradbot,idiag_fradbot)
        endif
!
      endif

    endsubroutine calc_0d_diagnostics_energy
!***********************************************************************
    subroutine calc_1d_diagnostics_energy(p)
!
!  1-D averages.
!
      use Diagnostics

      type(pencil_case) :: p

      if (l1davgfirst) then

        if (idiag_fradmz/=0) call xysum_mn_name_z(-hcond0*p%TT*p%glnTT(:,3),idiag_fradmz)
        if (idiag_fconvmz/=0) call xysum_mn_name_z(p%cp*p%rho*p%uu(:,3)*p%TT,idiag_fconvmz)
        call yzsum_mn_name_x(p%pp,idiag_ppmx)
        call xzsum_mn_name_y(p%pp,idiag_ppmy)
        call xysum_mn_name_z(p%pp,idiag_ppmz)
        call yzsum_mn_name_x(p%TT,idiag_TTmx)
        call xzsum_mn_name_y(p%TT,idiag_TTmy)
        call xysum_mn_name_z(p%TT,idiag_TTmz)
        if (idiag_ppuzmz/=0) call xysum_mn_name_z(p%pp*p%uu(:,3),idiag_ppuzmz)
        if (idiag_ethmz/=0) call xysum_mn_name_z(p%rho*p%ee,idiag_ethmz)
        if (idiag_ethuxmx/=0) call yzsum_mn_name_x(p%rho*p%ee*p%uu(:,1),idiag_ethuxmx)
        if (idiag_ethuxmz/=0) call xysum_mn_name_z(p%rho*p%ee*p%uu(:,1),idiag_ethuxmz)
        if (idiag_ethuymz/=0) call xysum_mn_name_z(p%rho*p%ee*p%uu(:,2),idiag_ethuymz)
        if (idiag_ethuzmz/=0) call xysum_mn_name_z(p%rho*p%ee*p%uu(:,3),idiag_ethuzmz)
        call xysum_mn_name_z(p%fpres(:,1),idiag_fpresxmz)
        call xysum_mn_name_z(p%fpres(:,2),idiag_fpresymz)
        call xysum_mn_name_z(p%fpres(:,3),idiag_fpreszmz)
!
        if (idiag_TT2mz/=0) call xysum_mn_name_z(p%TT**2,idiag_TT2mz)
        if (idiag_uxTmz/=0) call xysum_mn_name_z(p%uu(:,1)*p%TT,idiag_uxTmz)
        if (idiag_uyTmz/=0) call xysum_mn_name_z(p%uu(:,2)*p%TT,idiag_uyTmz)
        if (idiag_uzTmz/=0) call xysum_mn_name_z(p%uu(:,3)*p%TT,idiag_uzTmz)
!
      endif

    endsubroutine calc_1d_diagnostics_energy
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
        if (ltemperature_nolog) then
          f_target=0.0
        else
          f_target=1.0
        endif
      case ('constant')
        if (ltemperature_nolog) then
          f_target=TT_const
        else
          f_target=lnTT_const
        endif
      case ('initial-condition')
        if (ltemperature_nolog) then
          call set_border_initcond(f,iTT,f_target)
        else
          call set_border_initcond(f,ilnTT,f_target)
        endif
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
        if (ltemperature_nolog) then
          call border_driving(f,df,p,f_target,iTT)
        else
          call border_driving(f,df,p,f_target,ilnTT)
        endif
      endif
!
    endsubroutine set_border_entropy
!***********************************************************************
    subroutine energy_after_boundary(f)
!
!  Calculation of mean quantities.
!
!  17-apr-15/MR: coded
!
      use Sub, only: grad, finalize_aver
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      real :: fact
      real, dimension (nx,3):: gradTT
      real, dimension (nx)  :: temp
!
      integer :: nl
!
!  Calculate mean of temperature and its gradient.
!
      if (lcalc_TTmean) then
!
        fact=1./nxygrid
        do n=n1,n2
!
          nl = n-n1+1
          TTmz(nl)=0.; gTTmz(nl)=0.
!
          if (ltemperature_nolog) then
            do m=m1,m2
              TTmz(nl)=TTmz(nl)+sum(f(l1:l2,m,n,iTT))
              call grad(f,iTT,gradTT)
              gTTmz(nl)=gTTmz(nl)+sum(gradTT(:,3))
            enddo
          else
            do m=m1,m2
              temp = exp(f(l1:l2,m,n,ilnTT))
              TTmz(nl)=TTmz(nl)+sum(temp)
              call grad(f,ilnTT,gradTT)
              gTTmz(nl)=gTTmz(nl)+sum(gradTT(:,3)*temp)
            enddo
          endif
!
        enddo
!
        call finalize_aver(nprocxy,12,TTmz)
        call finalize_aver(nprocxy,12,gTTmz)
!
        TTmz  = fact*TTmz
        gTTmz = fact*gTTmz    ! simpler by deriving TTmz!!
!
      endif
!
      if (lenergy_slope_limited) &
        call fatal_error('energy_after_boundary', &
                         'Slope-limited diffusion not implemented')

    endsubroutine energy_after_boundary
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
      intent(inout) :: df
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
      if (.not. leos) &
        call fatal_error('rad_equil', &
            'EOS=noeos, but radiative equilibrium requires an EQUATION OF STATE')

      if (.not. ltemperature_nolog) &
          call fatal_error('temperature_idealgas', &
                           'rad_equil not implemented for lnTT')
      if (lroot) print*,'init_energy: rad_equil for kappa-mechanism pb'
!
      if (nzgrid == 1) &
          call fatal_error ('rad_equil', "not implemented for nzgrid=1")
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
!  Fill in the density and temperature f-arrays and each z-processor
!  writes its own setup in the file data/proc#/setup.dat
!
      open(unit=11,file=trim(directory)//'/setup.dat')
      write(11,'(5a14)') 'z','rho','temp','ss','hcond'
      do n=1,nz
        iz=ipz*nz+n
        f(:,:,nghost+n,ilnrho)=lnrho(iz)
        f(:,:,nghost+n,ilnTT)=temp(iz)
        call eoscalc(ilnrho_TT,lnrho(iz),temp(iz),ss=ss)
        call heatcond_TT(temp(iz), hcond)
        write(11,'(5e14.5)') z(nghost+n),exp(lnrho(iz)),temp(iz),ss,hcond
      enddo
      close(11)
!
!  Initialize cs2bot by taking into account the new bottom value of temperature
!  Note: cs2top=cs20 already defined in eos_idealgas.
!
      cs2bot=gamma_m1*temp(1)
      print*,'cs2top, cs2bot=', cs2top, cs2bot
!
    endsubroutine rad_equil
!***********************************************************************
    subroutine calc_heat_cool(f,df,p)
!
      use Diagnostics, only: sum_lim_mn_name
      use EquationOfState, only: cs20
      use Sub, only: step
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: tau, cooling, kappa, a1, a3, prof, heat
      real :: a2, kappa0, kappa0_cgs
!
!  Initialize
!
      intent(in) :: p
      intent(inout) :: df
!
      if (.not. leos) &
        call fatal_error('calc_heat_cool', &
            'EOS=noeos, but heatng/cooling requires an EQUATION OF STATE')

      if (headtt) print*,'enter calc_heat_cool', rcool, wcool, cool, cs20
!
      if (lgravr) then
        if (lboussinesq) then
          prof = step(p%r_mn,r_ext,wcool)
          heat = -cool*prof
          prof = 1.-step(p%r_mn,r_int,wcool)
          heat = heat-cool*prof
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + heat*f(l1:l2,m,n,ilnTT)
        else
          ! 2-D heating/cooling profiles
          prof = exp(-0.5*(p%r_mn/wheat)**2) * (2*pi*wheat**2)**(-1.)
          heat = luminosity*prof
          prof = step(p%r_mn,rcool,wcool)
          heat = heat - cool*prof*(p%cs2-cs20)/cs20
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + p%cv1*p%TT1*heat
        endif
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
          tau = 0.  ! to avoid compiler warnings
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
!  Calculate the radiative diffusion term for constant chi:
!  lnTT version: cp*chi*Div(rho*TT*glnTT)/(rho*cv*TT)
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
      ! [PAB]: Is the following correct for ltemperature_nolog=T?
      ! [PAB]: Should we then not use "iTT" instead of "ilnTT"?
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
    subroutine calc_heatcond_cubicstepchi(df,p)
!
!  Calculate the radiative diffusion term for chi with cubic step profile:
!
!  lnTT version: cp*chi*Div(rho*TT*glnTT)/(rho*cv*TT)
!           = gamma*[chi*((glnrho+glnTT).glnTT+g2lnTT) + glnTT.gradchi]
!
!    TT version: cp*chi*Div(rho*gTT)/(rho*cv)
!           = gamma*[chi*(glnrho.gTT+g2TT) + gT.gradchi]
!
!  21-oct-18/joern: adapted from calc_heatcond_constchi
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot, cubic_step, cubic_der_step
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: g2, chi_z, gradchi_z
!
      intent(in) :: p
      intent(inout) :: df
!
      if (chi_zwidth == 0.) chi_zwidth = 5.*dz
      chi_z= chi + chi*(chi_jump-1.)*cubic_step(p%z_mn,chi_z0,-chi_zwidth)
!
      gradchi_z=chi*(chi_jump-1.)*cubic_der_step(p%z_mn,chi_z0,-chi_zwidth)
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
      if (ltemperature_nolog) then
        df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) + gamma*(chi_z*g2 + gradchi_z*p%gTT(:,3))
      else
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + gamma*(chi_z*g2 + gradchi_z*p%glnTT(:,3))
      endif
!
!  Check maximum diffusion from thermal diffusion.
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chi_z*dxyz_2
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_cubicstepchi
!***********************************************************************
    subroutine calc_heatcond_constK(df,p)
!
!  Calculate the radiative diffusion term for constant K:
!
!  lnTT version: gamma*K/rho/TT/cp*div(TT*grad lnTT)
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
!  Note that rho does not in general need to be unity, even with Boussinesq.
!
      if (ldensity.or.lboussinesq.or.lanelastic) then
        chix=p%rho1*hcond0*p%cp1
      else
        chix=hcond0*p%cp1
      endif
      chix = gamma*chix
!
      if (ltemperature_nolog) then
        df(l1:l2,m,n,iTT)   = df(l1:l2,m,n,iTT)   + chix*p%del2TT
      else
        call dot(p%glnTT,p%glnTT,g2)
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + chix*(g2 + p%del2lnTT)
      endif
!
!  Check maximum diffusion from thermal diffusion.
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+chix*dxyz_2
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_constK
!***********************************************************************
    subroutine calc_heatcond_Ktherm(df,p)
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
      use EquationOfState, only: gamma, rho0
      use Sub, only: dot
!
      real, dimension(mx,my,mz,mvar) :: df
      type (pencil_case)  :: p
      real, dimension(nx) :: g2, chix,hcondTT
!
      intent(in) :: p
      intent(inout) :: df
!
      hcondTT=hcond0*sqrt(exp(p%lnTT))
!
!  Add heat conduction to RHS of temperature equation.
!  If ldensity=F, we need to divide by rho0, which can be /= 1.
!
      if (ldensity) then
        chix=p%rho1*hcondTT*p%cp1
      else
        chix=hcondTT*p%cp1/rho0
      endif
      chix = gamma*chix
!
      if (ltemperature_nolog) then
        df(l1:l2,m,n,iTT)   = df(l1:l2,m,n,iTT)   + chix*p%del2TT
      else
        call dot(p%glnTT,p%glnTT,g2)
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + chix*(g2 + p%del2lnTT)
      endif
!
!  Check maximum diffusion from thermal diffusion.
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+chix*dxyz_2
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_Ktherm
!***********************************************************************
    subroutine calc_heatcond_kramers(df,p)
!
!  Heat conduction using Kramers' opacity law
!
!  23-feb-11/pete: coded
!  24-aug-15/MR: bounds for chi introduced
!
      use Diagnostics
      use Debug_IO, only: output_pencil
      use Sub, only: dot
      use General, only: notanumber
      use EquationOfState, only: gamma
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: p
      intent(inout) :: df
!
      real, dimension(nx) :: thdiff,Krho1,g2
!
!  Diffusion of the form
!      cv*rho*D log(T)/Dt = ... + nab.(K*gradT)/T ,
!  where
!      K = K_0*(T**6.5/rho**2)**n.
!  In reality n=1, but we may need to use n\=1 for numerical reasons.
!
      Krho1 = hcond0_kramers*p%rho1**(2.*nkramers+1.)*p%TT**(6.5*nkramers)   ! = K/rho
      !Krho1 = hcond0_kramers*exp(-p%lnrho*(2.*nkramers+1.)+p%lnTT*(6.5*nkramers))   ! = K/rho
      !if (chimax_kramers>0.) &
      !  Krho1 = max(min(Krho1,chimax_kramers/p%cp1),chimin_kramers/p%cp1)
      call dot(-2.*nkramers*p%glnrho+(6.5*nkramers+1)*p%glnTT,p%glnTT,g2)
      thdiff = Krho1*(p%del2lnTT+g2)

      df(l1:l2,m,n,ilntt) = df(l1:l2,m,n,ilntt) + thdiff

      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+(gamma*p%cp1*Krho1)*dxyz_2
      endif

    endsubroutine calc_heatcond_kramers
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
      real, dimension (nx)   :: dhcond, g1, chix
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
      real, dimension(nx) :: g2, dhcond, chix
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
!  25-aug-09/bing: moved from denergy_dt to here
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
    subroutine get_AIA_tab_resp(channel,lgT,part_den2,respon)
!
!  Determine the indices for the temperature used in the tabulated
!  response functions of AIA on SDO.
!
!  8-may-19/joern: coded
!
      real, dimension  (nx) :: respon
      real, dimension  (nx) :: part_den2
      real, dimension  (nx) :: lgT
!
      real                  :: a0=1., a1=1., a2=1., a3=0., a4=0., a5=0.
!
      character (len=3) :: channel
!
!  the response functions are approximated by gaussian curves
!  using the amplitude (a0), center (a1), width (the standard deviation) (a2) of gaussian
!  following:
!        func(x) = a0 * exp((x-a1)**2/2*a2**2)
!
   select case (trim(channel))
!
!   AIA channels
!
     case ('094')
       a0=  3.68132d-27
       a1=      6.86406
       a2=     0.171703
!
     case ('131')
       a0=  5.69305d-26
       a1=      5.75804
       a2=     0.225604
!
     case ('171')
       a0=  1.29671d-24
       a1=      5.93474
       a2=     0.173903
!
     case ('193')
       a0=  5.42492d-25
       a1=      6.17694
       a2=     0.120203
!
     case ('211')
       a0=  1.71252d-25
       a1=      6.26945
       a2=     0.118102
!
     case ('304')
       a0=  4.26987d-25
       a1=      4.92302
       a2=     0.128603
!
     case ('335')
       a0=  6.38811d-27
       a1=      5.33373
       a2=     0.471810
!
!  XRT channel'
!  we use Al-poly
     case ('XRT')
       a0=  2.72928d-25
       a1=      6.93243
       a2=     0.317607
!
     case default
       a0= 1.
       a1= 1.
       a2=0.1
!
  endselect
!
!
!  Construct the reponse function
!
      respon=a0*(exp(-1.*((lgT-a1)**2.)/(2.*a2**2.)))
!
!  two seconds exposure time
!  mulitplying with the particle density squared
!  getting the units right
!
      respon=2*respon*part_den2*1e-4
!
    endsubroutine get_AIA_tab_resp
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
    subroutine rprint_energy(lreset,lwrite)
!
!  Reads and registers print parameters relevant to energy.
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Diagnostics, only: parse_name, set_type
      use FArrayManager, only: farray_index_append
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
        idiag_TTmax=0; idiag_TTzmask=0; idiag_TTmin=0; idiag_TTm=0; idiag_fradtop=0
        idiag_TugTm=0; idiag_Trms=0; idiag_fradbot=0
        idiag_uxTm=0; idiag_uyTm=0; idiag_uzTm=0; idiag_gT2m=0
        idiag_guxgTm=0; idiag_guygTm=0; idiag_guzgTm=0
        idiag_Tugux_uxugTm=0; idiag_Tuguy_uyugTm=0; idiag_Tuguz_uzugTm=0
        idiag_Tdxpm=0; idiag_Tdypm=0; idiag_Tdzpm=0
        idiag_yHmax=0; idiag_yHmin=0; idiag_yHm=0; idiag_gTmax=0
        idiag_ethm=0; idiag_ssm=0; idiag_thcool=0
        idiag_dtchi=0; idiag_dtc=0
        idiag_eem=0; idiag_ppm=0; idiag_csm=0; idiag_csmax=0
        idiag_ppmx=0; idiag_ppmy=0; idiag_ppmz=0; idiag_ppuzmz=0
        idiag_TTmx=0; idiag_TTmy=0; idiag_TTmz=0; idiag_ethuxmx=0
        idiag_TT2mz=0; idiag_uxTmz=0; idiag_uyTmz=0; idiag_uzTmz=0
        idiag_ethmz=0; idiag_ethuxmz=0; idiag_ethuymz=0; idiag_ethuzmz=0
        idiag_TTmxy=0; idiag_TTmxz=0
        idiag_fpresxmz=0; idiag_fpresymz=0; idiag_fpreszmz=0; idiag_fradmz=0
        idiag_ethtot=0; idiag_fconvmz=0; idiag_TT2m=0
        ivid_pp=0
        idiag_EmAIA94mxy=0; idiag_EmAIA94mxz=0
        idiag_EmAIA131mxy=0; idiag_EmAIA131mxz=0
        idiag_EmAIA171mxy=0; idiag_EmAIA171mxz=0
        idiag_EmAIA193mxy=0; idiag_EmAIA193mxz=0
        idiag_EmAIA211mxy=0; idiag_EmAIA211mxz=0
        idiag_EmAIA304mxy=0; idiag_EmAIA304mxz=0
        idiag_EmAIA335mxy=0; idiag_EmAIA335mxz=0
        idiag_EmXRTmxy=0; idiag_EmXRTmxz=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'TTmax',idiag_TTmax)
        call parse_name(iname,cname(iname),cform(iname),'TTzmask',idiag_TTzmask)
        call parse_name(iname,cname(iname),cform(iname),'gTmax',idiag_gTmax)
        call parse_name(iname,cname(iname),cform(iname),'TTmin',idiag_TTmin)
        call parse_name(iname,cname(iname),cform(iname),'TTm',idiag_TTm)
        call parse_name(iname,cname(iname),cform(iname),'TugTm',idiag_TugTm)
        call parse_name(iname,cname(iname),cform(iname),'Trms',idiag_Trms)
        call parse_name(iname,cname(iname),cform(iname),'TT2m',idiag_TT2m)
        call parse_name(iname,cname(iname),cform(iname),'uxTm',idiag_uxTm)
        call parse_name(iname,cname(iname),cform(iname),'uyTm',idiag_uyTm)
        call parse_name(iname,cname(iname),cform(iname),'uzTm',idiag_uzTm)
        call parse_name(iname,cname(iname),cform(iname),'gT2m',idiag_gT2m)
        call parse_name(iname,cname(iname),cform(iname),'guxgTm',idiag_guxgTm)
        call parse_name(iname,cname(iname),cform(iname),'guygTm',idiag_guygTm)
        call parse_name(iname,cname(iname),cform(iname),'guzgTm',idiag_guzgTm)
        call parse_name(iname,cname(iname),cform(iname),'Tugux_uxugT',idiag_Tugux_uxugTm)
        call parse_name(iname,cname(iname),cform(iname),'Tuguy_uyugT',idiag_Tuguy_uyugTm)
        call parse_name(iname,cname(iname),cform(iname),'Tuguz_uzugT',idiag_Tuguz_uzugTm)
        call parse_name(iname,cname(iname),cform(iname),'Tdxpm',idiag_Tdxpm)
        call parse_name(iname,cname(iname),cform(iname),'Tdypm',idiag_Tdypm)
        call parse_name(iname,cname(iname),cform(iname),'Tdzpm',idiag_Tdzpm)
        call parse_name(iname,cname(iname),cform(iname),'fradtop',idiag_fradtop)
        call parse_name(iname,cname(iname),cform(iname),'fradbot',idiag_fradbot)
        call parse_name(iname,cname(iname),cform(iname),'ethm',idiag_ethm)
        call parse_name(iname,cname(iname),cform(iname),'ethtot',idiag_ethtot)
        call parse_name(iname,cname(iname),cform(iname),'ssm',idiag_ssm)
        call parse_name(iname,cname(iname),cform(iname),'dtchi',idiag_dtchi)
        call parse_name(iname,cname(iname),cform(iname),'dtc',idiag_dtc)
        call parse_name(iname,cname(iname),cform(iname),'eem',idiag_eem)
        call parse_name(iname,cname(iname),cform(iname),'ppm',idiag_ppm)
        call parse_name(iname,cname(iname),cform(iname),'csmax',idiag_csmax)
        call parse_name(iname,cname(iname),cform(iname),'thcool',idiag_thcool)
      enddo

      if (idiag_fradbot/=0) call set_type(idiag_fradbot,lsurf=.true.)
      if (idiag_fradtop/=0) call set_type(idiag_fradtop,lsurf=.true.)
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
!
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TT2mz',idiag_TT2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uxTmz',idiag_uxTmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uyTmz',idiag_uyTmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzTmz',idiag_uzTmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fradmz',idiag_fradmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fconvmz',idiag_fconvmz)
!
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'TTmxy', &
            idiag_TTmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'EmAIA94mxy', &
           idiag_EmAIA94mxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'EmAIA131mxy', &
           idiag_EmAIA131mxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'EmAIA171mxy', &
           idiag_EmAIA171mxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'EmAIA193mxy', &
           idiag_EmAIA193mxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'EmAIA211mxy', &
           idiag_EmAIA211mxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'EmAIA304mxy', &
           idiag_EmAIA304mxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'EmAIA335mxy', &
           idiag_EmAIA335mxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'EmXRTmxy', &
           idiag_EmXRTmxy)
      enddo
!
!  Check for those quantities for which we want y-averages.
!
      do inamexz=1,nnamexz
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'TTmxz', &
            idiag_TTmxz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'EmAIA94mxz', &
           idiag_EmAIA94mxz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'EmAIA131mxz', &
           idiag_EmAIA131mxz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'EmAIA171mxz', &
           idiag_EmAIA171mxz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'EmAIA193mxz', &
           idiag_EmAIA193mxz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'EmAIA211mxz', &
           idiag_EmAIA211mxz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'EmAIA304mxz', &
           idiag_EmAIA304mxz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'EmAIA335mxz', &
           idiag_EmAIA335mxz)
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'EmXRTmxz', &
           idiag_EmXRTmxz)
      enddo
!
!  check for those quantities for which we want video slices
!     
      if (lwrite_slices) then 
        where(cnamev=='TT'.or.cnamev=='lnTT') cformv='DEFINED'
        do iname=1,nnamev
          call parse_name(iname,cnamev(iname),cformv(iname),'pp',ivid_pp)
        enddo
      endif
!
!  Write column where which variable is stored.
!
      if (lwr) then
        if (ltemperature_nolog) then
          call farray_index_append('ilnTT', 0)
        else
          call farray_index_append('iTT', 0)
        endif
      endif
!
    endsubroutine rprint_energy
!***********************************************************************
    subroutine get_slices_energy(f,slices)
!
      use Slices_methods, only: assign_slices_scal, process_slices, log2d, exp2d
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
            call assign_slices_scal(slices,f,iTT)
          else
            call assign_slices_scal(slices,f,ilnTT)
            call process_slices(slices,exp2d)
          endif
!  lnTT
        case ('lnTT')
          call assign_slices_scal(slices,f,ilnTT)
          if (iTT>0) call process_slices(slices,log2d)
!  Pressure
        case ('pp')
          call assign_slices_scal(slices,pp_xy,pp_xz,pp_yz,pp_xy2,pp_xy3,pp_xy4,pp_xz2)
!
      endselect
!
    endsubroutine get_slices_energy
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
      if (.not. leos) &
        call fatal_error('single_polytrope', &
            'EOS=noeos, but polytrope requires an EQUATION OF STATE')

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
      if (.not. leos) &
        call fatal_error('piecew_poly', &
            'EOS=noeos, but polytrope requires an EQUATION OF STATE')

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
!  Initialize energy for two superposed polytropes with a central heating
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

      if (.not. leos) &
        call fatal_error('star_heat', &
            'EOS=noeos, but star heating requires an EQUATION OF STATE')
!
!  Define the radial grid r=[0,r_max], luminosity and gravity
!
      r_max=sqrt(xyz1(1)**2+xyz1(2)**2+xyz1(3)**2)
      r(1)=0. ; lumi(1)=0. ; g(1)=0.
      do i=2,nr
        r(i)=r_max*real(i-1)/(nr-1)
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
      do iter=1,10
        rhotop=rt_old+(rt_new-rt_old)/(rb_new-rb_old)*(rbot-rb_old)
!
        crit=abs(rhotop/rt_new-1.)
        if (crit<=1e-4) exit
        call strat_heat(nr, r, lumi, g, hcond, temp, lnrho, rhotop, rhobot)
!
!  Update new estimates.
!
        rt_old=rt_new
        rb_old=rb_new
        rt_new=rhotop
        rb_new=rhobot
      enddo
      print*,'- iteration completed: rhotop,crit=',rhotop,crit
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
          call get_soundspeed(temp(i),cs2)
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
      if (.not. leos) &
        call fatal_error('strat_heat', &
            'EOS=noeos, but stratified heating requires an EQUATION OF STATE')

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
    subroutine impose_energy_floor(f)
!
!  Dummy subroutine; may not be necessary for lnTT
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine impose_energy_floor
!***********************************************************************
    subroutine dynamical_thermal_diffusion(uc)
!
!  Dummy subroutine
!
      real, intent(in) :: uc
!
      call keep_compiler_quiet(uc)
      call fatal_error('dynamical_thermal_diffusion', 'not implemented yet')
!
    endsubroutine dynamical_thermal_diffusion
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
    subroutine expand_shands_energy
!
!  Expands shorthand labels of temperature diagnostics.
!
!  16-may-12/MR: coded
!
      use Diagnostics, only : expand_cname
!
      if (nname>0) then
        call expand_cname(cname,nname,'uuTm','u',.true.)
        call expand_cname(cname,nname,'gugTm','gu',.true.)
        call expand_cname(cname,nname,'Tdpm','Td',.true.)
      endif
!
    endsubroutine expand_shands_energy
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
    subroutine update_char_vel_energy(f)
!
!  Updates characteristic veelocity for slope-limited diffusion.
!
!  9-jul-18/joern: adapted from update_char_vel_energy in entropy.f90
!
      use EquationOfState, only: eoscalc
!      use General, only: staggered_mean_scal
      use General, only: staggered_max_scal
!
      real, dimension(mx,my,mz,mfarray), intent(INOUT) :: f
!
      real, dimension(mx) :: cs2

      if (.not. leos) &
        call fatal_error('update_char_vel_energy', &
            'EOS=noeos, but sound speed requires an EQUATION OF STATE')
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
    endsubroutine update_char_vel_energy
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
