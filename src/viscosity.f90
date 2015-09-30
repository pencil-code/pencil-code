! $Id$
!
!  This modules implements viscous heating and diffusion terms
!  here for cases 1) nu constant, 2) mu = rho.nu 3) constant and
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lviscosity = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED fvisc(3); diffus_total; diffus_total2; diffus_total3
! PENCILS PROVIDED visc_heat; nu; gradnu(3)
!
!***************************************************************
module Viscosity
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'viscosity.h'
!
  integer, parameter :: nvisc_max=4
  character (len=labellen), dimension(nvisc_max) :: ivisc=''
  character (len=labellen) :: lambda_profile='uniform'
  real :: nu=0.0, nu_tdep=0.0, nu_tdep_exponent=0.0, nu_tdep_t0=0.0
  real :: zeta=0.0, nu_mol=0.0, nu_hyper2=0.0, nu_hyper3=0.0
  real :: nu_hyper3_mesh=5.0, nu_shock=0.0,nu_spitzer=0.0
  real :: nu_jump=1.0, xnu=1.0, xnu2=1.0, znu=1.0, widthnu=0.1, widthnu2=0.1
  real :: C_smag=0.0, gamma_smag=0.0, nu_jump2=1.0
  real :: znu_shock=1.0, widthnu_shock=0.1, nu_jump_shock=1.0
  real :: xnu_shock=1.0
  real :: pnlaw=0.0, Lambda_V0=0.,Lambda_V1=0.,Lambda_H1=0.
  real :: Lambda_V0t=0.,Lambda_V1t=0.,Lambda_V0b=0.,Lambda_V1b=0.
  real :: rzero_lambda=impossible,wlambda=0.,rmax_lambda=impossible
  real :: offamp_lambda=1.,r1_lambda=impossible,r2_lambda=impossible
  real :: lambda_jump=0.,roffset_lambda=0.
  real :: PrM_turb=0.0
  real :: meanfield_nuB=0.0
  real :: nu_infinity=0.,nu0=0.,non_newton_lambda=0.,carreau_exponent=0.
  character (len=labellen) :: nnewton_type='none'
  real :: nnewton_tscale=0.0,nnewton_step_width=0.0
  real, dimension(nx) :: xmask_vis=0
  real, dimension(2) :: vis_xaver_range=(/-max_real,max_real/)
  real, dimension(:), pointer :: etat_x, detat_x
  real, dimension(:), pointer :: etat_y, detat_y
  real, dimension(:), pointer :: etat_z, detat_z
  real, dimension(3) :: nu_aniso_hyper3=0.0
  real, dimension(mx) :: LV0_rprof,LV1_rprof,LH1_rprof,der_LV0_rprof,der_LV1_rprof
  real, pointer :: Pr ! get Prandtl number from hydro
  logical :: lvisc_first=.false.
  logical :: lvisc_simplified=.false.
  logical :: lvisc_nu_non_newtonian=.false.
  logical :: lvisc_rho_nu_const=.false.
  logical :: lvisc_rho_nu_const_bulk=.false.
  logical :: lvisc_sqrtrho_nu_const=.false.
  logical :: lvisc_nu_therm=.false.
  logical :: lvisc_mu_therm=.false.
  logical :: lvisc_nu_const=.false.
  logical :: lvisc_nu_tdep=.false.
  logical :: lvisc_nu_prof=.false.
  logical :: lvisc_nu_profx=.false.
  logical :: lvisc_nu_profr=.false.
  logical :: lvisc_nu_profr_powerlaw=.false.
  logical :: lvisc_nu_profr_twosteps=.false.
  logical :: lvisc_nut_from_magnetic=.false.
  logical :: lvisc_nu_shock=.false.
  logical :: lvisc_nu_shock_profz=.false.
  logical :: lvisc_nu_shock_profr=.false.
  logical :: lvisc_shock_simple=.false.
  logical :: lvisc_hyper2_simplified=.false.
  logical :: lvisc_hyper3_simplified=.false.
  logical :: lvisc_hyper3_polar=.false.
  logical :: lvisc_hyper3_mesh=.false.
  logical :: lvisc_hyper3_rho_nu_const=.false.
  logical :: lvisc_hyper3_mu_const_strict=.false.
  logical :: lvisc_hyper3_nu_const_strict=.false.
  logical :: lvisc_hyper3_rho_nu_const_symm=.false.
  logical :: lvisc_hyper3_rho_nu_const_aniso=.false.
  logical :: lvisc_hyper3_nu_const_aniso=.false.
  logical :: lvisc_hyper3_rho_nu_const_bulk=.false.
  logical :: lvisc_hyper3_nu_const=.false.
  logical :: lvisc_smag_simplified=.false.
  logical :: lvisc_smag_cross_simplified=.false.
  logical :: lvisc_snr_damp=.false.
  logical :: lvisc_heat_as_aux=.false.
  logical :: lvisc_mixture=.false.
  logical :: lvisc_spitzer=.false.
  logical :: lvisc_slope_limited=.false.
  logical :: limplicit_viscosity=.false.
  logical :: lmeanfield_nu=.false.
  logical :: lmagfield_nu=.false.
  logical :: llambda_effect=.false.
  logical :: luse_nu_rmn_prof=.false.
  logical, pointer:: lviscosity_heat
  logical :: lKit_Olem
  real :: damp_sound=0.
  real :: h_slope_limited=0.
  character (LEN=labellen) :: islope_limiter=''
!
  namelist /viscosity_run_pars/ &
      limplicit_viscosity, nu, nu_tdep_exponent, nu_tdep_t0, zeta, &
      nu_hyper2, nu_hyper3, ivisc, nu_mol, C_smag, gamma_smag, nu_shock, &
      nu_aniso_hyper3, lvisc_heat_as_aux,nu_jump,znu,xnu,xnu2,widthnu,widthnu2, &
      pnlaw,llambda_effect,Lambda_V0,Lambda_V1,Lambda_H1, nu_hyper3_mesh, &
      lambda_profile,rzero_lambda,wlambda,r1_lambda,r2_lambda,rmax_lambda,&
      offamp_lambda,lambda_jump,lmeanfield_nu,lmagfield_nu,meanfield_nuB, &
      PrM_turb, roffset_lambda, nu_spitzer, nu_jump2,&
      widthnu_shock, znu_shock, xnu_shock, nu_jump_shock, &
      nnewton_type,nu_infinity,nu0,non_newton_lambda,carreau_exponent,&
      nnewton_tscale,nnewton_step_width,lKit_Olem,damp_sound,luse_nu_rmn_prof, &
      lvisc_slope_limited, h_slope_limited, islope_limiter
!
! other variables (needs to be consistent with reset list below)
  integer :: idiag_nu_tdep=0    ! DIAG_DOC: time-dependent viscosity
  integer :: idiag_fviscm=0     ! DIAG_DOC: Mean value of viscous acceleration
  integer :: idiag_fviscmin=0   ! DIAG_DOC: Min value of viscous acceleration
  integer :: idiag_fviscmax=0   ! DIAG_DOC: Max value of viscous acceleration
  integer :: idiag_fviscrmsx=0  ! DIAG_DOC: Rms value of viscous acceleration
                                ! DIAG_DOC: for the vis_xaver_range
  integer :: idiag_num=0        ! DIAG_DOC: Mean value of viscosity
  integer :: idiag_nusmagm=0    ! DIAG_DOC: Mean value of Smagorinsky viscosity
  integer :: idiag_nusmagmin=0  ! DIAG_DOC: Min value of Smagorinsky viscosity
  integer :: idiag_nusmagmax=0  ! DIAG_DOC: Max value of Smagorinsky viscosity
  integer :: idiag_visc_heatm=0 ! DIAG_DOC: Mean value of viscous heating
  integer :: idiag_qfviscm=0    ! DIAG_DOC: $\left<\qv\cdot
                                ! DIAG_DOC: \fv_{\rm visc}\right>$
  integer :: idiag_Sij2m=0      ! DIAG_DOC: $\left<\Strain^2\right>$
  integer :: idiag_epsK=0       ! DIAG_DOC: $\left<2\nu\varrho\Strain^2\right>$
  integer :: idiag_epsK_LES=0   ! DIAG_DOC:
  integer :: idiag_dtnu=0       ! DIAG_DOC: $\delta t/[c_{\delta t,{\rm v}}\,
                                ! DIAG_DOC:   \delta x^2/\nu_{\rm max}]$
                                ! DIAG_DOC:   \quad(time step relative to
                                ! DIAG_DOC:   viscous time step;
                                ! DIAG_DOC:  see \S~\ref{time-step})
  integer :: idiag_meshRemax=0  ! DIAG_DOC: Max mesh Reynolds number
  integer :: idiag_Reshock=0    ! DIAG_DOC: Mesh Reynolds number at shock
  integer :: idiag_nuD2uxbxm=0  ! DIAG_DOC:
  integer :: idiag_nuD2uxbym=0  ! DIAG_DOC:
  integer :: idiag_nuD2uxbzm=0  ! DIAG_DOC:
!
! xy averaged diagnostics given in xyaver.in written every it1d timestep
!
  integer :: idiag_fviscmz=0    ! XYAVG_DOC: $\left<2\nu\varrho u_i
                                ! XYAVG_DOC: \mathcal{S}_{iz} \right>_{xy}$
                                ! XYAVG_DOC: ($z$-component of viscous flux)
  integer :: idiag_epsKmz=0     ! XYAVG_DOC: $\left<2\nu\varrho\Strain^2
                                ! XYAVG_DOC: \right>_{xy}$
!
! yz averaged diagnostics given in yzaver.in written every it1d timestep
!
  integer :: idiag_fviscmx=0    ! YZAVG_DOC: $\left<2\nu\varrho u_i
                                ! YZAVG_DOC: \mathcal{S}_{ix} \right>_{yz}$
                                ! YZAVG_DOC: ($x$-component of viscous flux)
  integer :: idiag_numx=0       ! YZAVG_DOC: $\left< \nu \right>_{yz}$
                                ! YZAVG_DOC: ($yz$-averaged viscosity)
!
! z averaged diagnostics given in zaver.in
!
  integer :: idiag_fviscmxy=0   ! ZAVG_DOC: $\left<2\nu\varrho u_i
                                ! ZAVG_DOC: \mathcal{S}_{ix} \right>_{z}$
                                ! ZAVG_DOC: ($x$-xomponent of viscous flux)
  integer :: idiag_fviscymxy=0  ! ZAVG_DOC: $\left<2\nu\varrho u_i
                                ! ZAVG_DOC: \mathcal{S}_{iy} \right>_{z}$
                                ! ZAVG_DOC: ($y$-xomponent of viscous flux)
!
! Module Variables
!
  real, dimension(mz) :: eth0z = 0.0
!
  contains
!***********************************************************************
    subroutine register_viscosity

    use FArrayManager, only: farray_register_auxiliary
!
!  19-nov-02/tony: coded
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  Needed for flux limited diffusion
!
      if (lvisc_slope_limited) then
        if (iFF_diff==0) then
          call farray_register_auxiliary('Flux_diff',iFF_diff,vector=dimensionality)
          iFF_diff1=iFF_diff; iFF_diff2=iFF_diff+dimensionality-1
        endif
        call farray_register_auxiliary('Div_flux_diff_uu',iFF_div_uu,vector=3)
        iFF_char_c=max(iFF_div_uu+2,iFF_div_ss)
        if (iFF_div_aa>0) iFF_char_c=max(iFF_char_c,iFF_div_aa+2)
      endif
!
    endsubroutine register_viscosity
!***********************************************************************
    subroutine initialize_viscosity()
!
!  20-nov-02/tony: coded
!
      use EquationOfState, only: get_stratz
      use FArrayManager, only: farray_register_auxiliary
      use Mpicomm, only: stop_it
      use SharedVariables, only: put_shared_variable,get_shared_variable
      use Sub, only: write_zprof, step
!
      integer :: i
      real :: nu_shock_jump1
!
!  Default viscosity.
!
      if ( (nu/=0.0).and.(ivisc(1)=='') ) ivisc(1)='nu-const'
!
!  Some viscosity types need the rate-of-strain tensor and grad(lnrho)
!
      lvisc_simplified=.false.
      lvisc_rho_nu_const=.false.
      lvisc_rho_nu_const_bulk=.false.
      lvisc_sqrtrho_nu_const=.false.
      lvisc_nu_therm=.false.
      lvisc_mu_therm=.false.
      lvisc_nu_const=.false.
      lvisc_nu_tdep=.false.
      lvisc_nu_prof=.false.
      lvisc_nu_profx=.false.
      lvisc_nu_profr=.false.
      lvisc_nu_profr_powerlaw=.false.
      lvisc_nu_profr_twosteps=.false.
      lvisc_nut_from_magnetic=.false.
      lvisc_nu_shock=.false.
      lvisc_nu_shock_profz=.false.
      lvisc_nu_shock_profr=.false.
      lvisc_shock_simple=.false.
      lvisc_hyper2_simplified=.false.
      lvisc_hyper3_simplified=.false.
      lvisc_hyper3_polar=.false.
      lvisc_hyper3_mesh=.false.
      lvisc_hyper3_rho_nu_const=.false.
      lvisc_hyper3_rho_nu_const_symm=.false.
      lvisc_hyper3_mu_const_strict=.false.
      lvisc_hyper3_nu_const_strict=.false.
      lvisc_hyper3_rho_nu_const_aniso=.false.
      lvisc_hyper3_nu_const_aniso=.false.
      lvisc_hyper3_rho_nu_const_bulk=.false.
      lvisc_hyper3_nu_const=.false.
      lvisc_smag_simplified=.false.
      lvisc_smag_cross_simplified=.false.
      lvisc_snr_damp=.false.
      lvisc_spitzer=.false.
!
      do i=1,nvisc_max
        select case (ivisc(i))
        case ('nu-simplified','simplified', '0')
          if (lroot) print*,'viscous force: nu*del2v'
          lvisc_simplified=.true.
        case ('nu-non-newtonian')
          if (lroot) print*,'viscous force: div(nu*Sij)'
          lvisc_nu_non_newtonian=.true.
        case ('rho-nu-const','rho_nu-const', '1')
          if (lroot) print*,'viscous force: mu/rho*(del2u+graddivu/3)'
          lvisc_rho_nu_const=.true.
        case ('rho-nu-const-bulk')
          if (lroot) print*,'viscous force: (zeta/rho)*graddivu'
          lvisc_rho_nu_const_bulk=.true.
        case ('sqrtrho-nu-const')
          if (lroot) print*,'viscous force: mu/sqrt(rho)*(del2u+graddivu/3+S.glnrho)'
          if (nu/=0.) lpenc_requested(i_sij)=.true.
          lvisc_sqrtrho_nu_const=.true.
        case ('nu-therm')
          if (lroot) print*,'viscous force: mu*sqrt(TT)*(del2u+graddivu/3+2S.glnrho)'
          if (nu/=0.) lpenc_requested(i_sij)=.true.
          lvisc_nu_therm=.true.
        case ('mu-therm')
          if (lroot) print*,'viscous force: nu*sqrt(TT)/rho*(del2u+graddivu/3)'
          if (nu/=0.) lpenc_requested(i_sij)=.true.
          lvisc_mu_therm=.true.
        case ('nu-const')
          if (lroot) print*,'viscous force: nu*(del2u+graddivu/3+2S.glnrho)'
          if (nu/=0.) lpenc_requested(i_sij)=.true.
!          if (meanfield_nuB/=0.) lpenc_requested(i_b2)=.true.
          lvisc_nu_const=.true.
        case ('nu-tdep')
          if (lroot) print*,'time-dependent nu*(del2u+graddivu/3+2S.glnrho)'
          if (nu/=0.) lpenc_requested(i_sij)=.true.
          lvisc_nu_tdep=.true.
        case ('nu-prof')
          if (lroot) print*,'viscous force with a vertical profile for nu'
          if (nu/=0.) lpenc_requested(i_sij)=.true.
          lvisc_nu_prof=.true.
        case ('nu-profx')
          if (lroot) print*,'viscous force with a horizontal profile for nu'
          if (nu/=0.) lpenc_requested(i_sij)=.true.
          lvisc_nu_profx=.true.
        case ('nu-profr')
          if (lroot) print*,'viscous force with a radial profile for nu'
          if (nu/=0.) lpenc_requested(i_sij)=.true.
          lvisc_nu_profr=.true.
       case ('nu-profr-twosteps')
          if (lroot) print*,'viscous force with a radial profile for nu'
          if (nu/=0.) lpenc_requested(i_sij)=.true.
          lvisc_nu_profr_twosteps=.true.
        case ('nu-profr-powerlaw','powerlaw','power-law')
          if (lroot) print*,'viscous force with a power law profile'
          if (nu/=0.) lpenc_requested(i_sij)=.true.
          lvisc_nu_profr_powerlaw=.true.
        case ('nut-from-magnetic')
          if (lroot) print*,'nut-from-magnetic via shared variables'
          if (PrM_turb/=0.) lpenc_requested(i_sij)=.true.
          lvisc_nut_from_magnetic=.true.
        case ('nu-shock','shock')
          if (lroot) print*,'viscous force: nu_shock*(XXXXXXXXXXX)'
          lvisc_nu_shock=.true.
          if (.not. lshock) &
           call stop_it('initialize_viscosity: shock viscosity'// &
                           ' but module setting SHOCK=noshock')
        case ('nu-shock-profz')
          if (lroot) print*,'viscous force: nu_shock*(XXXXXXXXXXX)  with a vertical profile'
          lvisc_nu_shock_profz=.true.
          if (.not. lshock) &
           call stop_it('initialize_viscosity: shock viscosity'// &
                           ' but module setting SHOCK=noshock')
          case ('nu-shock-profr')
          if (lroot) print*,'viscous force: nu_shock*(XXXXXXXXXXX)  with a radial profile'
          lvisc_nu_shock_profr=.true.
          if (.not. lshock) &
           call stop_it('initialize_viscosity: shock viscosity'// &
                           ' but module setting SHOCK=noshock')
        case ('shock_simple', 'shock-simple')
          if (lroot) print *, 'viscous force: div(nu_shock*grad(uu_i)))'
          lvisc_shock_simple = .true.
          if (.not. lshock) call stop_it('initialize_viscosity: a SHOCK module is required. ')
        case ('hyper2-simplified','hyper2_simplified', 'hyper4')
          if (lroot) print*,'viscous force: nu_hyper*del4v'
          lvisc_hyper2_simplified=.true.
        case ('hyper3-simplified','hyper3_simplified', 'hyper6')
          if (lroot) print*,'viscous force: nu_hyper*del6v'
          lvisc_hyper3_simplified=.true.
        case ('hyper3-cyl','hyper3_cyl','hyper3-sph','hyper3_sph')
          if (lroot) print*,'viscous force: nu_hyper3/pi^4 *(Deltav)^6/Deltaq^2'
          lvisc_hyper3_polar=.true.
        case ('hyper3-mesh','hyper3_mesh')
          if (lroot) print*,'viscous force: nu_hyper3_mesh/pi^5 *(Deltav)^6/Deltaq'
          lvisc_hyper3_mesh=.true.
        case ('hyper3-rho-nu-const','hyper3_rho_nu-const')
          if (lroot) print*,'viscous force: nu_hyper/rho*del6v'
          lvisc_hyper3_rho_nu_const=.true.
        case ('hyper3-rho-nu-const-symm','hyper3_rho_nu-const_symm')
          if (lroot) print*,'viscous force(i): nu_hyper/rho*(del6ui+der5(divu,i))'
          lvisc_hyper3_rho_nu_const_symm=.true.
        case ('hyper3-mu-const-strict','hyper3_mu-const_strict')
          if (lroot) print*, 'viscous force(i): '// &
              'nu_hyper/rho*(del2(del2(del2(u)))+del2(del2(grad(divu))))'
          if (.not.lhyperviscosity_strict) &
               call stop_it('initialize_viscosity: This viscosity type'//&
               ' cannot be used with HYPERVISC_STRICT=nohypervisc_strict')
          lvisc_hyper3_mu_const_strict=.true.
        case ('hyper3-nu-const-strict','hyper3_nu-const_strict')
          if (lroot) print*, 'viscous force(i): 1/rho*div[2*rho*nu_3*S^(3)]'
          if (.not.lhyperviscosity_strict) &
               call stop_it('initialize_viscosity: This viscosity type'//&
               ' cannot be used with HYPERVISC_STRICT=nohypervisc_strict')
           lvisc_hyper3_nu_const_strict=.true.
        case ('hyper3-rho-nu-const-aniso','hyper3_rho_nu-const_aniso')
          if (lroot) print*,&
               'viscous force(i): 1/rho*(nu.del6)ui'
          lvisc_hyper3_rho_nu_const_aniso=.true.
        case ('hyper3-nu-const-aniso','hyper3_nu-const_aniso')
          if (lroot) print*,&
               'viscous force(i): (nu.del6)ui  + ((nu.uij5).glnrho)'
          lpenc_requested(i_uij5)=.true.
          lpenc_requested(i_glnrho)=.true.
          lvisc_hyper3_nu_const_aniso=.true.
        case ('hyper3-rho-nu-const-bulk','hyper3_rho_nu-const_bulk')
          if (lroot) print*,'viscous force: duj/dt = nu_hyper/rho*d6uj/dxj6'
          lvisc_hyper3_rho_nu_const_bulk=.true.
        case ('hyper3-nu-const','hyper3_nu-const','hyper3-nu_const')
          if (lroot) print*,'viscous force: nu*(del6u+S.glnrho)'
          lpenc_requested(i_uij5)=.true.
          lvisc_hyper3_nu_const=.true.
        case ('smagorinsky-simplified','smagorinsky_simplified')
          if (lroot) print*,'viscous force: Smagorinsky_simplified'
          if (lroot) lvisc_LES=.true.
          lpenc_requested(i_sij)=.true.
          lvisc_smag_simplified=.true.
        case ('smagorinsky-cross-simplif','smagorinsky_cross_simplif')
          if (lroot) print*,'viscous force: Smagorinsky_simplified'
          if (lroot) lvisc_LES=.true.
          lvisc_smag_cross_simplified=.true.
!        case ('snr-damp','snr_damp')
!          if (lroot) print*,'viscous force: SNR damping'
!          lvisc_snr_damp=.true.
        case ('nu-mixture')
          if (lroot) print*,'viscous force: nu is calculated for a mixture'
          lvisc_mixture=.true.
        case ('nu-spitzer')
          if (lroot) print*,'viscous force: temperature dependent nu'
          lpenc_requested(i_sij)=.true.
          lvisc_spitzer=.true.
        case ('none',' ')
          ! do nothing
        case default
          if (lroot) print*, 'No such value for ivisc(',i,'): ', trim(ivisc(i))
          call stop_it('calc_viscous_forcing')
        endselect
      enddo
!
!  If we're timestepping, die or warn if the viscosity coefficient that
!  corresponds to the chosen viscosity type is not set.
!
      if (lrun) then
        if ((lvisc_simplified.or.lvisc_rho_nu_const.or. &
             lvisc_sqrtrho_nu_const.or.lvisc_nu_const.or. &
             lvisc_nu_tdep.or.lvisc_nu_therm.or.&
             lvisc_mu_therm).and.nu==0.0) &
            call warning('initialize_viscosity', &
            'Viscosity coefficient nu is zero!')
        if ((lvisc_rho_nu_const_bulk).and.zeta==0.0) &
            call warning('initialize_viscosity', &
            'Viscosity coefficient zeta is zero!')
        if (lvisc_hyper2_simplified.and.nu_hyper2==0.0) &
            call fatal_error('initialize_viscosity', &
            'Viscosity coefficient nu_hyper2 is zero!')
        if ( (lvisc_hyper3_simplified.or.lvisc_hyper3_rho_nu_const.or. &
              lvisc_hyper3_rho_nu_const_bulk.or.lvisc_hyper3_nu_const.or. &
              lvisc_hyper3_rho_nu_const_symm.or. &
              lvisc_hyper3_polar.or.&
              lvisc_hyper3_mu_const_strict .or. &
              lvisc_hyper3_nu_const_strict ).and. &
              nu_hyper3==0.0 ) &
            call fatal_error('initialize_viscosity', &
            'Viscosity coefficient nu_hyper3 is zero!')
        if (lvisc_hyper3_mesh.and.nu_hyper3_mesh==0.0) &
             call fatal_error('initialize_viscosity', &
            'Viscosity coefficient nu_hyper3_mesh is zero!')
        if (lvisc_spitzer.and.nu_spitzer==0.0) &
             call fatal_error('initialize_viscosity', &
            'Viscosity coefficient nu_spitzer is zero!')
        if ( (lvisc_hyper3_rho_nu_const_aniso.or.lvisc_hyper3_nu_const_aniso).and.&
             ((nu_aniso_hyper3(1)==0. .and. nxgrid/=1 ).or. &
              (nu_aniso_hyper3(2)==0. .and. nygrid/=1 ).or. &
              (nu_aniso_hyper3(3)==0. .and. nzgrid/=1 )) ) &
            call fatal_error('initialize_viscosity', &
             'A viscosity coefficient of nu_aniso_hyper3 is zero!')
        if ( (lvisc_smag_simplified.or.lvisc_smag_cross_simplified).and. &
             C_smag==0.0 ) &
            call fatal_error('initialize_viscosity', &
            'Viscosity coefficient C_smag is zero!')
        if (lvisc_nu_shock .and.nu_shock==0.0) &
            call fatal_error('initialize_viscosity', &
            'Viscosity coefficient nu_shock is zero!')
        if (lvisc_nu_shock_profz .and.nu_shock==0.0) &
            call fatal_error('initialize_viscosity', &
            'Viscosity coefficient nu_shock is zero!')
        if (lvisc_nu_shock_profr .and.nu_shock==0.0) &
            call fatal_error('initialize_viscosity', &
            'Viscosity coefficient nu_shock is zero!')
        if (lvisc_shock_simple .and. nu_shock == 0.0) call fatal_error('initialize_viscosity', 'nu_shock is zero. ')
!
!  Dynamical hyper-diffusivity operates only for mesh formulation of hyper-viscosity
!
        if (ldynamical_diffusion.and..not.lvisc_hyper3_mesh) then
          call fatal_error("initialize_viscosity",&
               "Dynamical diffusion requires mesh hyper-diffusion, switch ivisc='hyper3-mesh'")
        endif
!
      endif

      if (lvisc_nu_prof) then
!
!  Write out viscosity z-profile.
!  At present only correct for Cartesian geometry
!
        call write_zprof('visc', nu + (nu*(nu_jump-1.))*step(z(n1:n2),znu,-widthnu))
      endif

      if (lvisc_nu_shock_profz .or. lvisc_nu_shock_profr) then
        nu_shock_jump1=nu_shock*(nu_jump_shock-1.)
!
!  Write out shock viscosity z-profile.
!  At present only correct for Cartesian geometry
!
        if (lvisc_nu_shock_profz) &
          call write_zprof('visc_shock', &
                           nu_shock+nu_shock_jump1*step(z(n1:n2),znu_shock,-widthnu_shock))

!
!  Write out shock viscosity r-profile.
!  At present only correct for spherical and cylindrical geometry.
!
        if (lvisc_nu_shock_profr) &
          call write_zprof('visc_shock', &
                           nu_shock+nu_shock_jump1*step(x(l1:l2),xnu_shock,-widthnu_shock))
      endif
!
!  Register an extra aux slot for dissipation rate if requested (so
!  visc_heat is written sto snapshots and can be easily analyzed later).
!    NB: We are doing this here, rather than in register_viscosity, as the
!  register_XXX routines are called before read_{start,run}pars, so
!  lvisc_heat_as_aux isn't known there. This implies that we need to
!  append the ivisc_heat line to index.pro manually.
!
      if (lvisc_heat_as_aux) then
        call farray_register_auxiliary('visc_heat',ivisc_heat)
!
        if (lroot) then
          open(3,file=trim(datadir)//'/index.pro', POSITION='append')
          write(3,*) 'ivisc_heat=',ivisc_heat
          close(3)
          open(15,FILE=trim(datadir)//'/def_var.pro',position='append')
          write(15,*) 'visc_heat = fltarr(mx,my,mz)*one'
          close(15)
        endif
      endif
!
!  Shared variables.
!
      call put_shared_variable('lvisc_hyper3_nu_const_strict',lvisc_hyper3_nu_const_strict, &
                               caller='initialize_viscosity')
      call put_shared_variable('nu',nu)
      call get_shared_variable('lviscosity_heat',lviscosity_heat)
      call put_shared_variable('llambda_effect',llambda_effect)
!
!  For Boussinesq, density is constant, so must use lvisc_simplified.
!
      if (lrun.and.lboussinesq) then
        if (.not.lvisc_simplified) call fatal_error("initialize_viscosity", &
            "Boussinesq works only if ivisc='simplified' in run.in")
        if (lroot) print*,'viscosity: viscous force nu*del2v with nu=',nu
      endif
!
! DM Note:
! initialization of lambda should come after sharing the llambda_effect
! variable because even if there is no lambda_effect that has to be
! communicated to the boundary conditions too. In case llambda_effect is true
! more variables are shared inside initialize_lambda.
! initialize lambda effect.
!
      if (llambda_effect) call initialize_lambda
!
!  Check for possibility of getting etat profile and gradient from
!  the magnetic meanfield module.
!
      if (PrM_turb/=0.) then
        if (lmagn_mf) then
          call get_shared_variable('etat_x',etat_x)
          call get_shared_variable('etat_y',etat_y)
          call get_shared_variable('etat_z',etat_z)
          call get_shared_variable('detat_x',detat_x)
          call get_shared_variable('detat_y',detat_y)
          call get_shared_variable('detat_z',detat_z)
          do n=n1,n2
            print*,ipz,z(n),etat_z(n),detat_z(n)
          enddo
        endif
      endif
!
!
!  Compute mask for x-averaging where x is in vis_xaver_range.
!  Normalize such that the average over the full domain
!  gives still unity.
!
      if (l1 == l2) then
        xmask_vis = 1.
      else
        where (x(l1:l2) >= vis_xaver_range(1) .and. x(l1:l2) <= vis_xaver_range(2))
          xmask_vis = 1.
        elsewhere
          xmask_vis = 0.
        endwhere
        vis_xaver_range(1) = max(vis_xaver_range(1), xyz0(1))
        vis_xaver_range(2) = min(vis_xaver_range(2), xyz1(1))
        if (lspherical_coords) then
          xmask_vis = xmask_vis * (xyz1(1)**3 - xyz0(1)**3) &
              / (vis_xaver_range(2)**3 - vis_xaver_range(1)**3)
        elseif (lcylindrical_coords) then
          xmask_vis = xmask_vis * (xyz1(1)**2 - xyz0(1)**2) &
              / (vis_xaver_range(2)**2 - vis_xaver_range(1)**2)
        else
          xmask_vis = xmask_vis*Lxyz(1) &
              / (vis_xaver_range(2) - vis_xaver_range(1))
        endif
      endif
!
!  Get background energy stratification, if any.
!
      if (lstratz .and. lthermal_energy) call get_stratz(z, eth0z=eth0z)

      lslope_limit_diff = lslope_limit_diff .or. lvisc_slope_limited
!
!  debug output
!
      if (lroot.and.ip<14) then
        print*,'xmask_vis=',xmask_vis
      endif
!
    endsubroutine initialize_viscosity
!***********************************************************************
    subroutine initialize_lambda
!
! DM 26-05-2010: cut out of intialize viscosity
!
      use Sub, only: step,der_step,stepdown,der_stepdown
      use SharedVariables, only: put_shared_variable,get_shared_variable
!
      if ((Lambda_V0==0).and.(Lambda_V1==0).and.(Lambda_H1==0)) &
        call warning('initialize_lambda', &
            'You have chose llambda_effect=T but, all Lambda coefficients to be zero!')
      if ((Lambda_V0==0).and.((Lambda_V1/=0).or.(Lambda_H1==0))) &
        call warning('initialize_lambda', &
            'Lambda effect: V_zero=0 but V1 or H1 nonzero')
!
! Select the profile of Lambda, default is uniform. At present (May 2010) the
! only other coded profile is radial step.
!
      select case (lambda_profile)
      case ('radial_step_V0')
        if (lroot) print*,'lambda profile radial_step_V0, rzero_lambda, wlambda:',&
            rzero_lambda,wlambda
        LV0_rprof=step(x,rzero_lambda,wlambda)
        der_LV0_rprof=der_step(x,rzero_lambda,wlambda)
        LV1_rprof=1.;LH1_rprof=1.
        der_LV1_rprof=0.
      case ('radial_step')
        if (lroot) print*,'lambda profile radial_step, rzero_lambda, wlambda:',&
            rzero_lambda,wlambda
        LV0_rprof=step(x,rzero_lambda,wlambda)
        LV1_rprof=step(x,rzero_lambda,wlambda)
        LH1_rprof=step(x,rzero_lambda,wlambda)
        der_LV0_rprof=der_step(x,rzero_lambda,wlambda)
        der_LV1_rprof=der_step(x,rzero_lambda,wlambda)
      case ('top_hat')
        if (lroot) print*,'lambda profile top_hat, rzero_lambda, rmax_lambda,roffset_lambda, wlambda:',&
            rzero_lambda,rmax_lambda,roffset_lambda,wlambda
        LV0_rprof=step(x,rzero_lambda+roffset_lambda,wlambda) &
            -step(x,rmax_lambda+roffset_lambda,wlambda)
        LV1_rprof=step(x,rzero_lambda,wlambda)-step(x,rmax_lambda,wlambda)
        LH1_rprof=step(x,rzero_lambda,wlambda)-step(x,rmax_lambda,wlambda)
        der_LV0_rprof=der_step(x,rzero_lambda+roffset_lambda,wlambda) &
            -der_step(x,rmax_lambda+roffset_lambda,wlambda)
        der_LV1_rprof=der_step(x,rzero_lambda,wlambda) &
            -der_step(x,rmax_lambda,wlambda)
      case ('V1H1_roff')
        LV0_rprof=1.;der_LV0_rprof=0.
        LV1_rprof=1.+offamp_lambda*stepdown(x,rmax_lambda,wlambda)
        LH1_rprof=1.+offamp_lambda*stepdown(x,rmax_lambda,wlambda)
        der_LV1_rprof=offamp_lambda*der_stepdown(x,rmax_lambda,wlambda)
        if (lroot) print*,'LV1_rprof',LV1_rprof
        if (lroot) print*,'LH1_rprof',LH1_rprof
      case ('roff')
        LV0_rprof=1.+stepdown(x,rmax_lambda,wlambda)
        der_LV0_rprof=1.+der_stepdown(x,rmax_lambda,wlambda)
        LV1_rprof=1.+stepdown(x,rmax_lambda,wlambda)
        LH1_rprof=1.+stepdown(x,rmax_lambda,wlambda)
        der_LV1_rprof=der_stepdown(x,rmax_lambda,wlambda)
        if (lroot) print*,'LV1_rprof',LV1_rprof
        if (lroot) print*,'LH1_rprof',LH1_rprof
      case ('uniform')
        if (lroot) print*,'lambda profile :',lambda_profile
        LV0_rprof=1.;LV1_rprof=1;LH1_rprof=1.
        der_LV0_rprof=0.;der_LV1_rprof=0
      case ('V0jump')
        if (lroot) print*,'lambda_profile, radial_step, rzero_lambda, wlambda:',&
            lambda_profile,rzero_lambda,wlambda
        LV0_rprof=1.+(lambda_jump/Lambda_V0)*step(x,rzero_lambda,wlambda)
        der_LV0_rprof=(lambda_jump/Lambda_V0)*der_step(x,rzero_lambda,wlambda)
        LV1_rprof=1;LH1_rprof=1.
        der_LV1_rprof=0
      case default
        call fatal_error('initialize_lambda',&
            'default lambda_profile is uniform ! ')
      endselect

      lambda_V0t=lambda_V0*LV0_rprof(nx)
      lambda_V1t=lambda_V1*LV1_rprof(nx)
      lambda_V0b=lambda_V0*LV0_rprof(1)
      lambda_V1b=lambda_V1*LV1_rprof(1)
      call put_shared_variable('Lambda_V0t',Lambda_V0t,caller='initialize_lambda')
      call put_shared_variable('Lambda_V1t',Lambda_V1t)
      call put_shared_variable('Lambda_V0b',Lambda_V0b)
      call put_shared_variable('Lambda_V1b',Lambda_V1b)
      call put_shared_variable('Lambda_H1',Lambda_H1)
      call put_shared_variable('LH1_rprof',LH1_rprof)
!
    endsubroutine initialize_lambda
!***********************************************************************
    subroutine read_viscosity_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=viscosity_run_pars, IOSTAT=iostat)
!
    endsubroutine read_viscosity_run_pars
!***********************************************************************
    subroutine write_viscosity_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=viscosity_run_pars)
!
    endsubroutine write_viscosity_run_pars
!***********************************************************************
    subroutine rprint_viscosity(lreset,lwrite)
!
!  Writes ishock to index.pro file
!
!  24-nov-03/tony: adapted from rprint_ionization
!
      use Diagnostics, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
      integer :: iname,inamex,inamez,ixy
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtnu=0; idiag_nu_LES=0; idiag_Sij2m=0; idiag_epsK=0; idiag_epsK_LES=0
        idiag_visc_heatm=0; idiag_meshRemax=0; idiag_Reshock=0
        idiag_nuD2uxbxm=0; idiag_nuD2uxbym=0; idiag_nuD2uxbzm=0
        idiag_nu_tdep=0; idiag_fviscm=0 ; idiag_fviscrmsx=0
        idiag_fviscmz=0; idiag_fviscmx=0; idiag_fviscmxy=0
        idiag_epsKmz=0; idiag_numx=0; idiag_fviscymxy=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) print*,'rprint_viscosity: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'nu_tdep',idiag_nu_tdep)
        call parse_name(iname,cname(iname),cform(iname),'fviscm',idiag_fviscm)
        call parse_name(iname,cname(iname),cform(iname),'fviscmin',idiag_fviscmin)
        call parse_name(iname,cname(iname),cform(iname),'qfviscm',idiag_qfviscm)
        call parse_name(iname,cname(iname),cform(iname),'fviscmax',idiag_fviscmax)
        call parse_name(iname,cname(iname),cform(iname),'fviscrmsx',idiag_fviscrmsx)
        call parse_name(iname,cname(iname),cform(iname),'num',idiag_num)
        call parse_name(iname,cname(iname),cform(iname),'nusmagm',idiag_nusmagm)
        call parse_name(iname,cname(iname),cform(iname),'nusmagmin',idiag_nusmagmin)
        call parse_name(iname,cname(iname),cform(iname),'nusmagmax',idiag_nusmagmax)
        call parse_name(iname,cname(iname),cform(iname),'dtnu',idiag_dtnu)
        call parse_name(iname,cname(iname),cform(iname),'nu_LES',idiag_nu_LES)
        call parse_name(iname,cname(iname),cform(iname),'visc_heatm', &
            idiag_visc_heatm)
        call parse_name(iname,cname(iname),cform(iname),'Sij2m',idiag_Sij2m)
        call parse_name(iname,cname(iname),cform(iname),'epsK',idiag_epsK)
        call parse_name(iname,cname(iname),cform(iname),'epsK_LES',idiag_epsK_LES)
        call parse_name(iname,cname(iname),cform(iname),'meshRemax',idiag_meshRemax)
        call parse_name(iname,cname(iname),cform(iname),'Reshock',idiag_Reshock)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fviscmz',idiag_fviscmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'epsKmz',idiag_epsKmz)
      enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'fviscmx',idiag_fviscmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'numx',idiag_numx)
      enddo
!
!  Check for those quantities for which we want z-averages
!
      do ixy=1,nnamexy
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'fviscmxy',idiag_fviscmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'fviscymxy',idiag_fviscymxy)
      enddo
!
!  write column where which viscosity variable is stored
!
      if (present(lwrite)) then
        if (lwrite) then
          write(3,*) 'ihypvis=',ihypvis
        endif
      endif
!
    endsubroutine rprint_viscosity
!***********************************************************************
    subroutine pencil_criteria_viscosity()
!
!  All pencils that the Viscosity module depends on are specified here.
!
!  20-11-04/anders: coded
!  18-05-12/MR: request of sij2 for lvisc_simplified.and.lboussinesq added
!                       of graddivu for lboussinesq diasabled
!
      if ((lentropy.or.ltemperature) .and. &
          (lvisc_rho_nu_const .or. lvisc_rho_nu_const_bulk .or. &
           lvisc_sqrtrho_nu_const .or. lvisc_nu_therm .or.&
           lvisc_nu_const .or. lvisc_nu_tdep .or. lvisc_nu_shock .or. &
           lvisc_nu_prof .or. lvisc_nu_profx .or. lvisc_spitzer .or. &
           lvisc_nu_profr .or. lvisc_nu_profr_powerlaw .or. &
           lvisc_nu_profr_twosteps .or. lvisc_nu_shock_profz .or. &
           lvisc_nut_from_magnetic .or. lvisc_nu_shock_profr .or. lvisc_mu_therm))&
           lpenc_requested(i_TT1)=.true.
      if (lvisc_rho_nu_const .or. lvisc_rho_nu_const_bulk .or. &
          lvisc_sqrtrho_nu_const .or. &
          lvisc_nu_const .or. lvisc_nu_tdep .or. lvisc_nu_therm .or. &
          lvisc_nu_prof.or.lvisc_nu_profx.or.lvisc_spitzer .or. &
          lvisc_nu_profr.or.lvisc_nu_profr_powerlaw .or. &
          lvisc_nu_profr_twosteps .or. &
          lvisc_nut_from_magnetic.or.lvisc_mu_therm.or. &
          (lvisc_simplified.and.lboussinesq) ) then
        if (lenergy.and.lviscosity_heat) lpenc_requested(i_sij2)=.true.
        lpenc_requested(i_graddivu)=.true.
      endif
      if (lthermal_energy.and.lviscosity_heat) lpenc_requested(i_rho)=.true.
      if ((lentropy.or.ltemperature).and. &
          (lvisc_nu_therm.or.lvisc_mu_therm.or.lvisc_spitzer)) &
          lpenc_requested(i_lnTT)=.true.
      if (lvisc_smag_simplified .or. lvisc_smag_cross_simplified) &
          lpenc_requested(i_graddivu)=.true.
      if (lvisc_smag_simplified) lpenc_requested(i_sij2)=.true.
      if (lvisc_smag_cross_simplified) lpenc_requested(i_ss12)=.true.
      if (lvisc_nu_prof) lpenc_requested(i_z_mn)=.true.
      if (lvisc_nu_profx) lpenc_requested(i_x_mn)=.true.
      if (lvisc_nu_profr .or. lvisc_nu_profr_twosteps) then
        if (lsphere_in_a_box.or.lspherical_coords) then
          lpenc_requested(i_r_mn)=.true.
        else
          lpenc_requested(i_rcyl_mn)=.true.
          if (lcylinder_in_a_box) lpenc_requested(i_rcyl_mn1)=.true.
        endif
      endif
      if (lvisc_nu_non_newtonian) then
        lpenc_requested(i_nu)=.true.
        lpenc_requested(i_del2u)=.true.
        lpenc_requested(i_sij)=.true.
        lpenc_requested(i_sij2)=.true.
        lpenc_requested(i_uijk)=.true.
      endif
      if (lvisc_nu_profr_powerlaw) lpenc_requested(i_rcyl_mn)=.true.
      if (lvisc_nu_profr_powerlaw.and.luse_nu_rmn_prof) &
          lpenc_requested(i_r_mn)=.true.
      if (lvisc_rho_nu_const .or. &
          lvisc_sqrtrho_nu_const .or. lvisc_nu_const .or. lvisc_nu_tdep .or. &
          lvisc_smag_simplified .or. lvisc_smag_cross_simplified .or. &
          lvisc_nu_prof .or. lvisc_nu_profx .or. lvisc_spitzer .or. &
          lvisc_nu_profr_powerlaw .or. lvisc_nu_profr .or. &
          lvisc_nu_profr_twosteps .or. &
          lvisc_nut_from_magnetic .or. lvisc_nu_therm .or. lvisc_mu_therm) &
          lpenc_requested(i_del2u)=.true.
      if (.not. limplicit_viscosity .and. lvisc_simplified) lpenc_requested(i_del2u)=.true.
      if (lvisc_hyper3_simplified .or. lvisc_hyper3_rho_nu_const .or. &
          lvisc_hyper3_nu_const .or. lvisc_hyper3_rho_nu_const_symm) &
          lpenc_requested(i_del6u)=.true.
      if (lvisc_hyper3_rho_nu_const_symm) then
        lpenc_requested(i_grad5divu)=.true.
        if ((lentropy.or.ltemperature).and.lviscosity_heat) then
          lpenc_requested(i_uij5)=.true.
          lpenc_requested(i_uij)=.true.
        endif
      endif
      if (lvisc_hyper3_rho_nu_const_bulk) lpenc_requested(i_del6u_bulk)=.true.
      if (lvisc_hyper2_simplified) lpenc_requested(i_del4u)=.true.
      if (lvisc_rho_nu_const .or. lvisc_rho_nu_const_bulk .or. &
          lvisc_sqrtrho_nu_const .or. &
          lvisc_hyper3_rho_nu_const .or. lvisc_nu_therm .or.&
          lvisc_hyper3_rho_nu_const_bulk .or. &
          lvisc_hyper3_rho_nu_const_aniso .or. &
          lvisc_smag_simplified .or. lvisc_smag_cross_simplified .or. &
          lvisc_hyper3_rho_nu_const_symm .or. &
          lvisc_hyper3_mu_const_strict .or. lvisc_mu_therm .or. &
          lvisc_spitzer) &
          lpenc_requested(i_rho1)=.true.
!
      if (lvisc_nu_const .or. lvisc_nu_tdep .or. &
          lvisc_nu_prof .or. lvisc_nu_profx .or. &
          lvisc_smag_simplified .or. lvisc_smag_cross_simplified .or. &
          lvisc_nu_profr_powerlaw .or. lvisc_nu_profr .or. &
          lvisc_nu_profr_twosteps .or. lvisc_sqrtrho_nu_const .or. &
          lvisc_nut_from_magnetic.or.lvisc_nu_therm) &
          lpenc_requested(i_sglnrho)=.true.
!
      if (lvisc_spitzer .or. lvisc_mu_therm .or. lvisc_nu_therm) &
          lpenc_requested(i_sglnTT)=.true.
!
      if (lvisc_nu_const .and. lmagfield_nu) lpenc_requested(i_b2)=.true.
      if (lvisc_hyper3_nu_const) lpenc_requested(i_uij5glnrho)=.true.
      if (ldensity.and.lvisc_nu_shock) then
        lpenc_requested(i_graddivu)=.true.
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_gshock)=.true.
        lpenc_requested(i_divu)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
      if (ldensity.and.lvisc_nu_shock_profz) then
        lpenc_requested(i_graddivu)=.true.
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_gshock)=.true.
        lpenc_requested(i_divu)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
      if (ldensity.and.lvisc_nu_shock_profr) then
        lpenc_requested(i_graddivu)=.true.
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_gshock)=.true.
        lpenc_requested(i_divu)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
      shksmp: if (lvisc_shock_simple) then
        lpenc_requested(i_shock) = .true.
        lpenc_requested(i_gshock) = .true.
        lpenc_requested(i_uij) = .true.
        lpenc_requested(i_del2u) = .true.
      endif shksmp
      if (lvisc_slope_limited) then
        lpenc_requested(i_del2u)=.true.
      endif
      if (llambda_effect) then
        lpenc_requested(i_uij)=.true.
        lpenc_requested(i_glnrho)=.true.
        !if (lKit_Olem) lpenc_requested(i_dsdr) = .true.
        !if (lKit_Olem) lpenc_requested(i_gss) = .true.
      endif
      if (lvisc_mixture) then
        lpenc_requested(i_graddivu)=.true.
        lpenc_requested(i_del2u)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_sglnrho)=.true.
        lpenc_requested(i_nu)=.true.
        lpenc_requested(i_gradnu)=.true.
        lpenc_requested(i_sij)=.true.
        if ((lentropy.or.ltemperature).and.lviscosity_heat) &
            lpenc_requested(i_sij2)=.true.
      endif
!
      if (idiag_meshRemax/=0.or.idiag_Reshock/=0) lpenc_diagnos(i_u2)=.true.
      if (idiag_visc_heatm/=0) lpenc_diagnos(i_visc_heat)=.true.
      if (idiag_epsK/=0.or.idiag_epsK_LES/=0.or.idiag_epsKmz/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_sij2)=.true.
      endif
      if (idiag_Sij2m/=0.) lpenc_diagnos(i_sij2)=.true.
      if (idiag_epsK/=0.or.idiag_epsKmz/=0) then
        lpenc_diagnos(i_visc_heat)=.true.
        lpenc_diagnos(i_uu)=.true.
      endif
      if (lvisc_nu_shock.and.idiag_epsK/=0) then
        lpenc_diagnos(i_fvisc)=.true.
        lpenc_diagnos(i_diffus_total)=.true.
        lpenc_diagnos(i_shock)=.true.
        lpenc_diagnos(i_divu)=.true.
        lpenc_diagnos(i_rho)=.true.
      endif
      if (lvisc_nu_shock_profz.and.idiag_epsK/=0) then
        lpenc_diagnos(i_fvisc)=.true.
        lpenc_diagnos(i_diffus_total)=.true.
        lpenc_diagnos(i_shock)=.true.
        lpenc_diagnos(i_divu)=.true.
        lpenc_diagnos(i_rho)=.true.
      endif
      if (lvisc_nu_shock_profr.and.idiag_epsK/=0) then
        lpenc_diagnos(i_fvisc)=.true.
        lpenc_diagnos(i_diffus_total)=.true.
        lpenc_diagnos(i_shock)=.true.
        lpenc_diagnos(i_divu)=.true.
        lpenc_diagnos(i_rho)=.true.
      endif
      if (lvisc_heat_as_aux) then
        lpenc_diagnos(i_visc_heat)=.true.
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_sij2)=.true.
      endif
      if ((idiag_meshRemax/=0.or.idiag_dtnu/=0).and.(lvisc_nu_shock &
           .or.lvisc_nu_shock_profz.or.lvisc_nu_shock_profr)) then
        lpenc_diagnos(i_shock)=.true.
      endif
      if (idiag_qfviscm/=0) then
        lpenc_diagnos(i_curlo)=.true.
        lpenc_diagnos(i_fvisc)=.true.
      endif
      if (idiag_Reshock/=0) lpenc_diagnos(i_shock)=.true.
      if (idiag_fviscmz/=0.or.idiag_fviscmx/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_sij)=.true.
      endif
      if (idiag_numx/=0) then
        lpenc_diagnos(i_nu) = .true.
      endif
      if (idiag_fviscmxy/=0 .or. idiag_fviscymxy/=0) then
        lpenc_diagnos2d(i_rho)=.true.
        lpenc_diagnos2d(i_sij)=.true.
      endif
      if (lboussinesq) lpenc_requested(i_graddivu)=.false.
      if (damp_sound/=0.) lpenc_requested(i_divu)=.true.
!
    endsubroutine pencil_criteria_viscosity
!***********************************************************************
    subroutine pencil_interdep_viscosity(lpencil_in)
!
!  Interdependency among pencils from the Viscosity module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension (npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_viscosity
!***********************************************************************
    subroutine calc_pencils_viscosity(f,p)
!
!  Calculate Viscosity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-nov-04/anders: coded
!  18-may-12/MR: calculation of viscous heat for boussinesq added
!
      use Deriv, only: der5i1j,der6
      use Diagnostics, only: max_mn_name, sum_mn_name
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (nx,3) :: tmp,tmp2,gradnu,sgradnu, gradnu_shock
      real, dimension (nx) :: murho1,zetarho1,muTT,nu_smag,tmp3,tmp4,pnu, pnu_shock
      real, dimension (nx) :: lambda_phi,prof,prof2,derprof,derprof2,qfvisc
      real, dimension (nx) :: gradnu_effective,fac, nu_sld
      real, dimension (nx,3) :: deljskl2,fvisc_nnewton2
!
      integer :: i,j,ju,ii,jj,kk,ll
      real :: nu_shock_jump1
!
      intent(inout) :: f,p
!
!  Viscous force and viscous heating are calculated here (for later use).
!
      p%fvisc=0.0                               !!! not needed
      if (lpencil(i_visc_heat)) p%visc_heat=0.0
      if (lfirst .and. ldt) then
        p%diffus_total=0.0
        p%diffus_total2=0.0
        p%diffus_total3=0.0
      endif
!
!  viscous force: nu*del2v
!  -- not physically correct (no momentum conservation), but
!  numerically easy and in most cases qualitatively OK,
!  for boussinesq (divu=0) yet exact
!
      if (.not. limplicit_viscosity .and. lvisc_simplified) then
        p%fvisc=p%fvisc+nu*p%del2u
        if (lpencil(i_visc_heat)) then
          if (lboussinesq) then
            p%visc_heat=p%visc_heat+2.*nu*p%sij2
          else                                  ! Heating term not implemented
            if (headtt) then
              call warning('calc_pencils_viscosity', 'viscous heating term '// &
                'is not implemented for lvisc_simplified')
            endif
          endif
        endif
        if (lfirst .and. ldt) p%diffus_total=p%diffus_total+nu
      endif
!
! Non-newtonian viscosity
!
      if (lvisc_nu_non_newtonian) then
!            print*,'viscous force: div(nu*Sij); with'
!            print*,'nu a function of square of strain rate'
! fvisc_i = \partial_j \nu(S_{kl}S_{kl}) S_{ij}
!         = \nu(S_{kl}S_{kl}) \partial_j S_{ij} + S_{ij}\partial_j \nu(S_{kl}S_{kl})
!         = \nu(S_{kl}S_{kl}) \partial_j (uij+uji) + S_{ij}\nu^{\prime} \partial_j [S_{kl}S_{kl}]
!         = \nu(S_{kl}S_{kl}) \partial_jj (ui) + S_{ij}\nu^{\prime} \partial_j [(ukl+ulk)(ukl+ulk)]
!         = \nu(.) \del2 (ui) + S_{ij}\nu^{\prime} \partial_j [ukl ukl +ulk ukl + ukl ulk + ulk ulk]
!         = \nu(.) \del2 (ui) + S_{ij}\nu^{\prime} [(uklj ukl + ukl uklj) + (ulkj ukl + ulk uklj)
!                                                  +(uklj ulk + ukl ulkj) + (ulkj ulk + ulk ulkj)]
!         = \nu(.) \del2 (ui) + S_{ij}\nu^{\prime} 2[uklj S_{lk}+  ulkj S_{kl}]
!
        deljskl2=0.
        gradnu_effective=0.
        do jj=1,3
          do kk=1,3; do ll=1,3
            deljskl2 (:,jj) = deljskl2 (:,jj) + &
              p%uijk(:,kk,ll,jj)*p%sij(:,ll,kk) + p%uijk(:,ll,kk,jj)*p%sij(:,kk,ll)
          enddo;enddo
        enddo
        do ii=1,3
          fvisc_nnewton2(:,ii) = sum(p%sij(:,ii,:)*deljskl2,2)
        enddo
        call getnu_non_newtonian(p%sij2,p%nu,gradnu_effective)
        do ii=1,3
          p%fvisc(:,ii)=p%nu(:)*p%del2u(:,ii) + gradnu_effective(:)*fvisc_nnewton2(:,ii)
        enddo
!
        if (lfirst .and. ldt) p%diffus_total=p%diffus_total+nu
      endif
!
!  viscous force: mu/rho*(del2u+graddivu/3)
!  -- the correct expression for rho*nu=const
!
      if (lvisc_rho_nu_const) then
        murho1=nu*p%rho1  !(=mu/rho)
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i) + &
              murho1*(p%del2u(:,i)+1.0/3.0*p%graddivu(:,i))
        enddo
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*murho1*p%sij2
        if (lfirst .and. ldt) p%diffus_total=p%diffus_total+murho1
      endif
!
!  viscous force: zeta/rho*graddivu (constant dynamic bulk viscosity)
!
      if (lvisc_rho_nu_const_bulk) then
        zetarho1=zeta*p%rho1  !(=zeta/rho)
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i)+zetarho1*p%graddivu(:,i)
        enddo
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+zetarho1*p%divu**2
        if (lfirst .and. ldt) p%diffus_total=p%diffus_total+zetarho1
      endif
!
!  viscous force: mu/sqrt(rho)*(del2u+graddivu/3)
!  [Implemented by Fred Gent in r13626 of 13 April 2010]
!  [AB: this may not correspond to a symmetric stress tensor]
!  PJK: 21-09-13 modified to include the missing term
!  viscous force: mu/sqrt(rho)*(del2u+graddivu/3+S.glnrho)
!
      if (lvisc_sqrtrho_nu_const) then
        murho1=nu*sqrt(p%rho1)
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i) + &
              murho1*(p%del2u(:,i)+1.0/3.0*p%graddivu(:,i) &
              + p%sglnrho(:,i))
        enddo
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*murho1*p%sij2
        if (lfirst .and. ldt) p%diffus_total=p%diffus_total+murho1
      endif
!
!  viscous force: nu*sqrt(TT)/rho*(del2u+graddivu/3+S.glnTT)
!
      if (lvisc_mu_therm) then
        muTT=nu*p%rho1*exp(0.5*p%lnTT)
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i) + &
              muTT*(p%del2u(:,i)+1.0/3.0*p%graddivu(:,i)+p%sglnTT(:,i))
        enddo
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*muTT*p%sij2
        if (lfirst .and. ldt) p%diffus_total=p%diffus_total+muTT
      endif
!
!  viscous force: nu*TT^2.5/rho*(del2u+graddivu/3+5S.glnTT)
!
      if (lvisc_spitzer) then
        muTT=nu_spitzer*p%rho1*exp(2.5*p%lnTT)
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i) + &
              muTT*(p%del2u(:,i)+1.0/3.0*p%graddivu(:,i)+5.*p%sglnTT(:,i))
        enddo
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*muTT*p%sij2
        if (lfirst .and. ldt) p%diffus_total=p%diffus_total+muTT
      endif
!
!  viscous force: nu*sqrt(TT)*(del2u+graddivu/3+2S.glnrho)
!  -- for numerical stability viscous force propto soundspeed in interstellar run
!
      if (lvisc_nu_therm) then
        muTT=nu*exp(0.5*p%lnTT)
        if (ldensity) then
          do i=1,3
            p%fvisc(:,i) = p%fvisc(:,i) + 2*muTT*p%sglnrho(:,i)&
                +muTT*(p%del2u(:,i) + 1./3.*p%graddivu(:,i))
          enddo
          ! Tobi: This is not quite the full story in the presence of linear
          ! shear. In this case the rate-of-strain tensor S has xy and yx
          ! components
          !   S_xy = S_yx = 1/2 (du_x/dy + du_y/dx + Sshear)
          ! so that one must add
          !if (.false.) then
          !  p%fvisc(:,1) = p%fvisc(:,1) + Sshear*muTT*p%glnrho(:,2)
          !  p%fvisc(:,2) = p%fvisc(:,2) + Sshear*muTT*p%glnrho(:,1)
          !endif
        else
          if (lmeanfield_nu) then
            if (meanfield_nuB/=0.) then
              call multsv_mn(muTT,p%del2u+1./3.*p%graddivu,tmp)
              call multsv_mn_add(1./sqrt(1.+p%b2/meanfield_nuB**2),&
                  p%fvisc+tmp,p%fvisc)
            endif
          else
            do i=1,3
              p%fvisc(:,i)=p%fvisc(:,i)+muTT*(p%del2u(:,i)&
                  +1.0/3.0*p%graddivu(:,i))
            enddo
          endif
        endif
!
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*muTT*p%sij2
        if (lfirst .and. ldt) p%diffus_total=p%diffus_total+muTT
      endif
!
!  viscous force: nu*(del2u+graddivu/3+2S.glnrho)
!  -- the correct expression for nu=const
!
      if (lvisc_nu_const) then
        if (ldensity) then
          fac=nu
          if (damp_sound/=0.) fac = fac+damp_sound*abs(p%divu)
          do j=1,3
            p%fvisc(:,j) = p%fvisc(:,j) + fac*(2*p%sglnrho(:,j)+p%del2u(:,j) + 1./3.*p%graddivu(:,j))
          enddo
          ! Tobi: This is not quite the full story in the presence of linear
          ! shear. In this case the rate-of-strain tensor S has xy and yx
          ! components
          !   S_xy = S_yx = 1/2 (du_x/dy + du_y/dx + Sshear)
          ! so that one must add
          !if (.false.) then
          !  p%fvisc(:,1) = p%fvisc(:,1) + Sshear*nu*p%glnrho(:,2)
          !  p%fvisc(:,2) = p%fvisc(:,2) + Sshear*nu*p%glnrho(:,1)
          !endif
        else
          if (lmeanfield_nu) then
            if (meanfield_nuB/=0.) then
              call multsv_mn_add(1./sqrt(1.+p%b2/meanfield_nuB**2), &
                  p%fvisc+nu*(p%del2u+1./3.*p%graddivu),p%fvisc)
            endif
          else
            p%fvisc=p%fvisc+nu*(p%del2u+1.0/3.0*p%graddivu)
          endif
        endif
!
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*nu*p%sij2
        if (lfirst .and. ldt) p%diffus_total=p%diffus_total+nu
      endif
!
!  viscous force: nu(t)*(del2u+graddivu/3+2S.glnrho)
!  -- the correct expression for nu=const
!
      if (lvisc_nu_tdep) then
        nu_tdep=nu*(t-nu_tdep_t0)**nu_tdep_exponent         ! out of nm loop
        p%fvisc=p%fvisc+2*nu_tdep*p%sglnrho+nu_tdep*(p%del2u+1./3.*p%graddivu)
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*nu_tdep*p%sij2
        if (lfirst .and. ldt) p%diffus_total=p%diffus_total+nu_tdep
      endif
!
!  Viscous force: nu*(del2u+graddivu/3+2S.glnrho)+2S.gradnu
!
      if (lvisc_mixture) then
        call multmv(p%sij,p%gradnu,sgradnu)
        do i=1,3
          p%fvisc(:,i)=p%nu*(p%del2u(:,i)+1./3.*p%graddivu(:,i)) + 2*sgradnu(:,i)
        enddo
        if (ldensity) then
          do i=1,3
            p%fvisc(:,i)=p%fvisc(:,i) + 2*p%nu*p%sglnrho(:,i)
          enddo
        endif
!
!  Viscous heating and time step.
!
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*p%nu*p%sij2
        if (lfirst .and. ldt) p%diffus_total=p%diffus_total+p%nu
      endif
!
!  Viscous force: nu(x)*(del2u+graddivu/3+2S.glnrho)+2S.gnu
!  -- here the nu viscosity depends on x; nu_jump=nu2/nu1
!        pnu = nu + (nu*(nu_jump-1.))*step(abs(p%x_mn),xnu,widthnu)
!
      if (lvisc_nu_profx.or.lvisc_nu_profr) then
        if (lvisc_nu_profx) tmp3=p%x_mn
        if (lvisc_nu_profr) then
          if (lspherical_coords.or.lsphere_in_a_box) then
            tmp3=p%r_mn
          else
            tmp3=p%rcyl_mn
          endif
        endif
        if (lvisc_nu_profx.and.lvisc_nu_profr) then
          print*,'You are using both radial and horizontal '
          print*,'profiles for a viscosity jump. Are you sure '
          print*,'this is reasonable? Better stop and check.'
          call fatal_error("calc_pencils_viscosity","")
        endif
        pnu = nu + (nu*(nu_jump-1.))*(step(tmp3,xnu ,widthnu) - step(tmp3,xnu2,widthnu))
        tmp4 = (nu*(nu_jump-1.))*(der_step(tmp3,xnu ,widthnu)-der_step(tmp3,xnu2,widthnu))
!
!  Initialize gradnu before calculating it, otherwise gfortran complains.
!
        gradnu=0.
        call get_gradnu(tmp4,lvisc_nu_profx,&
                             lvisc_nu_profr,p,gradnu)
!  A routine for write_xprof should be written here
!       call write_xprof('visc',pnu)
        !gradnu(:,1) = (nu*(nu_jump-1.))*der_step(tmp3,xnu,widthnu)
        !gradnu(:,2) = 0.
        !gradnu(:,3) = 0.
        call multmv(p%sij,gradnu,sgradnu)
        call multsv(pnu,2*p%sglnrho+p%del2u+1./3.*p%graddivu,tmp)
        !tobi: The following only works with operator overloading for pencils
        !      (see sub.f90). Commented out because it seems to be slower.
        !p%fvisc=p%fvisc+2*pnu*p%sglnrho+pnu*(p%del2u+1./3.*p%graddivu) &
        !        +2*sgradnu
        p%fvisc=p%fvisc+tmp+2*sgradnu
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*pnu*p%sij2
        if (lfirst .and. ldt) p%diffus_total=p%diffus_total+pnu
      endif
!
!  Radial viscosity profile from power law.
!  Viscous force: nu(x)*(del2u+graddivu/3+2S.glnrho)+2S.gnu
!  -- here the nu viscosity depends on r; nu=nu_0*(r/r0)^(-pnlaw)
!
      if (lvisc_nu_profr_powerlaw) then
        if (.not.luse_nu_rmn_prof) then
          pnu = nu*(p%rcyl_mn/xnu)**(-pnlaw)
!
!  viscosity gradient
!
          if (lcylindrical_coords) then
            gradnu(:,1) = -pnlaw*nu*(p%rcyl_mn/xnu)**(-pnlaw-1)*1/xnu
            gradnu(:,2) = 0.
            gradnu(:,3) = 0.
          elseif (lspherical_coords) then
            gradnu(:,1) = -pnlaw*nu*p%rcyl_mn**(-pnlaw-1)*sinth(m)
            gradnu(:,2) = -pnlaw*nu*p%rcyl_mn**(-pnlaw-1)*costh(m)
            gradnu(:,3) = 0.
          else
            print*,'power-law viscosity only implemented '
            print*,'for spherical and cylindrical coordinates'
            call fatal_error("calc_pencils_viscosity","")
          endif
        else
          pnu = nu*(p%r_mn/xnu)**(-pnlaw)
!
!  viscosity gradient
!
          if (lspherical_coords) then
            gradnu(:,1) = -pnlaw*nu*(p%r_mn/xnu)**(-pnlaw-1)*1/xnu
            gradnu(:,2) = 0.
            gradnu(:,3) = 0.
          else
            print*,'power-law viscosity with luse_nu_rmn_prof=T only '
            print*,'implemented for spherical coordinates'
            call fatal_error("calc_pencils_viscosity","")
          endif
        endif
!
        call multmv(p%sij,gradnu,sgradnu)
        call multsv(pnu,2*p%sglnrho+p%del2u+1./3.*p%graddivu,tmp)
        p%fvisc=p%fvisc+tmp+2*sgradnu
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*pnu*p%sij2
        if (lfirst .and. ldt) p%diffus_total=p%diffus_total+pnu
      endif
!
!  Radial viscosity profile with two steps; very similar to nu_profr
!  Viscous force: nu(x)*(del2u+graddivu/3+2S.glnrho)+2S.gnu
!  -- here the nu viscosity depends on r; nu_jump=nu2/nu1
!
      if (lvisc_nu_profr_twosteps) then
        if (lspherical_coords.or.lsphere_in_a_box) then
          tmp3=p%r_mn
        else
          tmp3=p%rcyl_mn
        endif
        if (lvisc_nu_profx.and.lvisc_nu_profr_twosteps) then
          print*,'You are using both radial and horizontal '
          print*,'profiles for a viscosity jump. Are you sure '
          print*,'this is reasonable? Better stop and check.'
          call fatal_error("","")
        endif
        prof2    = step(tmp3,xnu2,widthnu2)
        prof     = step(tmp3,xnu,widthnu)-prof2
        derprof2 = der_step(tmp3,xnu2,widthnu2)
        derprof  = der_step(tmp3,xnu,widthnu)-derprof2
!
        pnu  = nu + (nu*(nu_jump-1.))*prof+(nu*(nu_jump2-1.))*prof2
        tmp4 = nu + (nu*(nu_jump-1.))*derprof+(nu*(nu_jump2-1.))*derprof2
!
!  Initialize gradnu before calculating it, otherwise gfortran complains.
!
        gradnu=0.
        call get_gradnu(tmp4,lvisc_nu_profx,lvisc_nu_profr_twosteps,p,gradnu)
        call multmv(p%sij,gradnu,sgradnu)
        call multsv(pnu,2*p%sglnrho+p%del2u+1./3.*p%graddivu,tmp)
        p%fvisc=p%fvisc+tmp+2*sgradnu
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*pnu*p%sij2
        if (lfirst .and. ldt) p%diffus_total=p%diffus_total+pnu
      endif
!
!  turbulent viscosity profile from magnetic
!  viscous force: nu(z)*(del2u+graddivu/3+2S.glnrho)+2S.gnu
!  The following can in future be generalized to read
!
      if (lvisc_nut_from_magnetic) then
        pnu=PrM_turb*etat_x*etat_y(m)*etat_z(n)
        gradnu(:,1)=PrM_turb*detat_x*etat_y(m)*etat_z(n)
        gradnu(:,2)=PrM_turb*etat_x*detat_y(m)*etat_z(n)
        gradnu(:,3)=PrM_turb*etat_x*etat_y(m)*detat_z(n)
        call multmv(p%sij,gradnu,sgradnu)
        call multsv(pnu,2*p%sglnrho+p%del2u+1./3.*p%graddivu,tmp)
        p%fvisc=p%fvisc+tmp+2*sgradnu
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*pnu*p%sij2
        if (lfirst .and. ldt) p%diffus_total=p%diffus_total+pnu
      endif
!
!  viscous force: nu(z)*(del2u+graddivu/3+2S.glnrho)+2S.gnu
!  -- here the nu viscosity depends on z; nu_jump=nu2/nu1
!
      if (lvisc_nu_prof) then
        pnu = nu + (nu*(nu_jump-1.))*step(p%z_mn,znu,-widthnu)
!
!  This indicates that the profile is meant to be applied in Cartesian
!  coordinates only. Why then z taken from pencil?
!
        gradnu(:,1) = 0.
        gradnu(:,2) = 0.
        gradnu(:,3) = (nu*(nu_jump-1.))*der_step(p%z_mn,znu,-widthnu)
        call multmv(p%sij,gradnu,sgradnu)
        call multsv(pnu,2*p%sglnrho+p%del2u+1./3.*p%graddivu,tmp)
        !tobi: The following only works with operator overloading for pencils
        !      (see sub.f90). Commented out because it seems to be slower.
        !p%fvisc=p%fvisc+2*pnu*p%sglnrho+pnu*(p%del2u+1./3.*p%graddivu) &
        !        +2*sgradnu
        p%fvisc=p%fvisc+tmp+2*sgradnu
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*pnu*p%sij2
        if (lfirst .and. ldt) p%diffus_total=p%diffus_total+pnu
      endif
!
!  viscous force: nu_shock
!
      if (lvisc_nu_shock) then
        if (ldensity) then
          !tobi: The following only works with operator overloading for pencils
          !      (see sub.f90). Commented out because it seems to be slower.
          !tmp=nu_shock*(p%shock*(p%divu*p%glnrho+p%graddivu)+p%divu*p%gshock)
          call multsv(p%divu,p%glnrho,tmp2)
          tmp=tmp2 + p%graddivu
          call multsv(nu_shock*p%shock,tmp,tmp2)
          call multsv_add(tmp2,nu_shock*p%divu,p%gshock,tmp)
          p%fvisc=p%fvisc+tmp
          if (lfirst .and. ldt) p%diffus_total=p%diffus_total+(nu_shock*p%shock)
          if (lpencil(i_visc_heat)) &
              p%visc_heat=p%visc_heat+nu_shock*p%shock*p%divu**2
        endif
      endif
!
!  viscous force: nu_shock with vertical profile
!
      if (lvisc_nu_shock_profz) then
        if (ldensity) then
          pnu_shock = nu_shock + nu_shock_jump1*step(p%z_mn,znu_shock,-widthnu_shock)
!
!  This indicates that the profile is meant to be applied in Cartesian
!  coordinates only. Why then z taken from pencil?
!
          gradnu_shock(:,1) = 0.
          gradnu_shock(:,2) = 0.
          gradnu_shock(:,3) = nu_shock_jump1*der_step(p%z_mn,znu_shock,-widthnu_shock)
!
          call multsv(p%divu,p%glnrho,tmp2)
          tmp=tmp2 + p%graddivu
          call multsv(pnu_shock*p%shock,tmp,tmp2)
          call multsv_add(tmp2,pnu_shock*p%divu,p%gshock,tmp)
          call multsv_mn_add(p%shock*p%divu,gradnu_shock,tmp)
          p%fvisc=p%fvisc+tmp
          if (lfirst .and. ldt) p%diffus_total=p%diffus_total+(pnu_shock*p%shock)
          if (lpencil(i_visc_heat)) &
              p%visc_heat=p%visc_heat+pnu_shock*p%shock*p%divu**2
        endif
      endif

!
!  viscous force: nu_shock with radial profile
!
      if (lvisc_nu_shock_profr) then
        if (ldensity) then
          if (lspherical_coords.or.lsphere_in_a_box) then
            tmp3=p%r_mn
          else
            tmp3=p%rcyl_mn
          endif
          pnu_shock = nu_shock + nu_shock_jump1*step(tmp3,xnu_shock,widthnu_shock)
!
!  This indicates that the profile is meant to be applied in spherical and
!  cylindrical coordinates only, referring to r and pomega, respectively.
!  Hence not correct for 'sphere in a box'!
!  Why then r taken from pencil?
!
          gradnu_shock(:,1) = nu_shock_jump1*der_step(tmp3,xnu_shock,widthnu_shock)
          gradnu_shock(:,2) = 0.
          gradnu_shock(:,3) = 0.
!
          call multsv(p%divu,p%glnrho,tmp2)
          tmp=tmp2 + p%graddivu
          call multsv(pnu_shock*p%shock,tmp,tmp2)
          call multsv_add(tmp2,pnu_shock*p%divu,p%gshock,tmp)
          call multsv_mn_add(p%shock*p%divu,gradnu_shock,tmp)
          p%fvisc=p%fvisc+tmp
          if (lfirst .and. ldt) p%diffus_total=p%diffus_total+(pnu_shock*p%shock)
          if (lpencil(i_visc_heat)) &
              p%visc_heat=p%visc_heat+pnu_shock*p%shock*p%divu**2
        endif
      endif
!
!  viscous force: div(nu_shock * grad(uu_i))
!
      shksmp: if (lvisc_shock_simple) then
        do i = 1, 3
          call dot(p%gshock, p%uij(:,i,:), tmp3)
          tmp(:,i) = tmp3 + p%shock * p%del2u(:,i)
        enddo
        p%fvisc = p%fvisc + nu_shock * tmp
        if (lfirst .and. ldt) p%diffus_total = p%diffus_total + nu_shock * p%shock
        if (lpencil(i_visc_heat) .and. headtt) call warning('calc_pencils_viscosity', 'shock heating does not exist. ')
      endif shksmp
!
!  viscous force: nu_hyper2*de46v (not momentum-conserving)
!
      if (lvisc_hyper2_simplified) then
        p%fvisc=p%fvisc-nu_hyper2*p%del4u
        if (lpencil(i_visc_heat)) then  ! Heating term not implemented
          if (headtt) then
            call warning('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_hyper2_simplified')
          endif
        endif
        if (lfirst .and. ldt) p%diffus_total2=p%diffus_total2+nu_hyper2
      endif
!
!  viscous force: nu_hyper3*del6v (not momentum-conserving)
!
      if (lvisc_hyper3_simplified) then
        p%fvisc=p%fvisc+nu_hyper3*p%del6u
        if (lpencil(i_visc_heat)) then  ! Heating term not implemented
          if (headtt) then
            call warning('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_hyper3_simplified')
          endif
        endif
        if (lfirst .and. ldt) p%diffus_total3=p%diffus_total3+nu_hyper3
      endif
!
! General way of coding an anisotropic hyperviscosity.
!
      if (lvisc_hyper3_polar) then
        do j=1,3
          ju=j+iuu-1
          do i=1,3
            call der6(f,ju,tmp3,i,IGNOREDX=.true.)
            p%fvisc(:,j) = p%fvisc(:,j) + &
                nu_hyper3*pi4_1*tmp3*dline_1(:,i)**2
          enddo
          if (lpencil(i_visc_heat)) then
            if (headtt) then
              call warning('calc_pencils_viscosity', 'viscous heating term '//&
                   'is not implemented for lvisc_hyper3_polar')
            endif
          endif
          if (lfirst .and. ldt) &
               p%diffus_total3=p%diffus_total3+nu_hyper3*pi4_1/dxyz_4
        enddo
      endif
!
! Following Axel's hyper3_mesh for density
!
      if (lvisc_hyper3_mesh) then
        do j=1,3
          ju=j+iuu-1
          do i=1,3
            call der6(f,ju,tmp3,i,IGNOREDX=.true.)
            if (ldynamical_diffusion) then
              p%fvisc(:,j) = p%fvisc(:,j) + nu_hyper3_mesh * tmp3 * dline_1(:,i)
            else
              p%fvisc(:,j) = p%fvisc(:,j) + &
                  nu_hyper3_mesh*pi5_1/60.*tmp3*dline_1(:,i)
            endif
          enddo
          if (lpencil(i_visc_heat)) then
            if (headtt) then
              call warning('calc_pencils_viscosity', 'viscous heating term '//&
                   'is not implemented for lvisc_hyper3_mesh')
            endif
          endif
          if (lfirst .and. ldt) then
            if (ldynamical_diffusion) then
              p%diffus_total3 = p%diffus_total3 + nu_hyper3_mesh
              advec_hypermesh_uu = 0.0
            else
              advec_hypermesh_uu=nu_hyper3_mesh*pi5_1*sqrt(dxyz_2)
            endif
          endif
        enddo
      endif
!
!  viscous force: mu/rho*del6u
!
      if (lvisc_hyper3_rho_nu_const) then
        murho1=nu_hyper3*p%rho1  ! (=mu_hyper3/rho)
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i)+murho1*p%del6u(:,i)
        enddo
        if (lpencil(i_visc_heat)) then  ! Heating term not implemented
          if (headtt) then
            call warning('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_hyper3_rho_nu_const')
          endif
        endif
        if (lfirst .and. ldt) p%diffus_total3=p%diffus_total3+murho1
      endif
!
!  For tau_ij=d^5u_i/dx_j^5 + d^5u_j/dx_i^5
!  Viscous force: du/dt = mu/rho*{del6(u) + grad5[div(u)]}
!
      if (lvisc_hyper3_rho_nu_const_symm) then
        murho1=nu_hyper3*p%rho1  ! (=mu_hyper3/rho)
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i)+murho1*(p%del6u(:,i) + p%grad5divu(:,i))
        enddo
        if (lpencil(i_visc_heat)) then
          do i=1,3; do j=1,3
!
!  Dissipation is *not* positive definite.
!
            p%visc_heat=p%visc_heat + &
                0.5*nu_hyper3*(p%uij5(:,i,j)+p%uij5(:,j,i))*p%uij(:,i,j)
          enddo; enddo
        endif
        if (lfirst .and. ldt) p%diffus_total3=p%diffus_total3+murho1
      endif
!
!  Viscous force:
!      du/dt = mu/rho*{del2(del2(del2(u))) + 1/3*del2(del2(grad(div(u))
!  This viscosity type requires HYPERVISC_STRICT =  hypervisc_strict_fft in
!  your Makefile.local.
!
      if (lvisc_hyper3_mu_const_strict) then
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i)+nu_hyper3*p%rho1*f(l1:l2,m,n,ihypvis-1+i)
        enddo
        if (lpencil(i_visc_heat)) then  ! Should be eps=2*mu*{del2[del2(S)]}^2
          if (headtt) then              ! (see Haugen & Brandenburg 2004 eq. 7)
            call warning('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_hyper3_mu_const_strict')
          endif
        endif
        if (lfirst .and. ldt) p%diffus_total3=p%diffus_total3+nu_hyper3
      endif
!
!  Viscous force:
!      du/dt = 1/rho*div[2*rho*nu_3*S^(3)]
!  This viscosity type requires HYPERVISC_STRICT =  hypervisc_strict_fft in
!  your Makefile.local.
!
      if (lvisc_hyper3_nu_const_strict) then
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i)+nu_hyper3*f(l1:l2,m,n,ihypvis-1+i)
        enddo
        if (lpencil(i_visc_heat)) then
          if (headtt) then
            call warning('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_hyper3_nu_const_strict')
          endif
        endif
        if (lfirst .and. ldt) p%diffus_total3=p%diffus_total3+nu_hyper3
      endif
!
!  viscous force: f_i = mu_i/rho*del6u
!  Used for non-cubic cells
!
      if (lvisc_hyper3_rho_nu_const_aniso) then
        call del6fjv(f,nu_aniso_hyper3,iuu,tmp)
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i)+tmp(:,i)*p%rho1
        enddo
!
        if (lpencil(i_visc_heat)) then  ! Heating term not implemented
          if (headtt) then
            call warning('calc_pencils_viscosity', 'viscous heating term '// &
                 'is not implemented for lvisc_hyper3_rho_nu_const_aniso')
          endif
        endif
!
        if (lfirst .and. ldt) p%diffus_total3=p%diffus_total3+&
             (nu_aniso_hyper3(1)*dx_1(l1:l2)**6 + &     !!! out
              nu_aniso_hyper3(2)*dy_1(  m  )**6 + &
              nu_aniso_hyper3(3)*dz_1(  n  )**6)/ dxyz_6
!
      endif
!
!  viscous force: f_i = (nu_j.del6)u_i + nu_j.uij5.glnrho
!  Used for non-cubic cells
!
      if (lvisc_hyper3_nu_const_aniso) then
        call del6fjv(f,nu_aniso_hyper3,iuu,tmp)
!
        do i=1,3
          tmp3=0.
          do j=1,3
            tmp3=tmp3+p%uij(:,i,j)*p%glnrho(:,j)*nu_aniso_hyper3(j)
          enddo
!
          p%fvisc(:,i)=p%fvisc(:,i)+tmp(:,i)+tmp3
        enddo
!
        if (lpencil(i_visc_heat)) then  ! Heating term not implemented
          if (headtt) then
            call warning('calc_pencils_viscosity', 'viscous heating term '// &
                 'is not implemented for lvisc_hyper3_nu_const_aniso')
          endif
        endif
!
! diffusion time: it will be multiplied by dxyz_2 again further down
!
         if (lfirst .and. ldt) p%diffus_total3=p%diffus_total3+&
                 (nu_aniso_hyper3(1)*dx_1(l1:l2)**6 + &        !!! out
                  nu_aniso_hyper3(2)*dy_1(  m  )**6 + &
                  nu_aniso_hyper3(3)*dz_1(  n  )**6)/ dxyz_6
!
      endif
!
!  viscous force: mu/rho*d6uj/dx6
!
      if (lvisc_hyper3_rho_nu_const_bulk) then
        murho1=nu_hyper3*p%rho1  ! (=mu_hyper3/rho)
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i)+murho1*p%del6u_bulk(:,i)
        enddo
        if (lpencil(i_visc_heat)) then  ! Heating term not implemented
          if (headtt) then
            call warning('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_hyper3_rho_nu_const_bulk')
          endif
        endif
        if (lfirst .and. ldt) p%diffus_total3=p%diffus_total3+murho1
      endif
!
!  viscous force: nu_hyper3*(del6u+S.glnrho), where S_ij=d^5 u_i/dx_j^5
!
      if (lvisc_hyper3_nu_const) then
        p%fvisc=p%fvisc+nu_hyper3*(p%del6u+p%uij5glnrho)
        if (lpencil(i_visc_heat)) then  ! Heating term not implemented
          if (headtt) then
            call warning('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_hyper3_nu_const')
          endif
        endif
        if (lfirst .and. ldt) p%diffus_total3=p%diffus_total3+nu_hyper3
      endif
!
!  viscous force: Handle damping at the core of SNRs
!
      if (lvisc_smag_simplified) then
!
!  viscous force: nu_smag*(del2u+graddivu/3+2S.glnrho)
!  ??? where nu_smag=(C_smag*dxmax)**2*sqrt(2*SS)
!  where nu_smag=(C_smag*dxmax)**2*sqrt(SS^2)
!
        if (ldensity) then
!
!  Find nu_smag
!
          nu_smag=(C_smag*dxmax)**2.*sqrt(2*p%sij2)
!
!  with quenching term
!
          if (gamma_smag/=0.) nu_smag=nu_smag/sqrt(1.+gamma_smag*p%sij2)
!
! Calculate viscous force
!
          call multsv_mn(nu_smag,p%sglnrho,tmp2)
          call multsv_mn(nu_smag,p%del2u+1./3.*p%graddivu,tmp)
          p%fvisc=p%fvisc+2*tmp2+tmp
          if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*nu_smag*p%sij2
          if (lfirst .and. ldt) p%diffus_total=p%diffus_total+nu_smag
        else
          if (lfirstpoint) print*, 'calc_viscous_force: '// &
              "ldensity better be .true. for ivisc='smagorinsky'"
!AB: not yet programmed
        endif
      endif
!
!  viscous force: nu_smag*(del2u+graddivu/3+2S.glnrho)
!  where nu_smag=(C_smag*dxmax)**2*sqrt(S:J)
!
      if (lvisc_smag_cross_simplified) then
        if (ldensity) then
          nu_smag=(C_smag*dxmax)**2.*p%ss12
!
! Calculate viscous force
!
          call multsv_mn(nu_smag,p%sglnrho,tmp2)
          call multsv_mn(nu_smag,p%del2u+1./3.*p%graddivu,tmp)
          p%fvisc=p%fvisc+2*tmp2+tmp
          if (lpencil(i_visc_heat)) then  ! Heating term not implemented
            if (headtt) then
              call warning('calc_pencils_viscosity','viscous heating term '//&
                'is not implemented for lvisc_smag_cross_simplified')
            endif
          endif
          if (lfirst .and. ldt) p%diffus_total=p%diffus_total+nu_smag
        else
          if (lfirstpoint) print*, 'calc_viscous_force: '// &
              "ldensity better be .true. for ivisc='smagorinsky'"
        endif
      endif
!
!  Calculate viscouse force form a slope limited diffusion
!  following Rempel (2014).
!  Here calculating the divergence of the flux.
!
      if (lvisc_slope_limited) then
        p%fvisc(:,iuu:iuu+2)=p%fvisc(:,iuu:iuu+2)-f(l1:l2,m,n,iFF_div_uu:iFF_div_uu+2)
        if (lfirst .and. ldt) then
          tmp=0.
          where (p%del2u /=0.)
!            tmp=f(l1:l2,m,n,iFF_div_uu:iFF_div_uu+2)/p%del2u
            tmp=f(l1:l2,m,n,iFF_div_uu:iFF_div_uu+2)/p%uu*dx**2
          endwhere
          nu_sld=maxval(abs(tmp),2)
!          print*,maxval(nu_sld)
          p%diffus_total=p%diffus_total+nu_sld
        endif
      endif
!
!  Calculate Lambda effect
!
      if (llambda_effect) then
        call calc_lambda(p,lambda_phi)
        p%fvisc(:,iuz)=p%fvisc(:,iuz) + lambda_phi
      endif
!
!  Store viscous heating rate in auxiliary variable if requested.
!  Just neccessary immediately before writing snapshots, but how would we
!  know we are?
!
      if (lvisc_heat_as_aux) f(l1:l2,m,n,ivisc_heat) = p%visc_heat
!
!  Do diagnostics related to viscosity.
!
      if (ldiagnos) then
        if (idiag_nusmagm/=0)   call sum_mn_name(nu_smag,idiag_nusmagm)
        if (idiag_nusmagmin/=0) call max_mn_name(-nu_smag,idiag_nusmagmin,lneg=.true.)
        if (idiag_nusmagmax/=0) call max_mn_name(nu_smag,idiag_nusmagmax)
        if (idiag_num/=0) call sum_mn_name(p%nu,idiag_num)
        if (idiag_qfviscm/=0) then
          call dot(p%curlo,p%fvisc,qfvisc)
          call sum_mn_name(qfvisc,idiag_qfviscm)
        endif
      endif
!
    endsubroutine calc_pencils_viscosity
!***********************************************************************
    subroutine viscosity_after_boundary(f)

      use Sub, only: div
      use Sub, only: notanumber

      real, dimension (mx,my,mz,mfarray) :: f

      integer :: ll,mm,nn,j,iff
      real, dimension(mx-1) :: tmpx
      real, dimension(my-1) :: tmpy
      real, dimension(mz-1) :: tmpz
      real, dimension(nx)   :: tmp
!
!  Slope limited diffusion following Rempel (2014)
!  First calculating the flux in a subroutine below
!  using a slope limiting procedure then storing in the
!  auxilaries variables in the f array (done above).
!
      if (lvisc_slope_limited) then
!
!  to avoid taking the sqrt several times
!
        f(:,:,:,iFF_char_c)=sqrt(f(:,:,:,iFF_char_c))
        f(:,:,:,iFF_diff1:iFF_diff2)=0.
!
       do j=1,3

          iff=iFF_diff

          if (nxgrid>1) then
            do nn=n1,n2; do mm=m1,m2
              tmpx = f(2:,mm,nn,iuu+j-1)-f(:mx-1,mm,nn,iuu+j-1)
!if (lroot.and.j==1.and.lfirst.and.ldiagnos) print*, 'tmpx=',tmpx
              call calc_diffusive_flux(tmpx,f(2:mx-2,mm,nn,iFF_char_c),f(2:mx-2,mm,nn,iff))
!if (lroot.and.j==1.and.lfirst.and.ldiagnos) print*,'flux=',f(2:mx-2,mm,nn,iff)
!if (notanumber(f(2:mx-2,mm,nn,iff))) print*, 'DIFFX:j,mm,nn=', j,mm,nn
!            if(lfirst.and.ldiagnos.and.j==1) print*,f(473:533,mm,nn,iff)
            enddo; enddo
            iff=iff+1
          endif

          if (nygrid>1) then
            do nn=n1,n2; do ll=l1,l2
              tmpy = f(ll,2:,nn,iuu+j-1)-f(ll,:my-1,nn,iuu+j-1)
              call calc_diffusive_flux(tmpy,f(ll,2:my-2,nn,iFF_char_c),f(ll,2:my-2,nn,iff))
!if (notanumber(f(ll,2:my-1,nn,iFF_diff+1))) print*, 'DIFFY:j,ll,nn=', j,ll,nn
            enddo; enddo
            iff=iff+1
          endif

          if (nzgrid>1) then
            do mm=m1,m2; do ll=l1,l2
              tmpz = f(ll,mm,2:,iuu+j-1)-f(ll,mm,:mz-1,iuu+j-1)
!if (notanumber(tmpz)) print*, 'DIFFZ:j,ll,mm=', j,ll,mm
            call calc_diffusive_flux(tmpz,f(ll,mm,2:mz-2,iFF_char_c),f(ll,mm,2:mz-2,iff))
!if (notanumber(f(ll,mm,2:mz-1,iFF_diff+2))) print*, 'DIFFZ:j,ll,mm=', j,ll,mm
            enddo; enddo
          endif

          do n=n1,n2; do m=m1,m2
!             if(lfirst.and.ldiagnos.and.j==1) print*,f(473:533,mm,nn,iFF_diff)
            call div(f,iFF_diff,f(l1:l2,m,n,iFF_div_uu+j-1),iorder=4)
!if (lroot.and.lfirst.and.ldiagnos.and.j==1) print*,'DIVflux=',f(508:511,m,n,iFF_diff1)
!            call div(f,iFF_diff,tmp,iorder=4)
!if (lroot.and.j==1.and.lfirst.and.ldiagnos) print*, tmp
          enddo; enddo

        enddo
      endif

    endsubroutine viscosity_after_boundary
!***********************************************************************
    subroutine calc_diffusive_flux(diffs,c_char,flux)
!
!  23-sep-15/MR,joern,fred,petri: coded
!
      use Sub, only: slope_limiter, diff_flux

      real, dimension(:),intent(in) :: diffs,c_char
      real, dimension(:),intent(out):: flux

      real, dimension(size(diffs)-1) :: slope
      real, dimension(size(diffs)-2) :: phi
      integer :: len

      len=size(diffs)

      call slope_limiter(diffs(2:),diffs(:len-1),slope,islope_limiter)

      flux = diffs(2:len-1) - 0.5*(slope(2:) + slope(1:len-2))

      call diff_flux(h_slope_limited, diffs(2:len-1), flux, phi)
      flux = -0.5*c_char*phi*flux

    endsubroutine calc_diffusive_flux
!***********************************************************************
    subroutine getnu_non_newtonian(gdotsqr,nu_effective,gradnu_effective)
!
! Outputs non-newtonian viscosity for an input strain-rate squared
!
      use Sub, only : step,der_step
      real, dimension(nx) :: gdotsqr
      real,dimension(nx) :: nu_effective,gradnu_effective
!
      select case(nnewton_type)
        case('carreau')
          nu_effective= nu_infinity+ &
            (nu0-nu_infinity)/(1+(non_newton_lambda**2)*gdotsqr)**((1.-carreau_exponent)/2.)
          gradnu_effective= (non_newton_lambda**2)*0.5*(carreau_exponent-1.)* &
            (nu0-nu_infinity)/(1+(non_newton_lambda**2)*gdotsqr)**(1.5-carreau_exponent/2.)
        case('step')
          nu_effective=nu0+(nu_infinity-nu0)*step(gdotsqr,nnewton_tscale**2,nnewton_step_width)
          gradnu_effective=(nu_infinity-nu0)*der_step(gdotsqr,nnewton_tscale**2,nnewton_step_width)
        case default
          write(*,*) 'nnewton_type=',nnewton_type
          call fatal_error('viscosity,getnu_non_newtonian:','select nnewton_type')
      endselect
!
    endsubroutine getnu_non_newtonian
!***********************************************************************
    subroutine calc_viscosity(f)
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_viscosity
!***********************************************************************
    subroutine get_gradnu(tmp,ljumpx,ljumpr,p,gradnu)
!
!  Calculate grad nu for vicosity jump in different
!  coordinate systems
!
!  20-jun-08/wlad :: coded
!
      real, dimension(nx)  ,intent(in) :: tmp
      real, dimension(nx,3),intent(out) :: gradnu
      logical, intent (in) :: ljumpx,ljumpr
      type (pencil_case) :: p
!
      if (ljumpx.or.                             &
         (ljumpr.and.(lcylindrical_coords.or.    &
                      lspherical_coords       ))) then
        gradnu(:,1)=tmp ; gradnu(:,2)=0 ; gradnu(:,3)=0
      elseif (ljumpr.and.lcylinder_in_a_box) then
        gradnu(:,1)=tmp*x(l1:l2)*p%rcyl_mn1
        gradnu(:,2)=tmp*y(  m  )*p%rcyl_mn1
      elseif (ljumpr.and.lsphere_in_a_box) then
        print*,'get_gradnu: not yet implemented for '
        print*,'embedded spheres'
        call fatal_error("","")
      endif
!
    endsubroutine get_gradnu
!***********************************************************************
    subroutine calc_viscous_heat(df,p,Hmax)
!
!  Calculate viscous heating term for right hand side of entropy equation.
!
!  20-nov-02/tony: coded
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: Hmax
!
!  Add viscous heat (which has units of energy/mass) to the RHS
!  of the entropy (both with and without pretend_lnTT), or of
!  the temperature equation. Divide by cv if pretend_lnTT.
!
      if (lentropy) then
!
        if (pretend_lnTT) then
          df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + p%cv1*p%TT1*p%visc_heat
        else
          df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + p%TT1*p%visc_heat
        endif
      else if (ltemperature) then
        if (ltemperature_nolog) then
          df(l1:l2,m,n,iTT)   = df(l1:l2,m,n,iTT)   + p%cv1*p%visc_heat
        else
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + p%cv1*p%TT1*p%visc_heat
        endif
      else if (lthermal_energy) then
        if (lstratz) then
          df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) + p%rho * p%visc_heat / eth0z(n)
        else
          df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) + p%rho * p%visc_heat
        endif
      endif
!
!  Calculate maximum heating (for time step constraint), so it is
!  only done on the first of the 3 substeps.
!
      if (lfirst .and. ldt) Hmax=Hmax+p%visc_heat
!
    endsubroutine calc_viscous_heat
!***********************************************************************
    subroutine calc_viscous_force(df,p)
!
!  Calculate viscous force term for right hand side of equation of motion.
!
!  20-nov-02/tony: coded
!   9-jul-04/nils: added Smagorinsky viscosity
!
      use Diagnostics, only: sum_mn_name, max_mn_name, xysum_mn_name_z, &
          yzsum_mn_name_x, zsum_mn_name_xy, max_mn_name
      use Sub, only: cross, dot2
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: nu_smag,Reshock,fvisc2
      real, dimension (nx,3) :: nuD2uxb
      type (pencil_case) :: p
      integer :: i
!
      intent (in) :: p
      intent (inout) :: df
!
!  Add viscosity to equation of motion
!
      if (lmagfield_nu) then
        do i=1,3
          df(l1:l2,m,n,iux+i-1) = df(l1:l2,m,n,iux+i-1) + p%fvisc(:,i)/(1.+p%b2/meanfield_nuB**2)
        enddo
      else
        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + p%fvisc
      endif
!
!  Calculate max total diffusion coefficient for timestep calculation etc.
!
      if (lfirst .and. ldt) then
        diffus_nu =p%diffus_total *dxyz_2
        diffus_nu2=p%diffus_total2*dxyz_4
        if (ldynamical_diffusion .and. lvisc_hyper3_mesh) then
          diffus_nu3 = p%diffus_total3 * (abs(dline_1(:,1)) + abs(dline_1(:,2)) + abs(dline_1(:,3)))    !!! out
        else
          diffus_nu3=p%diffus_total3*dxyz_6
        endif
      endif
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_nu_tdep/=0)  call sum_mn_name(spread(nu_tdep,1,nx),idiag_nu_tdep)
        !!!if (idiag_fviscm/=0)   call sum_mn_name(p%fvisc,idiag_fviscm)   !What is intended here? p%fvisc is a vector!!
        if (idiag_fviscm/=0)   call sum_mn_name(p%fvisc(:,1),idiag_fviscm)
        if (idiag_fviscmin/=0) call max_mn_name(-p%fvisc,idiag_fviscmin,lneg=.true.)
        if (idiag_fviscmax/=0) call max_mn_name(p%fvisc,idiag_fviscmax)
        if (idiag_fviscrmsx/=0) then
           call dot2(p%fvisc,fvisc2)
           call sum_mn_name(xmask_vis*fvisc2,idiag_fviscrmsx,lsqrt=.true.)
        endif
        if (lvisc_smag_simplified) then
          if (ldensity) nu_smag=(C_smag*dxmax)**2.*sqrt(2*p%sij2)
        endif
        if (lvisc_smag_cross_simplified) then
          if (ldensity) nu_smag=(C_smag*dxmax)**2.*p%ss12
        endif
        if (idiag_dtnu/=0) &
            call max_mn_name(diffus_nu/cdtv,idiag_dtnu,l_dt=.true.)
        if (idiag_nu_LES /= 0) call sum_mn_name(nu_smag,idiag_nu_LES)
        if (idiag_meshRemax/=0) call max_mn_name(sqrt(p%u2(:))*dxmax_pencil/p%diffus_total,idiag_meshRemax)
        if (idiag_Reshock/=0) then
          Reshock(:) = 0.
          where (abs(p%shock) > tini)
            Reshock = dxmax_pencil*sqrt(p%u2)/(nu_shock*p%shock)
          endwhere
          call max_mn_name(Reshock,idiag_Reshock)
        endif
        if (idiag_visc_heatm/=0) call sum_mn_name(p%visc_heat,idiag_visc_heatm)
        if (idiag_Sij2m/=0) call sum_mn_name(p%sij2,idiag_Sij2m)
        if (idiag_epsK/=0) call sum_mn_name(p%visc_heat*p%rho,idiag_epsK)
!
!  Viscous heating for Smagorinsky viscosity.
!
        if (idiag_epsK_LES/=0) then
          if (lvisc_smag_simplified) then
            call sum_mn_name(2*nu_smag*p%rho*p%sij2,idiag_epsK_LES)
          else if (lvisc_smag_cross_simplified) then
            call sum_mn_name(2*nu_smag*p%rho*p%sij2,idiag_epsK_LES)
          endif
        endif
!
!  correlation of viscous term with b-field; relevant for MTA
!  (MTA: minimal tau approximation)
!
        if (lmagnetic) then
          if (idiag_nuD2uxbxm/=0.or.idiag_nuD2uxbym/=0.or.idiag_nuD2uxbzm/=0) then
            call cross(p%fvisc,p%bb,nuD2uxb)
            call sum_mn_name(nuD2uxb(:,1),idiag_nuD2uxbxm)
            call sum_mn_name(nuD2uxb(:,2),idiag_nuD2uxbym)
            call sum_mn_name(nuD2uxb(:,3),idiag_nuD2uxbzm)
          endif
        endif
      endif
!
!  1D-averages.
!
      if (l1davgfirst) then
        if (idiag_fviscmz/=0) &
            call xysum_mn_name_z(-2.*p%rho*nu*( &
            p%uu(:,1)*p%sij(:,1,3)+ &
            p%uu(:,2)*p%sij(:,2,3)+ &
            p%uu(:,3)*p%sij(:,3,3)),idiag_fviscmz)
        if (idiag_epsKmz/=0) &
            call xysum_mn_name_z(p%visc_heat*p%rho,idiag_epsKmz)
        if (idiag_fviscmx/=0) &
            call yzsum_mn_name_x(-2.*p%rho*nu*( &
            p%uu(:,1)*p%sij(:,1,1)+ &
            p%uu(:,2)*p%sij(:,2,1)+ &
            p%uu(:,3)*p%sij(:,3,1)),idiag_fviscmx)
        if (idiag_numx/=0) &
            call yzsum_mn_name_x(p%nu,idiag_numx)
      endif
!
!  2D-averages.
!
      if (l2davgfirst) then
        if (idiag_fviscmxy/=0) call zsum_mn_name_xy(-2.*p%rho*nu*( &
            p%uu(:,1)*p%sij(:,1,1)+p%uu(:,2)*p%sij(:,2,1)+ &
            p%uu(:,3)*p%sij(:,3,1)),idiag_fviscmxy)
        if (idiag_fviscmxy/=0.and.lvisc_sqrtrho_nu_const) &
            call zsum_mn_name_xy(-2.*sqrt(p%rho)*nu*( &
            p%uu(:,1)*p%sij(:,1,1)+p%uu(:,2)*p%sij(:,2,1)+ &
            p%uu(:,3)*p%sij(:,3,1)),idiag_fviscmxy)
        if (idiag_fviscmxy/=0.and.lvisc_rho_nu_const) &
            call zsum_mn_name_xy(-2.*nu*( &
            p%uu(:,1)*p%sij(:,1,1)+p%uu(:,2)*p%sij(:,2,1)+ &
            p%uu(:,3)*p%sij(:,3,1)),idiag_fviscmxy)
        if (idiag_fviscymxy/=0) call zsum_mn_name_xy(-2.*p%rho*nu*( &
            p%uu(:,1)*p%sij(:,1,2)+p%uu(:,2)*p%sij(:,2,2)+ &
            p%uu(:,3)*p%sij(:,3,2)),idiag_fviscymxy)
      endif
!
    endsubroutine calc_viscous_force
!***********************************************************************
    subroutine calc_visc_heat_ppd(df,p)
!
!  Calculates viscous dissipation term from D'Angelo et al. (2003)
!
!  03/08 steveb
!
      use Sub, only: get_radial_distance
      use Gravity, only: acceleration
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: diss
      real, dimension (nx) :: rr_sph,rr_cyl,g_r,OO2
!
! 'Dissipative' heating term Y=9/4 \Sigma \nu \Omega_K^2
!  Need to get correct angular velocity first
!
      call get_radial_distance(rr_sph,rr_cyl)
      call acceleration(g_r)
      ! should generalize this to non-cartesian; see hydro L1937
      OO2=max(-g_r/rr_cyl,0.)
!
! only works for 2-D r-phi disks
! this goes into entropy equation, which divides by \rho, so no density here
! this is really surface density, though
!
      if (nzgrid==1) then
        diss = 9./4.*nu*OO2
        df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + p%TT1*diss
!
      else
        call fatal_error("calc_visc_heat_ppd", &
            "dissipation only implemented for 2d-disk")
      endif
!
    endsubroutine calc_visc_heat_ppd
!***********************************************************************
    subroutine getnu(p,nu_input,nu_pencil,ivis)
!
!  Return the viscosity (by value)
!
!  14-aug-08/kapelrud: coded
!
      type (pencil_case), optional :: p
      real, optional, intent(out) :: nu_input
      real, dimension(nx), optional, intent(out) :: nu_pencil
      character (len=labellen), optional :: ivis
!
      if (present(nu_input)) nu_input=nu
      if (present(ivis))     ivis=ivisc(1)
      if (present(nu_pencil)) then
        if (.not. present(p)) then
          call fatal_error('getnu',&
              'p must be present in call to getnu when nu_pencil is present!')
        endif
        if (lvisc_simplified) then
          nu_pencil=nu
        elseif(lvisc_mixture) then
          nu_pencil=nu
        elseif (lvisc_rho_nu_const) then
          nu_pencil=nu*p%rho1
        elseif (lvisc_rho_nu_const_bulk) then
          nu_pencil=zeta*p%rho1
        elseif (lvisc_nu_therm) then
          nu_pencil=nu*sqrt(exp(p%lnTT))
        elseif (lvisc_mu_therm) then
          nu_pencil=nu*sqrt(exp(p%lnTT))*p%rho1
        elseif (lvisc_nu_const) then
          nu_pencil=nu
        else
          call fatal_error('getnu','No such ivisc implemented!')
        endif
      endif
!
    endsubroutine getnu
!***********************************************************************
    subroutine dynamical_viscosity(urms)
!
!  Dynamically set viscosity coefficient given fixed mesh Reynolds number.
!
!  27-jul-11/ccyang: coded
!
      real, intent(in) :: urms
!
!  Hyper-viscosity coefficient
!
      if (nu_hyper3/=0.0) nu_hyper3 = pi5_1 * urms * dxmax**5 / re_mesh
      if (nu_hyper3_mesh/=0.0) nu_hyper3_mesh = pi5_1 * urms / re_mesh / sqrt(real(dimensionality))
!
    endsubroutine dynamical_viscosity
!***********************************************************************
    subroutine split_update_viscosity(f)
!
!  Update the velocity by integrating the operator split viscous terms.
!
!  22-aug-13/ccyang: coded.
!
      use ImplicitDiffusion, only: integrate_diffusion
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      if (limplicit_viscosity) &
        call integrate_diffusion(get_viscosity_implicit, f, iux, iuz)
!
    endsubroutine split_update_viscosity
!***********************************************************************
    subroutine get_viscosity_implicit(ndc, diffus_coeff, iz)
!
!  Gets the diffusion coefficient along a given pencil for the implicit algorithm.
!
!  22-aug-13/ccyang: coded.
!
      integer, intent(in) :: ndc
      real, dimension(ndc), intent(out) :: diffus_coeff
      integer, intent(in), optional :: iz
!
      if (lvisc_simplified) then
        diffus_coeff = nu
      else
        diffus_coeff = 0.0
      endif
!
    endsubroutine get_viscosity_implicit
!***********************************************************************
    subroutine calc_lambda(p,div_lambda)
!
!  Calculates the lambda effect
!
!  20-apr-10/dhruba: coded
! If lKit_Olem is true the lambda coefficients depends on dsdr which must be
! incorporated below.
!
!
      use cdata, only: Omega
!
      real,dimension(nx) :: div_lambda,lomega,dlomega_dr,dlomega_dtheta, &
                            lver,lhor,dlver_dr,dlhor_dtheta
      type (pencil_case) :: p
!
      lomega=p%uu(:,3)/(sinth(m)*x(l1:l2))+Omega
!
      dlomega_dr=(x(l1:l2)*p%uij(:,3,1)-p%uu(:,3))/ &
          (sinth(m)*x(l1:l2)*x(l1:l2))
      dlomega_dtheta=(p%uij(:,3,2)*x(l1:l2)-p%uu(:,3)*cotth(m))/ &
          (sinth(m)*x(l1:l2)*x(l1:l2))
!
      lver = -(Lambda_V0*LV0_rprof(l1:l2)+Lambda_V1*sinth(m)*sinth(m) &
          *LV1_rprof(l1:l2) )
      lhor = -Lambda_H1*sinth(m)*sinth(m)*LH1_rprof(l1:l2)
!
      dlver_dr = -(Lambda_V0*der_LV0_rprof(l1:l2)+Lambda_V1 &
          *sinth(m)*sinth(m)*der_LV1_rprof(l1:l2))
      dlhor_dtheta = -Lambda_H1*LH1_rprof(l1:l2)*2.*costh(m)*sinth(m)/x(l1:l2)
!
      div_lambda = lver*(sinth(m)*lomega*p%glnrho(:,1) &
          +3.*sinth(m)*lomega/x(l1:l2) + sinth(m)*dlomega_dr) &
          +lomega*sinth(m)*dlver_dr + lhor*(costh(m)*lomega*p%glnrho(:,2) &
          -sinth(m)*lomega/x(l1:l2) + 2.*cotth(m)*costh(m)*lomega/x(l1:l2) &
          +costh(m)*dlomega_dtheta) + lomega*costh(m)*dlhor_dtheta
!
    endsubroutine calc_lambda
!***********************************************************************
endmodule Viscosity
