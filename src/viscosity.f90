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
! PENCILS PROVIDED visc_heat; nu; gradnu(3); sgnu(3)
!
!***************************************************************
module Viscosity
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'viscosity.h'
!
  integer, parameter :: nvisc_max=4
  character (len=labellen), dimension(nvisc_max) :: ivisc=''
  character (len=labellen) :: lambda_profile='uniform'
  real :: nu=0.0, nu_tdep=0.0, nu_tdep_exponent=0.0, nu_tdep_t0=0.0
  real :: nu_mol=0.0, nu_hyper2=0.0, nu_hyper3=0.0
  real :: nu_hyper3_mesh=5.0, nu_shock=0.0,nu_spitzer=0.0
  real :: nu_jump=1.0, xnu=1.0, xnu2=1.0, znu=1.0, widthnu=0.1, C_smag=0.0
  real :: pnlaw=0.0, Lambda_V0=0.,Lambda_V1=0.,Lambda_H1=0.
  real :: Lambda_V0t=0.,Lambda_V1t=0.,Lambda_V0b=0.,Lambda_V1b=0.
  real :: rzero_lambda=impossible,wlambda=0.,rmax_lambda=impossible
  real :: offamp_lambda=1.,r1_lambda=impossible,r2_lambda=impossible
  real :: lambda_jump=0.,roffset_lambda=0.
  real :: PrM_turb=0.0
  real :: meanfield_nuB=0.0
  real, dimension(:), pointer :: etat_x, detat_x
  real, dimension(:), pointer :: etat_y, detat_y
  real, dimension(:), pointer :: etat_z, detat_z
  real, dimension(3) :: nu_aniso_hyper3=0.0
  real, dimension(mx) :: LV0_rprof,LV1_rprof,LH1_rprof,der_LV0_rprof,der_LV1_rprof
  logical :: lvisc_first=.false.
  logical :: lvisc_simplified=.false.
  logical :: lvisc_rho_nu_const=.false.
  logical :: lvisc_sqrtrho_nu_const=.false.
  logical :: lvisc_nu_therm=.false.
  logical :: lvisc_mu_therm=.false.
  logical :: lvisc_nu_const=.false.
  logical :: lvisc_nu_tdep=.false.
  logical :: lvisc_nu_prof=.false.
  logical :: lvisc_nu_profx=.false.
  logical :: lvisc_nu_profr=.false.
  logical :: lvisc_nu_profr_powerlaw=.false.
  logical :: lvisc_nut_from_magnetic=.false.
  logical :: lvisc_nu_shock=.false.
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
  logical :: lmeanfield_nu=.false.
  logical :: lmagfield_nu=.false.
  logical :: llambda_effect=.false.
  logical, pointer:: lviscosity_heat
!
  namelist /viscosity_run_pars/ &
      nu, nu_tdep_exponent, nu_tdep_t0, &
      nu_hyper2, nu_hyper3, ivisc, nu_mol, C_smag, nu_shock, &
      nu_aniso_hyper3, lvisc_heat_as_aux,nu_jump,znu,xnu,xnu2,widthnu, &
      pnlaw,llambda_effect,Lambda_V0,Lambda_V1,Lambda_H1, nu_hyper3_mesh, &
      lambda_profile,rzero_lambda,wlambda,r1_lambda,r2_lambda,rmax_lambda,&
      offamp_lambda,lambda_jump,lmeanfield_nu,lmagfield_nu,meanfield_nuB, &
      PrM_turb, roffset_lambda, nu_spitzer
!
! other variables (needs to be consistent with reset list below)
  integer :: idiag_nu_tdep=0    ! DIAG_DOC: time-dependent viscosity
  integer :: idiag_fviscm=0     ! DIAG_DOC: Mean value of viscous acceleration
  integer :: idiag_fviscmin=0   ! DIAG_DOC: Min value of viscous acceleration
  integer :: idiag_fviscmax=0   ! DIAG_DOC: Max value of viscous acceleration
  integer :: idiag_nusmagm=0    ! DIAG_DOC: Mean value of Smagorinsky viscosity
  integer :: idiag_nusmagmin=0  ! DIAG_DOC: Min value of Smagorinsky viscosity
  integer :: idiag_nusmagmax=0  ! DIAG_DOC: Max value of Smagorinsky viscosity
  integer :: idiag_visc_heatm=0 ! DIAG_DOC: Mean value of viscous heating
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
!
! z averaged diagnostics given in zaver.in
!
  integer :: idiag_fviscmxy=0   ! ZAVG_DOC: $\left<2\nu\varrho u_i
                                ! ZAVG_DOC: \mathcal{S}_{ix} \right>_{z}$
                                ! ZAVG_DOC: ($x$-xomponent of viscous flux)
!
  contains
!***********************************************************************
    subroutine register_viscosity()
!
!  19-nov-02/tony: coded
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_viscosity
!***********************************************************************
    subroutine initialize_viscosity(lstarting)
!
!  20-nov-02/tony: coded
!
      use FArrayManager, only: farray_register_auxiliary
      use Mpicomm, only: stop_it
      use SharedVariables, only: put_shared_variable,get_shared_variable
!
      logical, intent(in) :: lstarting
!
      integer :: i, ierr
!
!  Default viscosity.
!
      if ( (nu/=0.0).and.(ivisc(1)=='') ) ivisc(1)='nu-const'
!
!  Some viscosity types need the rate-of-strain tensor and grad(lnrho)
!
      lvisc_simplified=.false.
      lvisc_rho_nu_const=.false.
      lvisc_sqrtrho_nu_const=.false.
      lvisc_nu_therm=.false.
      lvisc_mu_therm=.false.
      lvisc_nu_const=.false.
      lvisc_nu_tdep=.false.
      lvisc_nu_prof=.false.
      lvisc_nu_profx=.false.
      lvisc_nu_profr=.false.
      lvisc_nu_profr_powerlaw=.false.
      lvisc_nut_from_magnetic=.false.
      lvisc_nu_shock=.false.
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
        case ('simplified', '0')
          if (lroot) print*,'viscous force: nu*del2v'
          lvisc_simplified=.true.
        case ('rho-nu-const','rho_nu-const', '1')
          if (lroot) print*,'viscous force: mu/rho*(del2u+graddivu/3)'
          lvisc_rho_nu_const=.true.
        case ('sqrtrho-nu-const')
          if (lroot) print*,'viscous force: mu/sqrt(rho)*(del2u+graddivu/3)'
          if (nu/=0.) lpenc_requested(i_sij)=.true.
          lvisc_sqrtrho_nu_const=.true.
        case ('nu-therm')
          if (lroot) print*,'viscous force: mu*sqrt(TT)*(del2u+graddivu/3)'
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
        case ('snr-damp','snr_damp')
          if (lroot) print*,'viscous force: SNR damping'
          lvisc_snr_damp=.true.
        case ('nu-mixture')
          if (lroot) print*,'viscous force: nu is calculated for a mixture'
          lvisc_mixture=.true.
        case ('nu-spitzer')
          if (lroot) print*,'viscous force: temperature dependent nu'
          lpenc_requested(i_sij)=.true.
          lvisc_spitzer=.true.
        case ('none')
          ! do nothing
        case ('')
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
        if ( (lvisc_simplified.or.lvisc_rho_nu_const.or.&
             lvisc_sqrtrho_nu_const.or.lvisc_nu_const.or.&
             lvisc_nu_tdep.or.lvisc_nu_therm.or.&
             lvisc_mu_therm).and.nu==0.0) &
            call warning('initialize_viscosity', &
            'Viscosity coefficient nu is zero!')
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
        if (lvisc_nu_shock.and.nu_shock==0.0) &
            call fatal_error('initialize_viscosity', &
            'Viscosity coefficient nu_shock is zero!')
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
        endif
      endif
!
!  Shared variables.
!
      call put_shared_variable('lvisc_hyper3_nu_const_strict',lvisc_hyper3_nu_const_strict,ierr)
      call put_shared_variable('nu',nu,ierr)
      call get_shared_variable('lviscosity_heat',lviscosity_heat,ierr)
      if (ierr/=0) call stop_it("initialize_viscosity: " &
          // "problem getting shared var lviscosity_heat")
      call put_shared_variable('llambda_effect',llambda_effect,ierr)
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
          call get_shared_variable('etat_x',etat_x,ierr)
          if (ierr/=0) call fatal_error("initialize_viscosity","shared etat_x")
          call get_shared_variable('etat_y',etat_y,ierr)
          if (ierr/=0) call fatal_error("initialize_viscosity","shared etat_y")
          call get_shared_variable('etat_z',etat_z,ierr)
          if (ierr/=0) call fatal_error("initialize_viscosity","shared etat_z")
          call get_shared_variable('detat_x',detat_x,ierr)
          if (ierr/=0) call fatal_error("initialize_viscosity","shared detat_x")
          call get_shared_variable('detat_y',detat_y,ierr)
          if (ierr/=0) call fatal_error("initialize_viscosity","shared detat_y")
          call get_shared_variable('detat_z',detat_z,ierr)
          if (ierr/=0) call fatal_error("initialize_viscosity","shared detat_z")
          print*,'ipz,z(n),etat_z(n),detat_z(n)'
          do n=n1,n2
            print*,ipz,z(n),etat_z(n),detat_z(n)
          enddo
        endif
      endif
!
      call keep_compiler_quiet(lstarting)
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
      integer :: ierr
!
      if ((Lambda_V0==0).and.(Lambda_V1==0).and.(Lambda_H1==0)) &
        call warning('initialize_viscosity', &
            'You have chose llambda_effect=T but, all Lambda coefficients to be zero!')
      if ((Lambda_V0==0).and.((Lambda_V1/=0).or.(Lambda_H1==0))) &
        call warning('initialize_viscosity', &
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
        call fatal_error('initialize_viscosity(lambda)',&
            'default lambda_profile is uniform ! ')
      endselect
      lambda_V0t=lambda_V0*LV0_rprof(nx)
      lambda_V1t=lambda_V1*LV1_rprof(nx)
      lambda_V0b=lambda_V0*LV0_rprof(1)
      lambda_V1b=lambda_V1*LV1_rprof(1)
      call put_shared_variable('Lambda_V0t',Lambda_V0t,ierr)
      call put_shared_variable('Lambda_V1t',Lambda_V1t,ierr)
      call put_shared_variable('Lambda_V0b',Lambda_V0b,ierr)
      call put_shared_variable('Lambda_V1b',Lambda_V1b,ierr)
      call put_shared_variable('Lambda_H1',Lambda_H1,ierr)
      call put_shared_variable('LH1_rprof',LH1_rprof,ierr)
!
        endsubroutine initialize_lambda
!***********************************************************************
    subroutine read_viscosity_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      call keep_compiler_quiet(present(iostat))
!
    endsubroutine read_viscosity_init_pars
!***********************************************************************
    subroutine write_viscosity_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_viscosity_init_pars
!***********************************************************************
    subroutine read_viscosity_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=viscosity_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=viscosity_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_viscosity_run_pars
!***********************************************************************
    subroutine write_viscosity_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=viscosity_run_pars)
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
        idiag_dtnu=0; idiag_nu_LES=0; idiag_epsK=0; idiag_epsK_LES=0
        idiag_visc_heatm=0; idiag_meshRemax=0; idiag_Reshock=0
        idiag_nuD2uxbxm=0; idiag_nuD2uxbym=0; idiag_nuD2uxbzm=0
        idiag_nu_tdep=0; idiag_fviscm=0
        idiag_fviscmz=0; idiag_fviscmx=0; idiag_fviscmxy=0
        idiag_epsKmz=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) print*,'rprint_viscosity: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'nu_tdep',idiag_nu_tdep)
        call parse_name(iname,cname(iname),cform(iname),'fviscm',idiag_fviscm)
        call parse_name(iname,cname(iname),cform(iname),'fviscmin',idiag_fviscmin)
        call parse_name(iname,cname(iname),cform(iname),'fviscmax',idiag_fviscmax)
        call parse_name(iname,cname(iname),cform(iname),'nusmagm',idiag_nusmagm)
        call parse_name(iname,cname(iname),cform(iname),'nusmagmin',idiag_nusmagmin)
        call parse_name(iname,cname(iname),cform(iname),'nusmagmax',idiag_nusmagmax)
        call parse_name(iname,cname(iname),cform(iname),'dtnu',idiag_dtnu)
        call parse_name(iname,cname(iname),cform(iname),'nu_LES',idiag_nu_LES)
        call parse_name(iname,cname(iname),cform(iname),'visc_heatm', &
            idiag_visc_heatm)
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
      enddo
!
!  Check for those quantities for which we want z-averages
!
      do ixy=1,nnamexy
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'fviscmxy',idiag_fviscmxy)
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
!
      if ((lentropy.or.ltemperature) .and. &
          (lvisc_rho_nu_const .or. &
           lvisc_sqrtrho_nu_const .or. lvisc_nu_therm .or.&
           lvisc_nu_const .or. lvisc_nu_tdep .or. lvisc_nu_shock .or. &
           lvisc_nu_prof .or. lvisc_nu_profx .or. lvisc_spitzer .or. &
           lvisc_nu_profr .or. lvisc_nu_profr_powerlaw .or. &
           lvisc_nut_from_magnetic .or. lvisc_mu_therm))&
           lpenc_requested(i_TT1)=.true.
      if (lvisc_rho_nu_const.or.lvisc_sqrtrho_nu_const.or. &
          lvisc_nu_const .or. lvisc_nu_tdep .or. lvisc_nu_therm .or. &
          lvisc_nu_prof.or.lvisc_nu_profx.or.lvisc_spitzer .or. &
          lvisc_nu_profr.or.lvisc_nu_profr_powerlaw .or. &
          lvisc_nut_from_magnetic.or.lvisc_mu_therm) then
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
      if (lvisc_nu_profr) then
        if (lsphere_in_a_box.or.lspherical_coords) then
          lpenc_requested(i_r_mn)=.true.
        else
          lpenc_requested(i_rcyl_mn)=.true.
          if (lcylinder_in_a_box) lpenc_requested(i_rcyl_mn1)=.true.
        endif
      endif
      if (lvisc_nu_profr_powerlaw) lpenc_requested(i_rcyl_mn)=.true.
      if (lvisc_simplified .or. lvisc_rho_nu_const .or. &
          lvisc_sqrtrho_nu_const .or. lvisc_nu_const .or. lvisc_nu_tdep .or. &
          lvisc_smag_simplified .or. lvisc_smag_cross_simplified .or. &
          lvisc_nu_prof .or. lvisc_nu_profx .or. lvisc_spitzer .or. &
          lvisc_nu_profr_powerlaw .or. lvisc_nu_profr .or. &
          lvisc_nut_from_magnetic .or. lvisc_nu_therm .or. lvisc_mu_therm) &
          lpenc_requested(i_del2u)=.true.
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
      if (lvisc_rho_nu_const .or. lvisc_sqrtrho_nu_const .or. &
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
          lvisc_nut_from_magnetic.or.lvisc_nu_therm .or.  &
          lvisc_mu_therm .or. lvisc_spitzer) &
          lpenc_requested(i_sglnrho)=.true.
      if (lvisc_nu_const .and. lmagfield_nu) lpenc_requested(i_b2)=.true.
      if (lvisc_hyper3_nu_const) lpenc_requested(i_uij5glnrho)=.true.
      if (ldensity.and.lvisc_nu_shock) then
        lpenc_requested(i_graddivu)=.true.
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_gshock)=.true.
        lpenc_requested(i_divu)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
      if (llambda_effect) then
        lpenc_requested(i_uij)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
      if (lvisc_mixture) then
        lpenc_requested(i_graddivu)=.true.
        lpenc_requested(i_del2u)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_sglnrho)=.true.
        lpenc_requested(i_nu)=.true.
        lpenc_requested(i_gradnu)=.true.
        lpenc_requested(i_sgnu)=.true.
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
      if (lvisc_heat_as_aux) then
        lpenc_diagnos(i_visc_heat)=.true.
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_sij2)=.true.
      endif
      if ((idiag_meshRemax/=0.or.idiag_dtnu/=0).and.lvisc_nu_shock) then
        lpenc_diagnos(i_shock)=.true.
      endif
      if (idiag_Reshock/=0) lpenc_diagnos(i_shock)=.true.
      if (idiag_fviscmz/=0.or.idiag_fviscmx/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_sij)=.true.
      endif
      if (idiag_fviscmxy/=0) then
        lpenc_diagnos2d(i_rho)=.true.
        lpenc_diagnos2d(i_sij)=.true.
      endif
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
!
      use Deriv, only: der5i1j,der6
      use Diagnostics, only: max_mn_name, sum_mn_name
      use Interstellar, only: calc_snr_damping
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (nx,3) :: tmp,tmp2,gradnu,sgradnu
      real, dimension (nx) :: murho1,muTT,nu_smag,tmp3,tmp4,pnu
      real, dimension (nx) :: lambda_phi
!
      integer :: i,j,ju
!
      intent(inout) :: f,p
!
!  Viscous force and viscous heating are calculated here (for later use).
!
      p%fvisc=0.0
      if (lpencil(i_visc_heat)) p%visc_heat=0.0
      if (lfirst.and.ldt) then
        p%diffus_total=0.0
        p%diffus_total2=0.0
        p%diffus_total3=0.0
      endif
!
      if (lvisc_simplified) then
!
!  viscous force: nu*del2v
!  -- not physically correct (no momentum conservation), but
!  numerically easy and in most cases qualitatively OK
!
        p%fvisc=p%fvisc+nu*p%del2u
        if (lpencil(i_visc_heat)) then  ! Heating term not implemented
          if (headtt) then
            call warning('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_simplified')
          endif
        endif
!
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+nu
      endif
!
      if (lvisc_rho_nu_const) then
!
!  viscous force: mu/rho*(del2u+graddivu/3)
!  -- the correct expression for rho*nu=const
!
        murho1=nu*p%rho1  !(=mu/rho)
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i) + &
              murho1*(p%del2u(:,i)+1.0/3.0*p%graddivu(:,i))
        enddo
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*nu*p%sij2*p%rho1
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+murho1
      endif
!
      if (lvisc_sqrtrho_nu_const) then
!
!  viscous force: mu/rho*(del2u+graddivu/3)
!
        murho1=nu*sqrt(p%rho1)
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i) + &
              murho1*(p%del2u(:,i)+1.0/3.0*p%graddivu(:,i))
        enddo
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*murho1*p%sij2
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+murho1
      endif
!
      if (lvisc_mu_therm) then
!
!  viscous force: nu*sqrt(TT)/rho*(del2u+graddivu/3+2S.glnrho)
!
        muTT=nu*p%rho1*sqrt(exp(p%lnTT))
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i) + &
              muTT*(p%del2u(:,i)+1.0/3.0*p%graddivu(:,i))&
              + 2*muTT*p%sglnrho(:,i)
        enddo
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*muTT*p%sij2
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+muTT
      endif
!
      if (lvisc_spitzer) then
!
!  viscous force: nu*TT^2.5/rho*(del2u+graddivu/3+2S.glnrho)
!
        muTT=nu_spitzer*p%rho1*exp(2.5*p%lnTT)
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i) + &
              muTT*(p%del2u(:,i)+1.0/3.0*p%graddivu(:,i))&
              + 2*muTT*p%sglnrho(:,i)
        enddo
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*muTT*p%sij2
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+muTT
      endif
!
      if (lvisc_nu_therm) then
!
!  viscous force: nu*sqrt(TT)*(del2u+graddivu/3+2S.glnrho)
!  -- for numerical stability viscous force propto soundspeed in interstellar run
!
        muTT=nu*sqrt(exp(p%lnTT))
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
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+muTT
      endif
!
      if (lvisc_nu_const) then
!
!  viscous force: nu*(del2u+graddivu/3+2S.glnrho)
!  -- the correct expression for nu=const
!
        if (ldensity) then
          p%fvisc = p%fvisc + 2*nu*p%sglnrho+nu*(p%del2u + 1./3.*p%graddivu)
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
              call multsv_mn_add(1./sqrt(1.+p%b2/meanfield_nuB**2),p%fvisc+nu*(p%del2u+1./3.*p%graddivu),p%fvisc)
            endif
          else
            p%fvisc=p%fvisc+nu*(p%del2u+1.0/3.0*p%graddivu)
          endif
        endif
!
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*nu*p%sij2
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+nu
      endif
!
      if (lvisc_nu_tdep) then
!
!  viscous force: nu(t)*(del2u+graddivu/3+2S.glnrho)
!  -- the correct expression for nu=const
!
        nu_tdep=nu*(t-nu_tdep_t0)**nu_tdep_exponent
        p%fvisc=p%fvisc+2*nu_tdep*p%sglnrho+nu_tdep*(p%del2u+1./3.*p%graddivu)
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*nu_tdep*p%sij2
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+nu_tdep
      endif
!
      if (lvisc_mixture) then
!
!  Viscous force: nu*(del2u+graddivu/3+2S.glnrho)+2S.gradnu
!
        if (lpencil(i_sgnu)) call multmv(p%sij,p%gradnu,p%sgnu)
!
        if (ldensity) then
          do i=1,3
            p%fvisc(:,i)=2*p%nu*p%sglnrho(:,i) &
            +p%nu*(p%del2u(:,i)+1./3.*p%graddivu(:,i)) &
            +2*p%sgnu(:,i)
!
          enddo
        endif
!
!  Viscous heating and time step.
!
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*p%nu*p%sij2
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+p%nu
      endif
!
      if (lvisc_nu_profx.or.lvisc_nu_profr) then
!
!  Viscous force: nu(x)*(del2u+graddivu/3+2S.glnrho)+2S.gnu
!  -- here the nu viscosity depends on x; nu_jump=nu2/nu1
!        pnu = nu + nu*(nu_jump-1.)*step(abs(p%x_mn),xnu,widthnu)
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
          call fatal_error("","")
        endif
        pnu = nu + nu*(nu_jump-1.)*(step(tmp3,xnu ,widthnu) - &
                                    step(tmp3,xnu2,widthnu))
        tmp4=nu*(nu_jump-1.)*(der_step(tmp3,xnu ,widthnu)- &
                              der_step(tmp3,xnu2,widthnu))
!
!  Initialize gradnu before calculating it, otherwise gfortran complains.
!
        gradnu=0.
        call get_gradnu(tmp4,lvisc_nu_profx,&
                             lvisc_nu_profr,p,gradnu)
!  A routine for write_xprof should be written here
!       call write_xprof('visc',pnu)
        !gradnu(:,1) = nu*(nu_jump-1.)*der_step(tmp3,xnu,widthnu)
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
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+pnu
      endif
!
!  Radial viscosity profile from power law.
!
      if (lvisc_nu_profr_powerlaw) then
!
!  Viscous force: nu(x)*(del2u+graddivu/3+2S.glnrho)+2S.gnu
!  -- here the nu viscosity depends on r; nu=nu_0*(r/r0)^(-pnlaw)
!
        pnu = nu*(p%rcyl_mn/xnu)**(-pnlaw)
!  viscosity gradient
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
!
        call multmv(p%sij,gradnu,sgradnu)
        call multsv(pnu,2*p%sglnrho+p%del2u+1./3.*p%graddivu,tmp)
        p%fvisc=p%fvisc+tmp+2*sgradnu
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*pnu*p%sij2
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+pnu
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
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+pnu
      endif
!
!  viscous force: nu(z)*(del2u+graddivu/3+2S.glnrho)+2S.gnu
!  -- here the nu viscosity depends on z; nu_jump=nu2/nu1
!
      if (lvisc_nu_prof) then
        pnu = nu + nu*(nu_jump-1.)*step(p%z_mn,znu,-widthnu)
!  Write out viscosity z-profile (during first time step only)
        call write_zprof('visc',pnu)
        gradnu(:,1) = 0.
        gradnu(:,2) = 0.
        gradnu(:,3) = nu*(nu_jump-1.)*der_step(p%z_mn,znu,-widthnu)
        call multmv(p%sij,gradnu,sgradnu)
        call multsv(pnu,2*p%sglnrho+p%del2u+1./3.*p%graddivu,tmp)
        !tobi: The following only works with operator overloading for pencils
        !      (see sub.f90). Commented out because it seems to be slower.
        !p%fvisc=p%fvisc+2*pnu*p%sglnrho+pnu*(p%del2u+1./3.*p%graddivu) &
        !        +2*sgradnu
        p%fvisc=p%fvisc+tmp+2*sgradnu
        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*pnu*p%sij2
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+pnu
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
          if (lfirst.and.ldt) p%diffus_total=p%diffus_total+(nu_shock*p%shock)
          if (lpencil(i_visc_heat)) &
              p%visc_heat=p%visc_heat+nu_shock*p%shock*p%divu**2
        endif
      endif
!
      if (lvisc_hyper2_simplified) then
!
!  viscous force: nu_hyper2*de46v (not momentum-conserving)
!
        p%fvisc=p%fvisc-nu_hyper2*p%del4u
        if (lpencil(i_visc_heat)) then  ! Heating term not implemented
          if (headtt) then
            call warning('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_hyper2_simplified')
          endif
        endif
        if (lfirst.and.ldt) p%diffus_total2=p%diffus_total2+nu_hyper2
      endif
!
      if (lvisc_hyper3_simplified) then
!
!  viscous force: nu_hyper3*del6v (not momentum-conserving)
!
        p%fvisc=p%fvisc+nu_hyper3*p%del6u
        if (lpencil(i_visc_heat)) then  ! Heating term not implemented
          if (headtt) then
            call warning('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_hyper3_simplified')
          endif
        endif
        if (lfirst.and.ldt) p%diffus_total3=p%diffus_total3+nu_hyper3
      endif
!
      if (lvisc_hyper3_polar) then
!
! General way of coding an anisotropic hyperviscosity.
!
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
          if (lfirst.and.ldt) &
               p%diffus_total3=p%diffus_total3+nu_hyper3*pi4_1/dxyz_4
        enddo
      endif
!
      if (lvisc_hyper3_mesh) then
!
! Following Axel's hyper3_mesh for density
!
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
          if (lfirst.and.ldt) then
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
      if (lvisc_hyper3_rho_nu_const) then
!
!  viscous force: mu/rho*del6u
!
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
        if (lfirst.and.ldt) p%diffus_total3=p%diffus_total3+murho1
      endif
!
      if (lvisc_hyper3_rho_nu_const_symm) then
!
!  For tau_ij=d^5u_i/dx_j^5 + d^5u_j/dx_i^5
!  Viscous force: du/dt = mu/rho*{del6(u) + grad5[div(u)]}
!
        murho1=nu_hyper3*p%rho1  ! (=mu_hyper3/rho)
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i)+murho1*(p%del6u(:,i) + p%grad5divu(:,i))
        enddo
        if (lpencil(i_visc_heat)) then
          do i=1,3; do j=1,3
!  Dissipation is *not* positive definite.
            p%visc_heat=p%visc_heat + &
                0.5*nu_hyper3*(p%uij5(:,i,j)+p%uij5(:,j,i))*p%uij(:,i,j)
          enddo; enddo
        endif
        if (lfirst.and.ldt) p%diffus_total3=p%diffus_total3+murho1
      endif
!
      if (lvisc_hyper3_mu_const_strict) then
!
!  Viscous force:
!      du/dt = mu/rho*{del2(del2(del2(u))) + 1/3*del2(del2(grad(div(u))
!  This viscosity type requires HYPERVISC_STRICT =  hypervisc_strict_fft in
!  your Makefile.local.
!
        murho1=nu_hyper3*p%rho1  ! (=mu_hyper3/rho)
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i)+murho1*f(l1:l2,m,n,ihypvis-1+i)
        enddo
        if (lpencil(i_visc_heat)) then  ! Should be eps=2*mu*{del2[del2(S)]}^2
          if (headtt) then              ! (see Haugen & Brandenburg 2004 eq. 7)
            call warning('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_hyper3_mu_const_strict')
          endif
        endif
        if (lfirst.and.ldt) p%diffus_total3=p%diffus_total3+nu_hyper3
      endif
!
      if (lvisc_hyper3_nu_const_strict) then
!
!  Viscous force:
!      du/dt = 1/rho*div[2*rho*nu_3*S^(3)]
!  This viscosity type requires HYPERVISC_STRICT =  hypervisc_strict_fft in
!  your Makefile.local.
!
        do i=1,3
          p%fvisc(:,i)=p%fvisc(:,i)+nu_hyper3*f(l1:l2,m,n,ihypvis-1+i)
        enddo
        if (lpencil(i_visc_heat)) then
          if (headtt) then
            call warning('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_hyper3_nu_const_strict')
          endif
        endif
        if (lfirst.and.ldt) p%diffus_total3=p%diffus_total3+nu_hyper3
      endif
!
      if (lvisc_hyper3_rho_nu_const_aniso) then
!
!  viscous force: f_i = mu_i/rho*del6u
!  Used for non-cubic cells
!
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
        if (lfirst.and.ldt) p%diffus_total3=p%diffus_total3+&
             (nu_aniso_hyper3(1)*dx_1(l1:l2)**6 + &
              nu_aniso_hyper3(2)*dy_1(  m  )**6 + &
              nu_aniso_hyper3(3)*dz_1(  n  )**6)/ dxyz_6
!
      endif
!
      if (lvisc_hyper3_nu_const_aniso) then
!
!  viscous force: f_i = (nu_j.del6)u_i + nu_j.uij5.glnrho
!  Used for non-cubic cells
!
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
         if (lfirst.and.ldt) p%diffus_total3=p%diffus_total3+&
                 (nu_aniso_hyper3(1)*dx_1(l1:l2)**6 + &
                  nu_aniso_hyper3(2)*dy_1(  m  )**6 + &
                  nu_aniso_hyper3(3)*dz_1(  n  )**6)/ dxyz_6
!
      endif
!
      if (lvisc_hyper3_rho_nu_const_bulk) then
!
!  viscous force: mu/rho*d6uj/dx6
!
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
        if (lfirst.and.ldt) p%diffus_total3=p%diffus_total3+murho1
      endif
!
      if (lvisc_hyper3_nu_const) then
!
!  viscous force: nu_hyper3*(del6u+S.glnrho), where S_ij=d^5 u_i/dx_j^5
!
        p%fvisc=p%fvisc+nu_hyper3*(p%del6u+p%uij5glnrho)
        if (lpencil(i_visc_heat)) then  ! Heating term not implemented
          if (headtt) then
            call warning('calc_pencils_viscosity', 'viscous heating term '// &
              'is not implemented for lvisc_hyper3_nu_const')
          endif
        endif
        if (lfirst.and.ldt) p%diffus_total3=p%diffus_total3+nu_hyper3
      endif
!
      if (linterstellar.and.lvisc_snr_damp) then
        call calc_snr_damping(p)
      endif
!
!  viscous force: Handle damping at the core of SNRs
!
      if (lvisc_smag_simplified) then
!
!  viscous force: nu_smag*(del2u+graddivu/3+2S.glnrho)
!  where nu_smag=(C_smag*dxmax)**2*sqrt(2*SS)
!
        if (ldensity) then
!
! Find nu_smag
!
          nu_smag=(C_smag*dxmax)**2.*sqrt(2*p%sij2)
!
! Calculate viscous force
!
          call multsv_mn(nu_smag,p%sglnrho,tmp2)
          call multsv_mn(nu_smag,p%del2u+1./3.*p%graddivu,tmp)
          p%fvisc=p%fvisc+2*tmp2+tmp
          if (lpencil(i_visc_heat)) then
              p%visc_heat=p%visc_heat+2*nu_smag*p%sij2
          endif
          if (lfirst.and.ldt) p%diffus_total=p%diffus_total+nu_smag
        else
          if (lfirstpoint) print*, 'calc_viscous_force: '// &
              "ldensity better be .true. for ivisc='smagorinsky'"
        endif
      endif
!
      if (lvisc_smag_cross_simplified) then
!
!  viscous force: nu_smag*(del2u+graddivu/3+2S.glnrho)
!  where nu_smag=(C_smag*dxmax)**2*sqrt(S:J)
!
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
          if (lfirst.and.ldt) p%diffus_total=p%diffus_total+nu_smag
        else
          if (lfirstpoint) print*, 'calc_viscous_force: '// &
              "ldensity better be .true. for ivisc='smagorinsky'"
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
      if (ldiagnos) then
        if (idiag_nusmagm/=0)   call sum_mn_name(nu_smag,idiag_nusmagm)
        if (idiag_nusmagmin/=0) call max_mn_name(-nu_smag,idiag_nusmagmin,lneg=.true.)
        if (idiag_nusmagmax/=0) call max_mn_name(nu_smag,idiag_nusmagmax)
      endif
!
    endsubroutine calc_pencils_viscosity
!***********************************************************************
    subroutine calc_viscosity(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
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
    subroutine calc_viscous_heat(f,df,p,Hmax)
!
!  Calculate viscous heating term for right hand side of entropy equation.
!
!  20-nov-02/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar)    :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: Hmax
!
!  Add viscous heat (which has units of energy/mass) to the RHS
!  of the entropy (both with and without pretend_lnTT), or of
!  the temperature equation. Divide by cv if pretend_lnTT.
!
      if (lentropy) then
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
        df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) + p%rho*p%visc_heat
      endif
!
!  Calculate maximum heating (for time step constraint), so it is
!  only done on the first of the 3 substeps.
!
      if (lfirst .and. ldt) Hmax=Hmax+p%visc_heat
!
      call keep_compiler_quiet(f)
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
      use Sub, only: cross
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: nu_smag,Reshock
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
      if (lfirst.and.ldt) then
        diffus_nu =p%diffus_total *dxyz_2
        diffus_nu2=p%diffus_total2*dxyz_4
        if (ldynamical_diffusion .and. lvisc_hyper3_mesh) then
          diffus_nu3 = p%diffus_total3 * (abs(dline_1(:,1)) + abs(dline_1(:,2)) + abs(dline_1(:,3)))
        else
          diffus_nu3=p%diffus_total3*dxyz_6
        endif
      endif
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_nu_tdep/=0)  call sum_mn_name(spread(nu_tdep,1,nx),idiag_nu_tdep)
        if (idiag_fviscm/=0)   call sum_mn_name(p%fvisc,idiag_fviscm)
        if (idiag_fviscmin/=0) call max_mn_name(-p%fvisc,idiag_fviscmin,lneg=.true.)
        if (idiag_fviscmax/=0) call max_mn_name(p%fvisc,idiag_fviscmax)
        if (lvisc_smag_simplified) then
          if (ldensity) then
            nu_smag=(C_smag*dxmax)**2.*sqrt(2*p%sij2)
          endif
        endif
        if (lvisc_smag_cross_simplified) then
          if (ldensity) then
            nu_smag=(C_smag*dxmax)**2.*p%ss12
          endif
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
        if (idiag_epsK/=0) call sum_mn_name(p%visc_heat*p%rho,idiag_epsK)
!  Viscous heating for Smagorinsky viscosity.
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
      endif
!
!  2D-averages.
!
      if (l2davgfirst) then
        if (idiag_fviscmxy/=0) call zsum_mn_name_xy(-2.*p%rho*nu*( &
            p%uu(:,1)*p%sij(:,1,1)+p%uu(:,2)*p%sij(:,2,1)+ &
            p%uu(:,3)*p%sij(:,3,1)),idiag_fviscmxy)
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
        elseif (lvisc_rho_nu_const) then
          nu_pencil=nu*p%rho1
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
    subroutine dynamical_viscosity(umax)
!
!  Dynamically set viscosity coefficient given fixed mesh Reynolds number.
!
!  27-jul-11/ccyang: coded
!
      real, intent(in) :: umax
!
!  Hyper-viscosity coefficient
!
      if (nu_hyper3/=0.0) nu_hyper3 = pi5_1 * umax * dxmax**5 / re_mesh
      if (nu_hyper3_mesh/=0.0) nu_hyper3_mesh = pi5_1 * umax / re_mesh
!
    endsubroutine dynamical_viscosity
!***********************************************************************
    subroutine calc_lambda(p,div_lambda)
!
!  Calculates the lambda effect
!
!  20-apr-10/dhruba: coded
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
