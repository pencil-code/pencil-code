! $Id$
!
!  This module provide a way for users to specify custom
!  (i.e. not in the standard Pencil Code) physics, diagnostics etc.
!
!  The module provides a set of standard hooks into the Pencil-Code and
!  currently allows the following customizations:
!
!  Description                                     | Relevant function call
!  ---------------------------------------------------------------------------
!  Special variable registration                   | register_special
!    (pre parameter read)                          |
!  Special variable initialization                 | initialize_special
!    (post parameter read)                         |
!  Special variable finalization                   | finalize_special
!    (deallocation, etc.)                          |
!                                                  |
!  Special initial condition                       | init_special
!   this is called last so may be used to modify   |
!   the mvar variables declared by this module     |
!   or optionally modify any of the other f array  |
!   variables.  The latter, however, should be     |
!   avoided where ever possible.                   |
!                                                  |
!  Special term in the mass (density) equation     | special_calc_density
!  Special term in the momentum (hydro) equation   | special_calc_hydro
!  Special term in the energy equation             | special_calc_energy
!  Special term in the induction (magnetic)        | special_calc_magnetic
!     equation                                     |
!                                                  |
!  Special equation                                | dspecial_dt
!    NOT IMPLEMENTED FULLY YET - HOOKS NOT PLACED INTO THE PENCIL-CODE
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 2
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED infl_phi; infl_dphi; gphi(3)
!***************************************************************
!
! HOW TO USE THIS FILE
! --------------------
!
! Change the line above to
!   lspecial = .true.
! to enable use of special hooks.
!
! The rest of this file may be used as a template for your own
! special module.  Lines which are double commented are intended
! as examples of code.  Simply fill out the prototypes for the
! features you want to use.
!
! Save the file with a meaningful name, eg. geo_kws.f90 and place
! it in the $PENCIL_HOME/src/special directory.  This path has
! been created to allow users ot optionally check their contributions
! in to the Pencil-Code SVN repository.  This may be useful if you
! are working on/using the additional physics with somebodyelse or
! may require some assistance from one of the main Pencil-Code team.
!
! To use your additional physics code edit the Makefile.local in
! the src directory under the run directory in which you wish to
! use your additional physics.  Add a line with all the module
! selections to say something like:
!
!   SPECIAL=special/geo_kws
!
! Where geo_kws it replaced by the filename of your new module
! upto and not including the .f90
!
module Special
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include '../special.h'
!
!
! Declare index of new variables in f array (if any).
!
!  integer :: iinfl_phi=0, iinfl_dphi=0, iinfl_hubble=0, iinfl_lna=0, Ndiv=100
  integer :: iinfl_phi=0, iinfl_dphi=0, iinfl_lna=0, Ndiv=100
  integer :: iinfl_rho_chi=0, iinfl_rho_rad=0
  real :: ncutoff_phi=1., infl_v=.1
  real :: axionmass=1.06e-6, axionmass2, ascale_ini=1.
  real :: phi0=.44, dphi0=-1.69e-7, c_light_axion=1., lambda_axion=0., eps=.01
  real :: amplphi=.1, ampldphi=.0, kx_phi=1., ky_phi=0., kz_phi=0., phase_phi=0., width=.1, offset=0.
  real :: initpower_phi=0.,  cutoff_phi=0.,  initpower2_phi=0.
  real :: initpower_dphi=0., cutoff_dphi=0., initpower2_dphi=0.
  real :: kgaussian_phi=0.,kpeak_phi=0., kgaussian_dphi=0., kpeak_dphi=0.
  real :: relhel_phi=0.
  real :: ddotam, a2rhopm, a2rhopm_all, a2rhom, a4rhom, a2rhom_all, a4rhom_all, rhom, rhom_all
  real :: edotbm, edotbm_all, e2m, e2m_all, b2m, b2m_all, a2rhophim, a4rhophim, a2rhophim_all, a4rhophim_all
  real :: a2rhopphim, a2rhopphim_all
  real :: sigE1m_all_nonaver, sigB1m_all_nonaver,sigEm_all,sigBm_all,sigEm_all_diagnos,sigBm_all_diagnos
  real :: a2rhogphim, a2rhogphim_all
  real :: a2, a21, Hscript
  real :: Hscript0=0., scale_rho_chi_Heqn=1., scale_rho_rad_Heqn=1., rho_chi_init=0.
  real :: cdt_rho_chi=1., cdt_phi=1e-2
  real :: amplee_BD_prefactor=0., deriv_prefactor_ee=-1.
  real :: echarge=.0, echarge_const=.303
  real :: count_eb0_all=0., rad_heating=0., ascale_heat=0., ascale_heat_off=0., heating
  real :: aphimax=0., aphimax2=0. !PAR_DOC: maximum a value above which the phi potential is quenched.
  real :: Gamma_phi0=impossible, Gamma_phi !PAR_DOC: damping factor for phi above aphimax
  real :: Gamma_phi_exp=3.        !PAR_DOC: scale factor exponent on Gamma_phi
  real :: rhophim_crit=1e-21      !PAR_DOC: minimum phi
  real :: wstate, wstate_aver     !PAR_DOC: critical w (EoS) value (slightly below 1/3)
  real :: wstate_crit=0.333333333 !PAR_DOC: critical w (EoS) value (1/3)
  real :: wstate_tolerance=0.     !PAR_DOC: tolerance w (EoS) value
  real :: wstate_prev=0.          !PAR_DOC: value of wstate in the previous iteration.
  real :: a4rhophim_crit=0.       !PAR_DOC: critical value of a4rhophim below lsolve_phi=F is set.
!
  real, target :: ddotam_all
  real, pointer :: alpf, eta
  real, pointer :: sigE_prefactor, sigB_prefactor, mass_chi
  real, dimension (nx) :: dt1_special
  logical :: lcompute_dphi0=.true.
  logical :: lbackreact_infl=.true., lem_backreact=.true., lzeroHubble=.false.
  logical :: lscale_tobox=.true., ldt_backreact_infl=.true., lconf_time=.true.
  logical :: lskip_projection_phi=.false., lvectorpotential=.false., lflrw=.false.
  logical :: lrho_chi=.false., lno_noise_phi=.false., lno_noise_dphi=.false.
  logical :: lrho_rad=.false.            !PAR_DOC: radiation from inflaton decay
  logical :: lrho_rad_apply=.true.       !PAR_DOC: radiation from inflaton decay, and also applied to df(l1:l2,m,n,iinfl_dphi)
  logical :: lrho_rad_apply2=.true.      !PAR_DOC: radiation from inflaton decay, and also applied to df(l1:l2,m,n,iinfl_dphi)
  logical :: lrho_chi_corrected=.true.   !PAR_DOC: for backward compatibility; when false, we use the wrong scale factor in the rho_chi equation
  logical :: lrho_chi_inhom=.false.      !PAR_DOC: inhomogeneous heating
  logical :: ldefine_a2rhophi_with_Vpotential=.true.  !PAR_DOC: define a2rhophi with Vpotential
  logical :: lsolve_for_phi=.true.       !PAR_DOC: whether we still want to solve for phi
  logical :: lsolve_for_phi2=.true.      !PAR_DOC: whether we still want to solve for phi
  logical :: lsolve_for_phi_switch=.true. !PAR_DOC: switch must be on for automatically switching off the phi solver.
  logical :: lsolve_for_phi_always=.true. !PAR_DOC: misnomer: is now used for switching off the phi solver: is now used for switching off the phi solver
  logical :: lwstate_crit=.true.         !PAR_DOC: lwstate_crit switch (would put phi=0, is false by default)
  logical :: lwstate_crit_old=.false.    !PAR_DOC: lwstate_crit_old (to restore the old wstate criterion used in the autotest)
  logical :: lheating=.false.            !PAR_DOC: heating criterion
  logical :: lheating_always=.false.     !PAR_DOC: heating criterion, set to true once lheating=T.
  logical :: lheating_keep_on=.false.    !PAR_DOC: heating criterion
  logical :: ldefine_a2rhopm_without_Vpotential=.false.    !PAR_DOC: should be false to have correct results
  logical :: la2rhop_wrong_factor=.false. !PAR_DOC: should be false to have correct results; kept for backward compatibility
  logical :: lappy_BD_k1D_factor=.false. !PAR_DOC: apply $k_1^D$ factor in the Bunch-Davies initial condition (NOTE typo in name!)
  logical :: lapply_BD_kNy_factor=.false. !PAR_DOC: apply $1/N^(D/2)$ factor in the Bunch-Davies initial condition.
  logical :: linv_BD=.true.              !PAR_DOC: apply forward transform in the Bunch-Davies initial condition.
  logical :: lswitch_toMHD_when_nophi=.true. !PAR_DOC: switch to MHD when phi evolution is turned off.
  logical, pointer :: lphi_hom, lphi_linear_regime, lnoncollinear_EB, lnoncollinear_EB_aver
  logical, pointer :: lcollinear_EB, lcollinear_EB_aver, lmass_suppression
  logical, pointer :: lallow_bprime_zero
  logical, pointer :: ladvance_ee
  !Whether the sums needed for the ODE and rhs advancement are done in the together in the same kernel as the rhs
  !advancement. Benchmarks seem to suggest that combining them is indeed more performant.
  !This approach is however strictly approximative since we effectively take the value of Hscript from the preceeding substep
  logical :: lcombine_prep_ode_right_with_rhs = .false.
  character (len=labellen) :: Vprime_choice='quadratic', Hscript_choice='default'
  character (len=labellen) :: heating_choice='Hscript_max'
  character (len=labellen) :: solve_phi_criterion='wstate'
  character (len=labellen), dimension(ninit) :: initspecial='nothing'
  character (len=50) :: echarge_type='const', init_rho_chi='zero', init_rho_rad='zero'
!
  namelist /special_init_pars/ &
      initspecial, phi0, dphi0, axionmass, eps, ascale_ini, &
      lcompute_dphi0, lem_backreact, &
      c_light_axion, lambda_axion, amplphi, ampldphi, lno_noise_phi, lno_noise_dphi, &
      kx_phi, ky_phi, kz_phi, phase_phi, width, offset, &
      initpower_phi, initpower2_phi, cutoff_phi, kgaussian_phi, kpeak_phi, &
      initpower_dphi, initpower2_dphi, cutoff_dphi, kpeak_dphi, &
      ncutoff_phi, lscale_tobox, Hscript0, Hscript_choice, infl_v, lflrw, &
      lrho_chi, scale_rho_chi_Heqn, scale_rho_rad_Heqn, amplee_BD_prefactor, deriv_prefactor_ee, &
      lrho_rad, init_rho_rad, Gamma_phi0, lconf_time, &
      echarge_type, init_rho_chi, rho_chi_init, lrho_chi_inhom, rhophim_crit, &
      wstate_crit, lwstate_crit, lwstate_crit_old, wstate_tolerance, &
      lsolve_for_phi_always, lsolve_for_phi_switch, &
      heating_choice, ldefine_a2rhopm_without_Vpotential, la2rhop_wrong_factor, &
      lappy_BD_k1D_factor, lapply_BD_kNy_factor, linv_BD
!
  namelist /special_run_pars/ &
      initspecial, phi0, dphi0, axionmass, eps, ascale_ini, &
      lbackreact_infl, lem_backreact, c_light_axion, lambda_axion, Vprime_choice, &
      lzeroHubble, ldt_backreact_infl, Ndiv, Hscript0, Hscript_choice, infl_v, &
      lflrw, lrho_chi, scale_rho_chi_Heqn, scale_rho_rad_Heqn, echarge_type, &
      cdt_rho_chi, cdt_phi, &
      lrho_rad, lrho_rad_apply, lrho_rad_apply2, lrho_chi_corrected, &
      lrho_chi_inhom, ldefine_a2rhophi_with_Vpotential, &
      rad_heating, ascale_heat, ascale_heat_off, aphimax, Gamma_phi0, lconf_time, rhophim_crit, &
      wstate_crit, lwstate_crit, lwstate_crit_old, wstate_tolerance, &
      lsolve_for_phi_always, lsolve_for_phi_switch, &
      heating_choice, lheating_keep_on, lcombine_prep_ode_right_with_rhs, &
      lswitch_toMHD_when_nophi, Gamma_phi_exp, a4rhophim_crit, solve_phi_criterion
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_phim=0       ! DIAG_DOC: $\left<\phi\right>$
  integer :: idiag_phi2m=0      ! DIAG_DOC: $\left<\phi^2\right>$
  integer :: idiag_phirms=0     ! DIAG_DOC: $\left<\phi^2\right>^{1/2}$
  integer :: idiag_dphim=0      ! DIAG_DOC: $\left<\phi'\right>$
  integer :: idiag_dphi2m=0     ! DIAG_DOC: $\left<(\phi')^2\right>$
  integer :: idiag_dphirms=0    ! DIAG_DOC: $\left<(\phi')^2\right>^{1/2}$
  integer :: idiag_dtphi=0      ! DIAG_DOC: $dt/cdtphi$
  integer :: idiag_Hscriptm=0   ! DIAG_DOC: $\left<{\cal a*H}\right>$
  integer :: idiag_ascale=0     ! DIAG_DOC: $a$
  integer :: idiag_lnam=0       ! DIAG_DOC: $\left<\ln a\right>$
  integer :: idiag_ddotam=0     ! DIAG_DOC: $a''/a$
  integer :: idiag_a2rhopm=0    ! DIAG_DOC: $a^2 (\rho+p)$
  integer :: idiag_a2rhom=0     ! DIAG_DOC: $a^2 \rho$
  integer :: idiag_a2rhophim=0  ! DIAG_DOC: $a^2 \rho_\phi$
  integer :: idiag_a4rhophim=0  ! DIAG_DOC: $a^4 \rho_\phi$
  integer :: idiag_a2rhogphim=0 ! DIAG_DOC: $0.5 <grad \phi^2>$
  integer :: idiag_rho_chi=0    ! DIAG_DOC: $\rho_\chi$
  integer :: idiag_rho_rad=0    ! DIAG_DOC: $\rho_\mathrm{rad}$
  integer :: idiag_sigEma=0     ! DIAG_DOC: $\rho_\chi$
  integer :: idiag_sigBma=0     ! DIAG_DOC: $\rho_\chi$
  integer :: idiag_count_eb0a=0 ! DIAG_DOC: $f_\mathrm{EB0}$
  integer :: idiag_heating=0    ! DIAG_DOC: $\theta_\mathrm{heat}$
  integer :: idiag_wstate=0     ! DIAG_DOC: $w_\mathrm{state}$
  integer :: idiag_wstate_aver=0 ! DIAG_DOC: $\langle w_\mathrm{state}\rangle$
  integer :: idiag_Gamma_phi=0  ! DIAG_DOC: $\langle w_\mathrm{state}\rangle$
!
  integer :: enum_hscript_choice = 0
  integer :: enum_vprime_choice = 0
  integer :: enum_echarge_type = 0

  real :: a2rhom_all_diagnos, a2rhopm_all_diagnos, a2rhophim_all_diagnos, a4rhophim_all_diagnos
  real :: a2rhogphim_all_diagnos, ddotam_all_diagnos

  contains
!****************************************************************************
    subroutine register_special
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
!
      use FArrayManager
      use SharedVariables, only: put_shared_variable
!
      if (lroot) call svn_id( &
          "$Id$")
!
      call farray_register_pde('infl_phi',iinfl_phi)
      call farray_register_pde('infl_dphi',iinfl_dphi)
!
     if (lflrw) call farray_register_ode('infl_lna',iinfl_lna)
     if (lrho_chi) call farray_register_ode('infl_rho_chi',iinfl_rho_chi)
     if (lrho_rad) call farray_register_ode('infl_rho_rad',iinfl_rho_rad)
!
!  for power spectra, it is convenient to use ispecialvar and
!
      ispecialvar=iinfl_phi
      ispecialvar2=iinfl_dphi
!
      call put_shared_variable('ddotam',ddotam_all,caller='register_backreact_infl')
      call put_shared_variable('Hscript',Hscript)
      call put_shared_variable('e2m_all',e2m_all)
      call put_shared_variable('b2m_all',b2m_all)
      call put_shared_variable('sigEm_all',sigEm_all,caller='register_backreact_infl')
      call put_shared_variable('sigBm_all',sigBm_all,caller='register_backreact_infl')
      call put_shared_variable('echarge',echarge,caller='register_backreact_infl')
      call put_shared_variable('lrho_chi',lrho_chi)
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!
      use SharedVariables, only: get_shared_variable
      use FArrayManager, only: farray_index_by_name_ode
      use Messages, only: warning
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: iLCDM_lna
!
      if (lflrw) then
        iLCDM_lna=farray_index_by_name_ode('iLCDM_lna')
        if (iLCDM_lna>0) call fatal_error('initialize_special', 'there is a conflict with iLCDM_lna')
      endif
!
!  set axionmass**2
!
      axionmass2=axionmass**2
      aphimax2=aphimax**2
!
      if (lmagnetic .and. lem_backreact) then
!
!  The following variables are defined in disp_current, but in principle,
!  the backreaction module can also run without it. In that case, these
!  variables need to be defined here.
!
        if (iex>0) then
          call get_shared_variable('alpf',alpf,caller='initialize_backreact_infl')
          call get_shared_variable('lphi_hom',lphi_hom)
          call get_shared_variable('lphi_linear_regime',lphi_linear_regime)
          call get_shared_variable('sigE_prefactor',sigE_prefactor)
          call get_shared_variable('sigB_prefactor',sigB_prefactor)
          call get_shared_variable('lcollinear_EB',lcollinear_EB)
          call get_shared_variable('lcollinear_EB_aver',lcollinear_EB_aver)
          call get_shared_variable('lnoncollinear_EB',lnoncollinear_EB)
          call get_shared_variable('lnoncollinear_EB_aver',lnoncollinear_EB_aver)
          call get_shared_variable('lmass_suppression',lmass_suppression)
          call get_shared_variable('lallow_bprime_zero',lallow_bprime_zero)
          call get_shared_variable('ladvance_ee',ladvance_ee)
          call get_shared_variable('mass_chi',mass_chi)
        else
          if (.not.associated(lphi_hom)) allocate(lphi_hom, lphi_linear_regime, ladvance_ee, &
            lcollinear_EB, lnoncollinear_EB, lcollinear_EB_aver, lnoncollinear_EB_aver)
          lphi_hom=.true.
          lphi_linear_regime=.true.
          ladvance_ee=.false.
          lcollinear_EB=.false.
          lnoncollinear_EB=.false.
          lcollinear_EB_aver=.false.
          lnoncollinear_EB_aver=.false.
          call get_shared_variable('eta',eta,caller='initialize_backreact_infl')
        endif
      else
        if (.not.associated(alpf)) allocate(alpf, lphi_hom, lphi_linear_regime, &
          sigE_prefactor, sigB_prefactor, lcollinear_EB, lcollinear_EB_aver, &
          lnoncollinear_EB, lnoncollinear_EB_aver, lmass_suppression, &
          lallow_bprime_zero, mass_chi)
        alpf=0.
        lphi_hom=.false.
        lphi_linear_regime=.false.
        sigE_prefactor=0.
        sigB_prefactor=0.
        lcollinear_EB=.false.
        lcollinear_EB_aver=.false.
        lnoncollinear_EB=.false.
        lnoncollinear_EB_aver=.false.
        lmass_suppression=.false.
        lallow_bprime_zero=.false.
        mass_chi=0.
      endif
!
!  Redundancy checks. To be removed when these logicals are removed.
!
      if (lwstate_crit_old) call warning("initialize_special","lwstate_crit_old ONLY for 0d-tests/reheating")
      if (.not.lwstate_crit) call fatal_error("initialize_special","lwstate_crit disabled")
      if (.not.lsolve_for_phi2) call fatal_error("initialize_special","lsolve_for_phi2 disabled")
      if (.not.lsolve_for_phi_always) call fatal_error("initialize_special","lsolve_for_phi_always disabled")
      if (.not.lsolve_for_phi_switch) call fatal_error("initialize_special","lsolve_for_phi_switch disabled")
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Initcond, only: gaunoise, sinwave_phase, hat, power_randomphase_hel, power_randomphase, bunch_davies
      use Mpicomm, only: mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: Vpotential, Hubble_ini, infl_gam, amplphi_BD, amplee_BD, deriv_prefactor
      integer :: j
      real :: lnascale
!
      intent(inout) :: f
!
      do j=1,ninit
        select case (initspecial(j))
          case ('nothing'); if (lroot) print*,'init_special: nothing'
          case ('phi=sinkx')
            f(:,:,:,iinfl_phi)=f(:,:,:,iinfl_phi) &
              +spread(spread(amplphi*sin(kx_phi*x),2,my),3,mz)
          case ('phi=tanhkx')
            f(:,:,:,iinfl_phi)=f(:,:,:,iinfl_phi) &
              +spread(spread(.5*amplphi*(1.+tanh(kx_phi*(x-offset))),2,my),3,mz)
          ! sine-Gordon solution
          case ('phi=atan_exp_kx')
            infl_gam=1./sqrt(1.-infl_v**2)
            f(:,:,:,iinfl_phi)=f(:,:,:,iinfl_phi) &
              +spread(spread(4.*amplphi*atan(exp(infl_gam*kx_phi*(x-offset))),2,my),3,mz)
            f(:,:,:,iinfl_dphi)=f(:,:,:,iinfl_dphi)+spread(spread( &
              -4.*amplphi*kx_phi*infl_gam*infl_v*exp(infl_gam*kx_phi*(x-offset)) &
              /(exp(2.*infl_gam*kx_phi*(x-offset))+1.) &
              ,2,my),3,mz)
          case ('nophi')
            Vpotential=.5*axionmass2*phi0**2
            dphi0=0.
            tstart=-sqrt(3./(8.*pi))/(ascale_ini*sqrt(Vpotential))
            t=tstart
            Hubble_ini=sqrt(8.*pi/3.*(.5*dphi0**2+.5*axionmass2*phi0**2*ascale_ini**2))
            lnascale=log(ascale_ini)
            if (lroot .and. lflrw) f_ode(iinfl_lna)=lnascale
!
          case ('default')
            Vpotential=.5*axionmass2*phi0**2
!
!  Hubble_ini is here based on the standard (non-reduced) Planck mass.
!
            Hubble_ini=sqrt(8.*pi/3.*(.5*axionmass2*phi0**2*ascale_ini**2))
            ! dphi0=-ascale_ini*sqrt(2*eps/3.*Vpotential)
            if (lcompute_dphi0) dphi0=-sqrt(1/(12.*pi))*axionmass*ascale_ini
            ! dphi0=-sqrt(1/(12.*pi))*axionmass*ascale_ini
            ! dphi0=-sqrt(16*pi/3)*axionmass*ascale_ini
!
!  Initial time.
!
            if (lconf_time) &
              tstart=-1./(ascale_ini*Hubble_ini)
            t=tstart
!
            lnascale=log(ascale_ini)
            f(:,:,:,iinfl_phi)   =f(:,:,:,iinfl_phi)   +phi0
            f(:,:,:,iinfl_dphi)  =f(:,:,:,iinfl_dphi)  +dphi0
            if (lroot .and. lflrw) then
              f_ode(iinfl_lna)   =lnascale
              a2                 =exp(f_ode(iinfl_lna))**2
              Hscript            =Hubble_ini/exp(lnascale)
!
!  Should not be needed.
!
!              f(iinfl_hubble)   =Hscript
            endif
          case ('gaussian-noise')
            call gaunoise(amplphi,f,iinfl_phi)
          case ('sinwave-phase')
            !call sinwave_phase(f,iinfl_phi,amplphi,kx_phi,ky_phi,kz_phi,phase_phi)
            !f(:,:,:,iinfl_phi)=tanh(f(:,:,:,iinfl_phi)/width)
            call hat(amplphi,f,iinfl_phi,width,kx_phi,ky_phi,kz_phi)
            f(:,:,:,iinfl_phi)=f(:,:,:,iinfl_phi)+offset
          case ('phi_power_randomphase')
            call power_randomphase_hel(amplphi,initpower_phi,initpower2_phi, &
              cutoff_phi,ncutoff_phi,kpeak_phi,f,iinfl_phi,iinfl_phi, &
              relhel_phi,kgaussian_phi, lskip_projection_phi, lvectorpotential, &
              lscale_tobox, lpower_profile_file=.false., lno_noise=lno_noise_phi)
          case ('dphi_power_randomphase')
            call power_randomphase_hel(ampldphi,initpower_dphi,initpower2_dphi, &
              cutoff_dphi,ncutoff_phi,kpeak_dphi,f,iinfl_dphi,iinfl_dphi, &
              relhel_phi,kgaussian_phi, lskip_projection_phi, lvectorpotential, &
              lscale_tobox, lpower_profile_file=.false., lno_noise=lno_noise_dphi)
!
!  For Bunch-Davies, the amplitude Hubble_ini is used.
!  We apply this optionally here also to the gauge field.
!  The amplitudes amplphi and amplee_BD_prefactor should be unity in theory.
!
          case ('Bunch-Davies')
            if (lroot) print*,'Hubble_ini=',Hubble_ini
            amplphi_BD=amplphi*Hubble_ini
            deriv_prefactor=1.
            call bunch_davies(f,iinfl_phi,iinfl_phi,iinfl_dphi,iinfl_dphi, &
              amplphi_BD,kpeak_phi,deriv_prefactor, &
              lappy_BD_k1D_factor=lappy_BD_k1D_factor, &
              lapply_BD_kNy_factor=lapply_BD_kNy_factor, linv=linv_BD)
            if (amplee_BD_prefactor/=0.) then
              amplee_BD=amplee_BD_prefactor*Hubble_ini
              if (iex>0) then
                deriv_prefactor=deriv_prefactor_ee
                call bunch_davies(f,iax,iaz,iex,iez, &
                  amplee_BD,kpeak_phi,deriv_prefactor, &
                  lappy_BD_k1D_factor=lappy_BD_k1D_factor, &
                  lapply_BD_kNy_factor=lapply_BD_kNy_factor, linv=linv_BD)
              else
                deriv_prefactor=0.
              endif
            endif
          case default
            call fatal_error("init_special","no such initspecial: "//trim(initspecial(j)))
        endselect
      enddo
!
!  initial condition for energy density of charged particles
!
      if (lroot .and. lrho_chi) then
        select case (init_rho_chi)
          case ('zero'); f_ode(iinfl_rho_chi)=0.
          case ('given'); f_ode(iinfl_rho_chi)=rho_chi_init
          case default
            call fatal_error("init_special","no such init_rho_chi: "//trim(init_rho_chi))
        endselect
      endif
!
!  initial condition for energy density of radiation
!
      if (lroot .and. lrho_rad) then
        select case (init_rho_rad)
          case ('zero'); f_ode(iinfl_rho_rad)=0.
          !case ('given'); f_ode(iinfl_rho_rad)=rho_rad_init
          case default
            call fatal_error("init_special: No such init_rho_rad: ", trim(init_rho_rad))
        endselect
      endif
!
!  Better default value based on alpf and axionmass
!  Particle physics group Chin. Phys. C 38, 090001, (2014).
!  But we take larger values of 10^{-8} or so.
!
      if (Gamma_phi0==impossible) then
        Gamma_phi0=alpf**2*axionmass**3/(64.*pi)
      endif
!
      call mpibcast_real(a2)
      call mpibcast_real(Hscript)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
!  pencil for gradient of phi
!
      ! alberto: gphi pencil does not seem to be used anywhere, right?
      ! lpenc_requested(i_gphi)=.true.
!
!  Magnetic field needed for Maxwell stress
!
      if (lmagnetic .and. lem_backreact) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_el)=.true.
        ! alberto: pencil p%e2 does not seem to be used
        ! if (lrho_chi .or. lnoncollinear_EB .or. lnoncollinear_EB_aver .or. &
        !   lcollinear_EB .or. lcollinear_EB_aver) lpenc_requested(i_e2)=.true.
      endif
      !  Call pencils phi and dphi
      lpenc_requested(i_infl_phi)=.true.
      lpenc_requested(i_infl_dphi)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Special pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  24-nov-04/tony: coded
!
      use Sub, only: grad
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
! infl_phi
      if (lpencil(i_infl_phi)) p%infl_phi=f(l1:l2,m,n,iinfl_phi)
!
! infl_dphi
      if (lpencil(i_infl_dphi)) p%infl_dphi=f(l1:l2,m,n,iinfl_dphi)
!
! gphi
      if (lpencil(i_gphi)) call grad(f,iinfl_phi,p%gphi)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine get_Hscript_and_a2(Hscript,a2rhom_all)
!
      real, intent(out) :: Hscript
      real, intent(in), optional :: a2rhom_all
!
!  Choice of prescription for Hscript (or H for cosmic time).
!  alberto: to be changed, default to 'set' with Hscript0=0 and remove lzeroHubble
!           as it trivially corresponds to new default choice
!           old 'default' should correspond to 'friedmann'
!
      select case (Hscript_choice)
        case ('default')
          if (lconf_time) then
            Hscript=sqrt((8.*pi/3.)*a2rhom_all)
          else
            Hscript=sqrt((8.*pi/3.)*a2rhom_all*a21)
          endif
          if (lgpu) call get_a2
        case ('set')
          Hscript=Hscript0
          a2=1.
          a21=1./a2
        case ('friedmann')
          Hscript=sqrt((8.*pi/3.)*a2rhom_all)
          if (lgpu) call get_a2
        case default
          call fatal_error("get_Hscript_and_a2","no such Hscript_choice: "//trim(Hscript_choice))
      endselect

!  Possibility of turning off evolution of scale factor and Hubble parameter
!  By default, lzeroHubble=F, so we use the calculation from above.
!
      if (lzeroHubble) then
        a2=1.
        a21=1./a2
        Hscript=0.
      endif
!
    endsubroutine get_Hscript_and_a2
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  The entire module could be renamed to Klein-Gordon or Scalar field equation.
!  Calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  Several precalculated Pencils of information are passed for
!  efficiency.
!
!  06-oct-03/tony: coded
!   2-nov-21/axel: first set of equations coded
!
      use Diagnostics, only: sum_mn_name, max_mn_name, save_name
      use Sub, only: dot_mn, dot2_mn, del2, grad
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: gphi
      real, dimension (nx) :: Vprime, Vpotential, a2rhophi, a4rhophi
      real, dimension (nx) :: tmp, del2phi, gphi2
      real :: pref_Vprime=1., pref_Hubble=2., pref_del2=1., pref_alpf, pref_Gamma=impossible
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!
!  Choice of different potentials.
!  For the 1-cos profile, -Vprime (on the rhs) enters with -sin().
!
! alberto: changed to use the pencils p%infl_phi and p%infl_dphi
      select case (Vprime_choice)
        case ('quadratic')
          Vprime=axionmass2*p%infl_phi
          Vpotential=.5*axionmass2*p%infl_phi**2
        case ('quartic')
          Vprime=axionmass2*p%infl_phi+(lambda_axion/6.)*p%infl_phi**3
          Vpotential=axionmass2*p%infl_phi+(lambda_axion/6.)*p%infl_phi**3  !(to be corrected)
        case ('cos-profile')
          Vprime=axionmass2*lambda_axion*sin(lambda_axion*p%infl_phi)
          Vpotential=axionmass2*lambda_axion*sin(lambda_axion*p%infl_phi)  !(to be corrected)
        case default
          call fatal_error("dspecial_dt","no such Vprime_choice: "//trim(Vprime_choice))
      endselect
!
!  Current choice of temporal form of Gamma_phi. Heating is currently instantaenous.
!  to rename lheating -> lheating_phi. The switch lheating_always is false by default
!  and set to true after the first time lheating is true and if lheating_keep_on is true.
!
      if (lheating .or. lheating_always) then
        Gamma_phi=Gamma_phi0
      else
        Gamma_phi=0.
      endif
!
!  Update df.
!  dphi/dt = psi
!  dpsi/dt = - ...
!
!  alberto: determine prefactors for the different terms beforehand
!  Allowed for quenching factor on pref_Vprime to limit excessive oscillations
!  (not yet done for conformal time).
!
      if (lconf_time) then
        pref_alpf=a21
        pref_Vprime=a2
        pref_Gamma=ascale
      ! alberto: for cosmic time, coefficient of Hscript should be 3
      else
        pref_Hubble=3.; pref_Vprime=1.; pref_del2=a21
        pref_alpf=a21**2
        pref_Gamma=1.
      endif
!
!  In the following, Hscript means H if cosmic time is used.
!  dphi/dt = dphi. So lrho_rad_apply is needed in the inhomogeneous case (but not lrho_rad_apply2).
!  We have d2phi/dt^2 = ... (2*Hscript+Gamma_phi)*dphi/dt, so the timestep constraint is
!  dt < 1/(2*Hscript+Gamma_phi).
!
      if (lsolve_for_phi) then
        df(l1:l2,m,n,iinfl_phi)=df(l1:l2,m,n,iinfl_phi)+p%infl_dphi
        if (lrho_rad .and. lrho_rad_apply) then
          df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi) - &
              (pref_Hubble*Hscript+pref_Gamma*Gamma_phi)*p%infl_dphi-pref_Vprime*Vprime
        else
          df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi) - &
              (pref_Hubble*Hscript)*p%infl_dphi-pref_Vprime*Vprime
        endif
      else
!
!  Here we use the possibility of switching off the rhs of the phi evolution
!
        df(l1:l2,m,n,iinfl_phi)=0.
        df(l1:l2,m,n,iinfl_dphi)=0.
      endif
!
!  Here we include the Gamma_phi contribution to the rhs of the lnrho equation.
!  drho/dt = ... + a*Gam_phi*a2rhophi/a^6.
!  This is independent of whether or not lrho_rad=T; which is just for testing.
!
      if (lrho_chi_inhom .and. lheating) then
        if (ldensity) then
          call grad(f,iinfl_phi,gphi)
          call dot2_mn(gphi,gphi2)
          a2rhophi=0.5*p%infl_dphi**2+0.5*gphi2+a2*Vpotential
          tmp=Gamma_phi*a2rhophi*ascale**Gamma_phi_exp
          if (ldensity_nolog) then
            df(l1:l2,m,n,irho)=df(l1:l2,m,n,irho)+tmp
          else
            df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)+tmp*p%rho1
          endif
        else
          call fatal_error("dspecial_dt: ","density must be true")
        endif
      endif
!
!  speed of light term
!
      if (lsolve_for_phi) then
        if (c_light_axion/=0. .and. .not. lphi_hom) then
          call del2(f,iinfl_phi,del2phi)
          df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi) + &
              c_light_axion**2*pref_del2*del2phi
          ! if (lconf_time) then
          !   df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)+c_light_axion**2*del2phi
          ! else
          !   df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)+c_light_axion**2*a21*del2phi
          !endif
        endif
!
!  magnetic terms, add (alpf/a^2)*(E.B) to dphi'/dt equation
!
        if (lmagnetic .and. lem_backreact) then
          if (lphi_hom .and. .not. lphi_linear_regime) then
            df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)+pref_alpf*alpf*edotbm_all
          endif
!
!  Compute E.B only when displacement current is included.
!  Note that alpf does not (currently) exist in MHD.
!
          if (.not. lphi_hom) then
            if (iex>0) then
              call dot_mn(p%el,p%bb,tmp)
            else
              call dot_mn(p%jj,p%bb,tmp)
              tmp=eta*tmp
            endif
            df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)+pref_alpf*alpf*tmp
          endif
          ! if (lconf_time) then
          !   if (lphi_hom) then
          !     if (.not. lphi_linear_regime) &
          !       df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)+alpf*edotbm_all*a21
          !   else
          !     call dot_mn(p%el,p%bb,tmp)
          !     df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)+alpf*tmp*a21
          !   endif
          ! else
          !   if (lphi_hom) then
          !     if (.not. lphi_linear_regime) &
          !       df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)+alpf*edotbm_all*a21**2
          !   else
          !     call dot_mn(p%el,p%bb,tmp)
          !     df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)+alpf*tmp*a21**2
          !   endif
          ! endif
        endif
      endif
!
!  Total contribution to the timestep.
!  If Ndiv=0 is set, we compute the timestep based on the frequency of phi and
!  based on the advective timestep from the Alfven speed for rho_chi.
!  vA=B/sqrt(rho_chi), so dt=C_M*dx/vA. In practice, C_M (=cdt_rho_chi) can be 20.
!  This timestep constraint should not be used if hydrodynamics is evolved.
!  If it is used, it should be based on cobformal time using comoving B and rho.
!
      if (lfirst .and. ldt .and. ldt_backreact_infl) then
        if (Ndiv==0.) then
          if (lsolve_for_phi) then
            if (dimensionality==0) then
              dt1_special=axionmass*sqrt(a2)/cdt_phi
            else
              advec2=max(advec2,axionmass2*a2*dxyz_2/cdt_phi**2)
              dt1_special=0.
            endif
          endif
!
!  Additional constraint from vA=B/sqrt(rho_chi), but this is only relevant
!  when ldensity=F, because otherwise the standard Alfven constaint applies.
!
          if (lrho_chi) then
            advec2=max(advec2,a21**2*b2m_all/f_ode(iinfl_rho_chi)*dxyz_2/cdt_rho_chi**2)
          else
            call fatal_error("dspecial_dt", "lrho_chi must be .true. when Ndiv=0")
          endif
        else
          dt1_special = Ndiv*abs(Hscript)
        endif
        dt1_max=max(dt1_max,dt1_special)
      endif
!
!  Diagnostics
!
      call calc_diagnostics_special(f,p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine dspecial_dt_ode
!
      use SharedVariables, only: get_shared_variable
      use Diagnostics , only: 
!
!  Here is the Friedmann equation. For cosmic time, Hscript means H.
!
      if (lgpu) call read_sums_from_device
      call get_Hscript_and_a2(Hscript,a2rhom_all)
      if (lflrw) df_ode(iinfl_lna)=df_ode(iinfl_lna)+Hscript
!
!  Energy density of the charged particles.
!  This is currently only done for <sigE>*<E^2>, and not for <sigE*E^2>.
!
      if (lrho_chi) then
        if (lnoncollinear_EB .or. lnoncollinear_EB_aver .or. &
          lcollinear_EB .or. lcollinear_EB_aver) then
          if (lrho_chi_corrected) then
            df_ode(iinfl_rho_chi)=df_ode(iinfl_rho_chi)-4.*Hscript*f_ode(iinfl_rho_chi) &
              +(sigEm_all*e2m_all+sigBm_all*edotbm_all)/ascale**4
          else
            df_ode(iinfl_rho_chi)=df_ode(iinfl_rho_chi)-4.*Hscript*f_ode(iinfl_rho_chi) &
              +(sigEm_all*e2m_all+sigBm_all*edotbm_all)/ascale**3
          endif
!
!  In a test example, assume that a certain level of heating, rad_heating,
!  is turned on after ascale has exceeded a critical level, ascale_heat.
!
     !  elseif (ascale_heat>0) then
     !    if (ascale_heat_off>0) then
     !      heating=rad_heating*.25*(1.+tanh(ascale-ascale_heat))*(1.-tanh(ascale-ascale_heat_off))
     !    else
     !      heating=rad_heating*.5*(1.+tanh(ascale-ascale_heat))
     !    endif
     !    df_ode(iinfl_rho_chi)=df_ode(iinfl_rho_chi)-4.*Hscript*f_ode(iinfl_rho_chi)+heating/ascale**4
        else
          df_ode(iinfl_rho_chi)=df_ode(iinfl_rho_chi)-4.*Hscript*f_ode(iinfl_rho_chi)
        endif
      endif
!
!  Energy density of radiation from decaying inflaton field.
!  Since a2rhophim_all includes a^2, and we want to have the physical
!  energy density, so we divide by a^2, but in conformal time,
!  we also need to multiply by a, so we divide by a.
!  This "heating" enters the physical (not comoving) rho_rad!
!
      if (lrho_rad) then
        if (lconf_time) then
          heating=Gamma_phi*a2rhophim_all/ascale
        else
          heating=Gamma_phi*a2rhophim_all*a21
        endif
        df_ode(iinfl_rho_rad)=df_ode(iinfl_rho_rad)-4.*Hscript*f_ode(iinfl_rho_rad)+heating
      endif
!
!  Diagnostics
!
      if (.not. lmultithread) then
        sigEm_all_diagnos = sigEm_all
        sigBm_all_diagnos = sigBm_all
        call calc_ode_diagnostics_special(f_ode)
      endif
!
    endsubroutine dspecial_dt_ode
!***********************************************************************
    subroutine calc_ode_diagnostics_special(f_ode)
!
      use Diagnostics 
      
      real, dimension(n_odevars), intent(in) :: f_ode
      real :: rho_chi, rho_rad, lnascale
      real :: Hscript_diagnos
!
!  Set rho_chi for diagnostics
!
      if (lrho_chi) then
        rho_chi=f_ode(iinfl_rho_chi)
      else
        rho_chi=0.
      endif
!
!  Set rho_rad for diagnostics
!
      if (lrho_rad) then
        rho_rad=f_ode(iinfl_rho_rad)
      else
        rho_rad=0.
      endif
!
      if (ldiagnos) then
        call get_Hscript_and_a2(Hscript_diagnos,a2rhom_all_diagnos)
        if (lflrw) lnascale=f_ode(iinfl_lna)
        call save_name(Hscript_diagnos,idiag_Hscriptm)
        call save_name(ddotam_all_diagnos,idiag_ddotam)
        call save_name(lnascale,idiag_lnam)
        call save_name(ascale,idiag_ascale)
        call save_name(a2rhopm_all_diagnos,idiag_a2rhopm)
        call save_name(a2rhom_all_diagnos,idiag_a2rhom)
        call save_name(a2rhophim_all_diagnos,idiag_a2rhophim)
        call save_name(a4rhophim_all_diagnos,idiag_a4rhophim)
        call save_name(a2rhogphim_all_diagnos,idiag_a2rhogphim)
        call save_name(rho_chi,idiag_rho_chi)
        call save_name(rho_rad,idiag_rho_rad)
        call save_name(sigEm_all_diagnos,idiag_sigEma)
        call save_name(sigBm_all_diagnos,idiag_sigBma)
        if (lnoncollinear_EB_aver .or. lcollinear_EB_aver) &
          call save_name(count_eb0_all,idiag_count_eb0a)
        call save_name(heating,idiag_heating)
        call save_name(wstate,idiag_wstate)
        call save_name(wstate_aver,idiag_wstate_aver)
        call save_name(Gamma_phi,idiag_Gamma_phi)
!
      endif
!
    endsubroutine calc_ode_diagnostics_special
!***********************************************************************
    subroutine calc_diagnostics_special(f,p)
!
      use Diagnostics

      real, dimension(mx,my,mz,mfarray) :: f
      type(pencil_case) :: p

      call keep_compiler_quiet(f)
! alberto: changed to use the pencils p%infl_phi and p%infl_dphi
      if (ldiagnos) then
        call sum_mn_name(p%infl_phi,idiag_phim)
        if (idiag_phi2m/=0) call sum_mn_name(p%infl_phi**2,idiag_phi2m)
        if (idiag_phirms/=0) call sum_mn_name(p%infl_phi**2,idiag_phirms,lsqrt=.true.)
        call sum_mn_name(p%infl_dphi,idiag_dphim)
        if (idiag_dphi2m/=0) call sum_mn_name(p%infl_dphi**2,idiag_dphi2m)
        if (idiag_dphirms/=0) call sum_mn_name(p%infl_dphi**2,idiag_dphirms,lsqrt=.true.)
        call max_mn_name(sqrt(advec2)/cdt_phi,idiag_dtphi,l_dt=.true.)
      endif
!
    endsubroutine calc_diagnostics_special
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  Reads and registers print parameters relevant to special.
!
!  06-oct-03/tony: coded
!
      use Diagnostics, only: parse_name
!
      integer :: iname
      logical :: lreset,lwrite

      call keep_compiler_quiet(lwrite)
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_phim=0; idiag_phi2m=0; idiag_phirms=0
        idiag_dphim=0; idiag_dphi2m=0; idiag_dphirms=0; idiag_dtphi=0
        idiag_Hscriptm=0; idiag_lnam=0; idiag_ascale=0; idiag_ddotam=0
        idiag_a2rhopm=0; idiag_a2rhom=0; idiag_a4rhophim=0; idiag_a2rhophim=0
        idiag_a2rhogphim=0; idiag_rho_chi=0; idiag_rho_rad=0; idiag_sigEma=0
        idiag_sigBma=0; idiag_count_eb0a=0; idiag_heating=0; idiag_wstate=0
        idiag_wstate_aver=0; idiag_Gamma_phi=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'phim',idiag_phim)
        call parse_name(iname,cname(iname),cform(iname),'phi2m',idiag_phi2m)
        call parse_name(iname,cname(iname),cform(iname),'phirms',idiag_phirms)
        call parse_name(iname,cname(iname),cform(iname),'dphim',idiag_dphim)
        call parse_name(iname,cname(iname),cform(iname),'dphi2m',idiag_dphi2m)
        call parse_name(iname,cname(iname),cform(iname),'dphirms',idiag_dphirms)
        call parse_name(iname,cname(iname),cform(iname),'dtphi',idiag_dtphi)
        call parse_name(iname,cname(iname),cform(iname),'Hscriptm',idiag_Hscriptm)
        call parse_name(iname,cname(iname),cform(iname),'lnam',idiag_lnam)
        call parse_name(iname,cname(iname),cform(iname),'ascale',idiag_ascale)
        call parse_name(iname,cname(iname),cform(iname),'ddotam',idiag_ddotam)
        call parse_name(iname,cname(iname),cform(iname),'a2rhopm',idiag_a2rhopm)
        call parse_name(iname,cname(iname),cform(iname),'a2rhom',idiag_a2rhom)
        call parse_name(iname,cname(iname),cform(iname),'a2rhophim',idiag_a2rhophim)
        call parse_name(iname,cname(iname),cform(iname),'a4rhophim',idiag_a4rhophim)
        call parse_name(iname,cname(iname),cform(iname),'a2rhogphim',idiag_a2rhogphim)
        call parse_name(iname,cname(iname),cform(iname),'rho_chi',idiag_rho_chi)
        call parse_name(iname,cname(iname),cform(iname),'rho_rad',idiag_rho_rad)
        call parse_name(iname,cname(iname),cform(iname),'sigEma',idiag_sigEma)
        call parse_name(iname,cname(iname),cform(iname),'sigBma',idiag_sigBma)
        call parse_name(iname,cname(iname),cform(iname),'count_eb0a',idiag_count_eb0a)
        call parse_name(iname,cname(iname),cform(iname),'heating',idiag_heating)
        call parse_name(iname,cname(iname),cform(iname),'wstate',idiag_wstate)
        call parse_name(iname,cname(iname),cform(iname),'wstate_aver',idiag_wstate_aver)
        call parse_name(iname,cname(iname),cform(iname),'Gamma_phi',idiag_Gamma_phi)
      enddo
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_echarge
!
      real :: energy_scale
!
!  Choice of echarge prescription.
!
      if (lnoncollinear_EB .or. lnoncollinear_EB_aver &
        .or. lcollinear_EB .or. lcollinear_EB_aver) then
        select case (echarge_type)
          case ('const')
            echarge=echarge_const
          case ('erun')
            energy_scale=(.5*e2m_all+.5*b2m_all)**.25/ascale
            echarge=1./sqrt(1./.35**2+41./(48.*pi**2)*log(mass_zboson/energy_scale))
        endselect
      else
        echarge=echarge_const
      endif
!
    endsubroutine get_echarge
!***********************************************************************
    subroutine get_sigE_and_B
!
      real :: boost, gam_EB, eprime, bprime, jprime1
      real :: mass_suppression_fact
      real :: sigE1m_all,sigB1m_all
!
!  Compute sigE and sigB from sigE1 and sigB1.
!  The mean conductivities are also needed in the local cases,
!  because they are used in the calculation of rho_chi.
!
      if (lnoncollinear_EB .or. lnoncollinear_EB_aver &
        .or. lcollinear_EB .or. lcollinear_EB_aver) then
        if (lnoncollinear_EB_aver) then
          boost=sqrt((e2m_all-b2m_all)**2+4.*edotbm_all**2)
          gam_EB=sqrt21*sqrt(1.+(e2m_all+b2m_all)/boost)
          eprime=sqrt21*sqrt(e2m_all-b2m_all+boost)
          bprime=sqrt21*sqrt(b2m_all-e2m_all+boost)*sign(1.,edotbm_all)
          if (eprime/=0. .and. bprime/=0.) then
            jprime1=1./(6.*pi**2)*eprime*abs(bprime)/tanh(pi*abs(bprime)/eprime)
            sigE1m_all=abs(jprime1)*eprime/(gam_EB*boost)
            sigB1m_all=abs(jprime1)*edotbm_all/(eprime*gam_EB*boost)
            count_eb0_all=0.
          else
            sigE1m_all=0.
            sigB1m_all=0.
            count_eb0_all=1.
          endif
!
!  Similarly for collinear case.
!  Mass suppression is currently defined for the collinear case only.
!
        elseif (lcollinear_EB_aver) then
          eprime=sqrt(e2m_all)
          bprime=sqrt(b2m_all)
          if (eprime/=0. .and. bprime/=0.) then
            sigE1m_all=1./(6.*pi**2)*bprime/tanh(pi*abs(bprime)/eprime)
            sigB1m_all=0.
            count_eb0_all=0.
          else
            sigE1m_all=0.
            sigB1m_all=0.
            count_eb0_all=1.
          endif
          if (lmass_suppression) then
            mass_suppression_fact=exp(-pi*mass_chi**2/(Chypercharge**onethird*echarge*eprime))
            sigE1m_all=sigE1m_all*mass_suppression_fact
          endif
        else
          sigE1m_all = sigE1m_all_nonaver
          sigB1m_all = sigB1m_all_nonaver
        endif
!
!  Apply Chypercharge, echarge, and Hscript universally for aver and nonaver.
!
        sigEm_all=sigE_prefactor*Chypercharge*echarge**3*sigE1m_all/Hscript
        sigBm_all=sigB_prefactor*Chypercharge*echarge**3*sigB1m_all/Hscript
      endif
!
    endsubroutine get_sigE_and_B
!***********************************************************************
    subroutine prep_rhs_special
!
!  1-aug-25/TP: coded
!
      call get_Hscript_and_a2(Hscript,a2rhom_all)
      call get_echarge
      call get_sigE_and_B
!
    endsubroutine prep_rhs_special
!***********************************************************************
    subroutine get_a2
!
      real :: lnascale

      if (lflrw) then
        lnascale=f_ode(iinfl_lna)
        ascale=exp(lnascale)
      endif
      a2=ascale**2
      a21=1./a2
!
    endsubroutine get_a2
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  Possibility to modify the f array after the boundaries are
!  communicated.
!
!  06-jul-06/tony: coded
!
      use Mpicomm, only: mpireduce_sum, mpiallreduce_sum, mpibcast_real
      use Sub, only: dot2_mn, grad, curl, dot_mn
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real :: sigE1m, sigB1m, rho_rad, Hscript_prev=0.
!
! TP: to avoid code duplication could this function not be combined with the copy of it in
!     klein_gordon.f90? We could make an appropriate module and call it from there
!
!  If requested, calculate here <dphi**2+gphi**2+(4./3.)*(E^2+B^2)/a^2>.
!  This needs to be done on all processors, because otherwise ascale
!  is not known on all processors.
!
      call get_a2
      call mpibcast_real(a2)
      call mpibcast_real(a21)
!
!  Here we use the possibility of switching off the phi evolution by setting phi=dphi=0.
!
      if (.not. lsolve_for_phi) then
        f(:,:,:,iinfl_phi)=0.
        f(:,:,:,iinfl_dphi)=0.
        if (lswitch_toMHD_when_nophi .and. iex>0) ladvance_ee=.false.
      endif
!
!  In the following loop, go through all penciles and add up results to get e2m, etc.
!
      ddotam=0.; a2rhopm=0.; a2rhom=0.; rhom=0; e2m=0; b2m=0; edotbm=0
      a2rhophim=0.; a4rhophim=0.; a2rhopphim=0.; a2rhogphim=0.; sigE1m=0.; sigB1m=0.
!
!  In the following, sum over all mn pencils.
!
      do n=n1,n2
      do m=m1,m2
        call prep_ode_right(f,sigE1m,sigB1m)
      enddo
      enddo
!
      rhom=rhom/nwgrid
      a2rhopm=a2rhopm/nwgrid
      a2rhom=a2rhom/nwgrid
      a2rhophim=a2rhophim/nwgrid
      a4rhophim=a4rhophim/nwgrid
      a2rhopphim=a2rhopphim/nwgrid
      a2rhogphim=a2rhogphim/nwgrid
      ddotam=(four_pi_over_three/nwgrid)*ddotam
      if (lphi_hom .or. lrho_chi .or. lnoncollinear_EB .or. lnoncollinear_EB_aver) then
        edotbm=edotbm/nwgrid
        call mpiallreduce_sum(edotbm,edotbm_all)
      endif
!
!  Schwinger conductivities:
!
      if (lrho_chi .or. lnoncollinear_EB .or. lnoncollinear_EB_aver &
        .or. lcollinear_EB .or. lcollinear_EB_aver) then
        e2m=e2m/nwgrid
        b2m=b2m/nwgrid
        call mpiallreduce_sum(e2m,e2m_all)
        call mpiallreduce_sum(b2m,b2m_all)
!
!  The following is not done for averages. This is because sigE1m and sigB1m
!  are calculated in their own section further below in the same routine.
!
        if (lrho_chi .or. lnoncollinear_EB .or. lcollinear_EB) then
          sigE1m=sigE1m/nwgrid
          sigB1m=sigB1m/nwgrid
          call mpiallreduce_sum(sigE1m,sigE1m_all_nonaver)
          call mpiallreduce_sum(sigB1m,sigB1m_all_nonaver)
        endif
      endif
!
!  Some use mpiallreduce_sum and others mpireduce_sum.
!  When the result is needed on all processors, use mpiallreduce_sum.
!  But we additionally use "call mpibcast_real(rhom_all)",
!  so maybe "all" is not needed?
!
      call mpiallreduce_sum(rhom,rhom_all)
      call mpireduce_sum(a2rhopm,a2rhopm_all)
      call mpiallreduce_sum(a2rhom,a2rhom_all)
      call mpireduce_sum(a2rhophim,a2rhophim_all)
      call mpireduce_sum(a4rhophim,a4rhophim_all)
      call mpireduce_sum(a2rhopphim,a2rhopphim_all)
      call mpireduce_sum(a2rhogphim,a2rhogphim_all)
      call mpiallreduce_sum(ddotam,ddotam_all)
      a2rhom_all_diagnos     = a2rhom_all
      a2rhopm_all_diagnos    = a2rhopm_all
      a2rhophim_all_diagnos  = a2rhophim_all
      a4rhophim_all_diagnos  = a4rhophim_all
      a2rhogphim_all_diagnos = a2rhogphim_all
      ddotam_all_diagnos     = ddotam_all
!
!  Set rho_rad for diagnostics
!
      if (lrho_rad) then
        rho_rad=f_ode(iinfl_rho_rad)
      else
        rho_rad=0.
      endif
!
!  Get Hscript and a2rhom_all.
!
      if (lroot .and. lflrw) call get_Hscript_and_a2(Hscript,a2rhom_all)
!
!  Broadcast to other processors, and each processor uses put_shared_variable
!  to get the values to other subroutines.
!
      call mpibcast_real(Hscript)
      call mpibcast_real(e2m_all)
      call mpibcast_real(b2m_all)
      call mpibcast_real(a2rhophim_all)
      call mpibcast_real(a4rhophim_all)
      call mpibcast_real(a2rhopphim_all)
!
!  Compute rhom, which is needed for wstate.
!
      if (ldensity) then
        call mpibcast_real(rhom_all)
        wstate=(a2rhopphim_all*a21+onethird*rhom)/(a2rhophim_all*a21+rhom)
      else
        wstate=(a2rhopphim_all*a21+onethird*rho_rad)/(a2rhophim_all*a21+rho_rad)
      endif
!
!  Alternatitives for deciding when to solve for phi: either when
!  wstate is not yet reached, or when a4rhophim is still big enough.
!
      if (lsolve_for_phi) then
        select case (solve_phi_criterion)
          case ('rhophi')
            lsolve_for_phi=a4rhophim > a4rhophim_crit
            if (lroot .and. .not. lsolve_for_phi) print*,'rhophi criterion activated'
          case ('wstate')
            if (lwstate_crit_old) then
              lsolve_for_phi=(wstate<wstate_crit)
              if (lroot .and. .not. lsolve_for_phi) print*,'OLD wstate criterion activated'
            else 
              if (lsolve_for_phi_always) then
                lsolve_for_phi=abs(wstate-wstate_prev) > wstate_tolerance
                wstate_prev=wstate
                lsolve_for_phi_always=lsolve_for_phi
              endif
              if (lroot .and. .not. lsolve_for_phi) print*,'wstate criterion activated'
            endif
          case default
            call fatal_error("special_after_boundary: No such solve_phi_criterion: ",trim(solve_phi_criterion))
        endselect
      endif
!
!  Alternatitives for deciding when to turn on heating,
!  i.e., when the end of inflation occurs.
!
      select case (heating_choice)
        case ('Hscript_max')
          if (.not. lheating_always) then
            lheating=Hscript<Hscript_prev
            Hscript_prev=Hscript
            lheating_always=lheating .and. lheating_keep_on
          endif
        case ('aphimax2')
          lheating=a2>aphimax2
        case default
          call fatal_error("special_after_boundary: No such heating_choice: ",trim(heating_choice))
      endselect
!
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine prep_ode_right(f,sigE1m,sigB1m)
!
      use Sub, only: dot2_mn, grad, curl, dot_mn, cross_mn
!
!  11-mar-26/axel: added accumulation of a2rho_phi (V added here a bit later)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, intent(inout) :: sigE1m,sigB1m
      real, dimension (nx,3) :: el, bb, gphi
      real, dimension (nx) :: e2, b2, gphi2, dphi, a2rhop, a2rho, a2rhophi, a4rhophi
      real, dimension (nx) :: a2rhopphi
      real, dimension (nx) :: ddota, phi, Vpotential, edotb, sigE1, sigB1
      real, dimension (nx) :: boost, gam_EB, eprime, bprime, jprime1
!
!  if requested, calculate here <dphi**2+gphi**2+(4./3.)*(E^2+B^2)/a^2>
!  rhop is purely an output quantity
!  It is called a2rhom because rhom=a2rhom/a^2.
!
      phi=f(l1:l2,m,n,iinfl_phi)
      dphi=f(l1:l2,m,n,iinfl_dphi)
!
!  For conformal time, rho is defined by .5*dphi^2/a^2,
!  but for cosmic time, we have just .5*dphi^2, so a2rho=.5*dphi^2*a^2.
!
      if (lconf_time) then
        a2rho=0.5*dphi**2
        if (la2rhop_wrong_factor) then
          a2rhop=dphi**2
        else
          a2rhop=0.5*dphi**2
        endif
      else
        a2rho=0.5*dphi**2*a2
        a2rhop=0.5*dphi**2*a2
      endif
!
!  a2rhop is for pressure.
!
!AB: removed following line
      !a2rhop=dphi**2
      ! if (lphi_hom) then
      !   a2rhop=dphi**2
      !   a2rho=0.5*dphi**2
      !   a2rhophim=a2rhophim+sum(a2rho)
      ! else
      if (.not. lphi_hom) then
        call grad(f,iinfl_phi,gphi)    !MR: the ghost zones are not necessarily updated!!!
        ! alberto: this function is called from special_after_boundary so shouldn't have the ghost zones updated?
        !TP: Not necessarily, since if learly_finalize=.false. it can be that the halo exchange has not finished
        !    here yet.
        call dot2_mn(gphi,gphi2)
        a2rhogphim=a2rhogphim+sum(0.5*gphi2)
!
!AB: probably mistake
!
        if (la2rhop_wrong_factor) then
          a2rhop=a2rhop+onethird*gphi2
        else
          a2rhop=a2rhop-0.5*onethird*gphi2
        endif
        a2rho=a2rho+0.5*gphi2
      endif
!
!  Set a2rhophim for later accummulation.
!  This is still local.
!
      a2rhophi=a2rho
      a4rhophi=a2rho*a2
      a2rhopphi=a2rhop
!
!  Note the .5*fourthird factor in front of (e2+b2)*a21, but that is
!  just for rhop, which is output quantity.
!
      if (lem_backreact) then
        call curl(f,iaa,bb)          !MR: the ghost zones are not necessarily updated!!! AB: Yes, but should be ok within the domain)
        call dot2_mn(bb,b2)
!
!  In MHD, when we don't have the electric field, we can use -uxB for now.
!  This method is currently not used.
!
        if (iex/=0) then
          el=f(l1:l2,m,n,iex:iez)
        else
          el=0.
        endif
        call dot2_mn(el,e2)
!
!AB: probably mistake
!
        if (la2rhop_wrong_factor) then
          a2rhop=a2rhop+(.5*fourthird)*(e2+b2)*a21
        else
          a2rhop=a2rhop+(.5*onethird)*(e2+b2)*a21
        endif
        if (.not. lphi_linear_regime) a2rho=a2rho+.5*(e2+b2)*a21
      endif
!
!  option to take the inhomogeneous rho instead
!  Here, in the expression a2rho, rho is not comoving, but the rho from f(l1:l2,m,n,ilnrho) is comoving.
!
      if (lrho_chi) then
        if (lrho_chi_inhom) then
          if (ldensity) then
            if (ldensity_nolog) then
              a2rho=a2rho+scale_rho_chi_Heqn/a2*f(l1:l2,m,n,irho)
            else
              a2rho=a2rho+scale_rho_chi_Heqn/a2*exp(f(l1:l2,m,n,ilnrho))
            endif
          else
            call fatal_error("backreact_infl special_after_boundary: No such Vprime_choice: ","density must be true")
          endif
        else
          a2rho=a2rho+scale_rho_chi_Heqn*a2*f_ode(iinfl_rho_chi)
        endif
      endif
!
!  Do the same for rho_rad. But this should not be applied in the inhomogeneous case when heating is on.
!
      if (lrho_rad .and. lrho_rad_apply2) then
        if (lrho_chi_inhom) then
          call fatal_error("prep_ode_right","lrho_rad_apply2 must be true if ldensity")
        else
          a2rho=a2rho+scale_rho_rad_Heqn*a2*f_ode(iinfl_rho_rad)
        endif
      endif
!
!  ldefine_a2rhopm_without_Vpotential would be *true* for reproducing the old runs.
!
      if (ldefine_a2rhopm_without_Vpotential) &
        a2rhopm=a2rhopm+sum(a2rhop)
!
!  Choice of different potentials
!
      select case (Vprime_choice)
        case ('quadratic')  ; Vpotential=.5*axionmass2*phi**2
        case ('quartic')    ; Vpotential=axionmass2*phi+(lambda_axion/6.)*phi**3  !(to be corrected)
        case ('cos-profile'); Vpotential=axionmass2*lambda_axion*sin(lambda_axion*phi)  !(to be corrected)
        case default
          call fatal_error("special_after_boundary: No such Vprime_choice: ",trim(Vprime_choice))
      endselect
!
!  compute ddotam = a"/a (needed for GW module)
!
      ddota=-dphi**2+4.*a2*Vpotential
      if (.not. lphi_hom) ddota=ddota-gphi2
      ! if (lphi_hom) then
      !   ddota=-dphi**2+4.*a2*Vpotential
      ! else
      !   ddota=-dphi**2-gphi2+4.*a2*Vpotential
      ! endif
      ddotam=ddotam+sum(ddota)
      a2rho=a2rho+a2*Vpotential
      if (.not. ldefine_a2rhopm_without_Vpotential) &
        a2rhop=a2rhop-a2*Vpotential
      if (ldefine_a2rhophi_with_Vpotential) &
        a2rhophi=a2rhophi+a2*Vpotential
!
!  Compute pressure for phi field (without gauge field):
!
      a2rhopphi=a2rhopphi-a2*Vpotential
!
!  Compute rhom, which is needed for wstate.
!
      if (ldensity) then
        if (ldensity_nolog) then
          rhom=rhom+sum(f(l1:l2,m,n,irho))
        else
          rhom=rhom+sum(exp(f(l1:l2,m,n,ilnrho)))
        endif
      endif
!
!  Compute volume average of both a2rho and a2rho_phi
!  Here, we begin to sum over a2rhophi, so this the last place where the local one was updated.
!
      a2rhom=a2rhom+sum(a2rho)
      if (.not. ldefine_a2rhopm_without_Vpotential) &
        a2rhopm=a2rhopm+sum(a2rhop)
      a2rhophim=a2rhophim+sum(a2rhophi)
      a4rhophim=a4rhophim+sum(a4rhophi)
      a2rhopphim=a2rhopphim+sum(a2rhopphi)
!
!  Compute electromagnetic averages.
!
      if (lmagnetic .and. lem_backreact) then
        if (lphi_hom .or. lrho_chi .or. lnoncollinear_EB .or. lnoncollinear_EB_aver &
                                   .or. lcollinear_EB .or. lcollinear_EB_aver) then
          call dot_mn(el,bb,edotb)
          edotbm=edotbm+sum(edotb)
!
!  Repeat calculation of sigE and sigB. Do this first without
!  echarge and Hscript and apply those factors later.
!  Do the following block only when lnoncollinear_EB, but not when lnoncollinear_EB_aver.
!
          if (lnoncollinear_EB) then
            boost=sqrt((e2-b2)**2+4.*edotb**2)
            gam_EB=sqrt21*sqrt(1.+(e2+b2)/boost)
            eprime=sqrt21*sqrt(e2-b2+boost)
            bprime=sqrt21*sqrt(b2-e2+boost)*sign(1.,edotb)
            if (lallow_bprime_zero) then
              where (eprime/=0.)
                where (bprime/=0.)
                  jprime1=1./(6.*pi**2)*eprime*abs(bprime)/tanh(pi*abs(bprime)/eprime)
                elsewhere
                  jprime1=1./(6.*pi**3)*eprime**2
                endwhere
                sigE1=abs(jprime1)*eprime/(gam_EB*boost)
                sigB1=abs(jprime1)*edotb/(eprime*gam_EB*boost)
              elsewhere
                sigE1=0.
                sigB1=0.
              endwhere
            else
              where (eprime/=0. .and. bprime/=0.)
                jprime1=1./(6.*pi**2)*eprime*abs(bprime)/tanh(pi*abs(bprime)/eprime)
                sigE1=abs(jprime1)*eprime/(gam_EB*boost)
                sigB1=abs(jprime1)*edotb/(eprime*gam_EB*boost)
              elsewhere
                sigE1=0.
                sigB1=0.
              endwhere
            endif
          endif
!
!  Repeat calculation of sigE and sigB. Do this first without
!  echarge and Hscript and apply those factors later.
!  Do the following block only when lnoncollinear_EB, but not when lnoncollinear_EB_aver.
!  Note: no multiplication by mass suppression factor is or can be done here.
!
          if (lcollinear_EB) then
            eprime=sqrt(e2)
            bprime=sqrt(b2)
            where (eprime/=0. .and. bprime/=0.)
              sigE1=1./(6.*pi**2)*bprime/tanh(pi*bprime/eprime)
              sigB1=0.
            elsewhere
              sigE1=0.
              sigB1=0.
            endwhere
          endif
        endif
      endif
!
!  Compute e2m per pencil. It becomes the total e2m after calling prep_ode_right
!  for each pencil. Also require either lrho_chi or lnoncollinear_EB
!  In case of lnoncollinear_EB_aver or lcollinear_EB_aver, do the computation outside prep_ode_right.
!
      if ((lmagnetic .and. lem_backreact) .and. (lrho_chi)) then
        if (lnoncollinear_EB .or. lnoncollinear_EB_aver .or. &
          lcollinear_EB .or. lcollinear_EB_aver) then
          e2m=e2m+sum(e2)
          b2m=b2m+sum(b2)
          if ((lnoncollinear_EB .or. lcollinear_EB)) then
            sigE1m=sigE1m+sum(sigE1)
            sigB1m=sigB1m+sum(sigB1)
          endif
        endif
      endif
!
    endsubroutine prep_ode_right
!********************************************************************
! Subroutines below needed only for GPUs, if you do not care about GPUs don't worry about them
!***********************************************************************
    subroutine read_sums_from_device
!
!  The sums needed for ODE and PDE advancement are computed on the GPUs in before_boundary.
!  Then we need them back on the host to advance the ODEs which this function does.
!
      use GPU, only: get_gpu_reduced_vars

      real, dimension(10) :: tmp

      call get_gpu_reduced_vars(tmp)

      a2rhom_all     = tmp(1)
      a2rhopm_all    = tmp(2)
      a2rhophim_all  = tmp(3)
      a2rhogphim_all = tmp(4)
      sigE1m_all_nonaver = tmp(5)
      sigB1m_all_nonaver = tmp(6)
      ddotam_all         = tmp(7)
      e2m_all            = tmp(8)
      b2m_all            = tmp(9)
      call get_Hscript_and_a2(Hscript,a2rhom_all)
      call get_echarge
      call get_sigE_and_B
      if (lfirst .and. lout) then
        a2rhom_all_diagnos     = a2rhom_all 
        a2rhopm_all_diagnos    = a2rhopm_all 
        a2rhophim_all_diagnos  = a2rhophim_all 
        a2rhogphim_all_diagnos = a2rhogphim_all
        ddotam_all_diagnos     = ddotam_all
        sigEm_all_diagnos      = sigEm_all
        sigBm_all_diagnos      = sigBm_all
      endif
!
    endsubroutine read_sums_from_device
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr
    use General , only: string_to_enum

    integer, parameter :: n_pars=50
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call string_to_enum(enum_hscript_choice,hscript_choice)
    call string_to_enum(enum_vprime_choice,vprime_choice)
    call copy_addr(enum_hscript_choice,p_par(1)) ! int
    call copy_addr(enum_vprime_choice,p_par(2)) ! int
    call copy_addr(iinfl_phi,p_par(3)) ! int
    call copy_addr(iinfl_dphi,p_par(4)) ! int
    call copy_addr(lnoncollinear_eb,p_par(5)) ! bool
    call copy_addr(lnoncollinear_eb_aver,p_par(6)) ! bool
    call copy_addr(lcollinear_eb,p_par(7)) ! bool
    call copy_addr(lcollinear_eb_aver,p_par(8)) ! bool
    call string_to_enum(enum_echarge_type,echarge_type)
    call copy_addr(enum_echarge_type,p_par(9)) ! int
    call copy_addr(echarge_const,p_par(10))
    call copy_addr(hscript0,p_par(11))
    call copy_addr(lzerohubble,p_par(12)) ! bool
    call copy_addr(axionmass2,p_par(13))
    call copy_addr(lambda_axion,p_par(14))
    call copy_addr(lconf_time,p_par(15)) ! bool
    call copy_addr(c_light_axion,p_par(16)) 
    call copy_addr(lem_backreact,p_par(17)) ! bool
    call copy_addr(ldt_backreact_infl,p_par(18)) ! bool
    call copy_addr(ndiv,p_par(19)) ! int
    call copy_addr(lrho_chi,p_par(20)) ! bool
    call copy_addr(iinfl_rho_chi,p_par(21)) ! int
    call copy_addr(cdt_rho_chi,p_par(22))
    call copy_addr(lflrw,p_par(23)) ! bool
    call copy_addr(iinfl_lna,p_par(24)) ! int
    call copy_addr(scale_rho_chi_heqn,p_par(25))
    call copy_addr(aphimax2,p_par(26))
    call copy_addr(gamma_phi0,p_par(27)) 
    call copy_addr(phi0,p_par(28))
    call copy_addr(lrho_rad,p_par(29)) ! bool

    call copy_addr(lrho_rad_apply,p_par(30)) ! bool
    call copy_addr(heating,p_par(31))
    call copy_addr(lsolve_for_phi,p_par(32)) ! bool
    call copy_addr(lheating,p_par(33)) ! bool

    call keep_compiler_quiet(rad_heating)
    call keep_compiler_quiet(rhophim_crit)
    call keep_compiler_quiet(phase_phi)
    call keep_compiler_quiet(lbackreact_infl)
    call keep_compiler_quiet(eps)
    call keep_compiler_quiet(kgaussian_dphi)
    call keep_compiler_quiet(ascale_heat_off)
    call keep_compiler_quiet(ascale_heat)

    call copy_addr(lrho_chi_inhom,p_par(34)) ! bool
    call copy_addr(lheating_always,p_par(35)) ! bool
    call copy_addr(lcombine_prep_ode_right_with_rhs,p_par(36)) ! bool
    call copy_addr(lsolve_for_phi_always,p_par(37)) ! bool

    call copy_addr(rhom_all,p_par(38))
    call copy_addr(axionmass,p_par(39))
    call copy_addr(cdt_phi,p_par(40))
    call copy_addr(gamma_phi_exp,p_par(41))
    call copy_addr(lsolve_for_phi2,p_par(42)) ! bool
    endsubroutine pushpars2c
!********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING        *************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!***********************************************************************
endmodule Special
