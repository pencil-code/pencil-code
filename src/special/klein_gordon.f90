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
! PENCILS PROVIDED phi; dphi; gphi(3); cov_der(4,4)
! PENCILS PROVIDED phi_doublet(3); dphi_doublet(3); phi_doublet_mod
! PENCILS EXPECTED GammaY, GammaW1, GammaW2, GammaW3
! PENCILS EXPECTED W1(3); W2(3); W3(3), aa(3)
! PENCILS PROVIDED psi; dpsi; gpsi(3)
!
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
  integer :: iphi=0, idphi=0, ilna=0, Ndiv=100, iinfl_rho_chi=0
  integer :: ipsi=0, idpsi=0
  integer :: iphi_up_re=0, iphi_up_im=0, iphi_down_re=0, iphi_down_im=0
  integer :: idphi_up_re=0, idphi_up_im=0, idphi_down_re=0, idphi_down_im=0
  real :: ncutoff_phi=1., phi_v=.1
  real :: phimass=1.06e-6, phimass2, ascale_ini=1.
  real :: psimass=1., psimass2
  real :: phi0=.44, dphi0=-1.69e-7, c_phi=1., lambda_phi=0., eps=.01
  real :: lambda_psi=0., coupl_phipsi=0., c_psi=1.
  real :: amplphi=.1, ampldphi=.0, kx_phi=1., ky_phi=0., kz_phi=0., phase_phi=0., width_phi=.1, offset=0.
  real :: amplpsi=0., ampldpsi=0.
  real :: initpower_phi=0.,  cutoff_phi=0.,  initpower2_phi=0.
  real :: initpower_dphi=0., cutoff_dphi=0., initpower2_dphi=0.
  real :: kgaussian_phi=0.,kpeak_phi=0., kpeak_dphi=0.
  real :: relhel_phi=0.
  real :: ddotam, a2rhopm, a2rhopm_all, a2rhom, a2rhom_all
  real :: edotbm, edotbm_all, e2m, e2m_all, b2m, b2m_all, a2rhophim, a2rhophim_all
  real :: sigE1m_all_nonaver, sigB1m_all_nonaver,sigEm_all,sigBm_all,sigEm_all_diagnos,sigBm_all_diagnos
  real :: a2rhogphim, a2rhogphim_all
  real :: a2rhopsim, a2rhopsim_all, a2rhogpsim, a2rhogpsim_all
  real :: a2, a21, Hscript
  real :: Hscript0=0., scale_rho_chi_Heqn=1., rho_chi_init=0., cdt_rho_chi=1.
  real :: amplee_BD_prefactor=0., deriv_prefactor_ee=-1.
  real :: echarge=.0, echarge_const=.303
  real :: count_eb0_all=0., eta_phi=0.
  real, pointer :: coupl_gw, coupl_gy
  ! logical, pointer :: llongitudinalE, llongitudinalW
  real, target :: ddotam_all
  real, pointer :: alpf, alpfpsi
  real, pointer :: sigE_prefactor, sigB_prefactor, mass_chi
  real, dimension (nx) :: dt1_special
  real, dimension (nx, 4, 3) :: dfdxs=0.
  logical :: lcompute_dphi0=.true., lem_backreact=.false.
  logical :: lscale_tobox=.true., ldt_klein_gordon=.true., lconf_time=.true.
  logical :: lskip_projection_phi=.false., lvectorpotential=.false., lflrw=.false.
  logical :: lrho_chi=.false., lno_noise_phi=.false., lno_noise_dphi=.false.
  logical, pointer :: lphi_hom, lpsi_hom
  logical, pointer :: lphi_linear_regime, lnoncollinear_EB, lnoncollinear_EB_aver
  logical, pointer :: lcollinear_EB, lcollinear_EB_aver, lmass_suppression
  logical, pointer :: lallow_bprime_zero
  logical :: lhiggs_friction=.false., lwaterfall=.false.
  real :: higgs_friction=0.
  logical :: lphi_doublet=.false., lphi_weakcharge=.false., lphi_hypercharge=.false.
  character (len=labellen) :: Vprime_choice='quadratic', Hscript_choice='set'
  character (len=labellen), dimension(ninit) :: initspecial='nothing'
  character (len=50) :: echarge_type='const', init_rho_chi='zero'
!
  namelist /special_init_pars/ &
      initspecial, phi0, dphi0, phimass, eps, ascale_ini, &
      lcompute_dphi0, lem_backreact, &
      c_phi, lambda_phi, amplphi, ampldphi, lno_noise_phi, lno_noise_dphi, &
      kx_phi, ky_phi, kz_phi, phase_phi, width_phi, offset, &
      initpower_phi, initpower2_phi, cutoff_phi, kgaussian_phi, kpeak_phi, &
      initpower_dphi, initpower2_dphi, cutoff_dphi, kpeak_dphi, &
      ncutoff_phi, lscale_tobox, Hscript0, Hscript_choice, phi_v, lflrw, &
      lrho_chi, scale_rho_chi_Heqn, amplee_BD_prefactor, deriv_prefactor_ee, &
      echarge_type, init_rho_chi, rho_chi_init, eta_phi, lphi_doublet, &
      lphi_weakcharge, lphi_hypercharge, lhiggs_friction, higgs_friction, &
      lwaterfall, lambda_psi, coupl_phipsi, c_psi, amplpsi, ampldpsi, psimass
!
  namelist /special_run_pars/ &
      initspecial, phi0, dphi0, phimass, eps, ascale_ini, &
      lem_backreact, c_phi, lambda_phi, Vprime_choice, &
      ldt_klein_gordon, Ndiv, Hscript0, Hscript_choice, &
      lflrw, lrho_chi, scale_rho_chi_Heqn, echarge_type, cdt_rho_chi, &
      phi_v, lhiggs_friction, higgs_friction, lwaterfall, lambda_psi, &
      coupl_phipsi, c_psi
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_phim=0       ! DIAG_DOC: $\left<\phi\right>$
  integer :: idiag_phi2m=0      ! DIAG_DOC: $\left<\phi^2\right>$
  integer :: idiag_phirms=0     ! DIAG_DOC: $\left<\phi^2\right>^{1/2}$
  integer :: idiag_dphim=0      ! DIAG_DOC: $\left<\phi'\right>$
  integer :: idiag_dphi2m=0     ! DIAG_DOC: $\left<(\phi')^2\right>$
  integer :: idiag_dphirms=0    ! DIAG_DOC: $\left<(\phi')^2\right>^{1/2}$
  integer :: idiag_psim=0       ! DIAG_DOC: $\left<\psi\right>$
  integer :: idiag_psi2m=0      ! DIAG_DOC: $\left<\psi^2\right>$
  integer :: idiag_psirms=0     ! DIAG_DOC: $\left<\psi^2\right>^{1/2}$
  integer :: idiag_dpsim=0      ! DIAG_DOC: $\left<\psi'\right>$
  integer :: idiag_dpsi2m=0     ! DIAG_DOC: $\left<(\psi')^2\right>$
  integer :: idiag_dpsirms=0    ! DIAG_DOC: $\left<(\psi')^2\right>^{1/2}$
  integer :: idiag_Hscriptm=0   ! DIAG_DOC: $\left<{\cal a*H}\right>$
  integer :: idiag_lnam=0       ! DIAG_DOC: $\left<\ln a\right>$
  integer :: idiag_ddotam=0     ! DIAG_DOC: $a''/a$
  integer :: idiag_a2rhopm=0    ! DIAG_DOC: $a^2 (rho+p)$
  integer :: idiag_a2rhom=0     ! DIAG_DOC: $a^2 rho$
  integer :: idiag_a2rhophim=0  ! DIAG_DOC: $a^2 rho$
  integer :: idiag_a2rhogphim=0 ! DIAG_DOC: $0.5 <grad phi^2>$
  integer :: idiag_a2rhopsim=0  ! DIAG_DOC: $a^2 rho$
  integer :: idiag_a2rhogpsim=0 ! DIAG_DOC: $0.5 <grad psi^2>$
  integer :: idiag_rho_chi=0    ! DIAG_DOC: $\rho_\chi$
  integer :: idiag_sigEma=0     ! DIAG_DOC: $\rho_\chi$
  integer :: idiag_sigBma=0     ! DIAG_DOC: $\rho_\chi$
  integer :: idiag_count_eb0a=0 ! DIAG_DOC: $f_\mathrm{EB0}$
!
  integer :: enum_hscript_choice = 0
  integer :: enum_vprime_choice = 0
  integer :: enum_echarge_type = 0

  real :: a2rhom_all_diagnos, a2rhopm_all_diagnos, a2rhophim_all_diagnos
  real :: a2rhogphim_all_diagnos, ddotam_all_diagnos
  real :: a2rhopsim_all_diagnos, a2rhogpsim_all_diagnos

  integer :: ia0 = 0, iW0 = 0
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
      lklein_gordon=.true.
!
      if (lphi_doublet) then
      ! alberto: register 4 components for the Higgs doublet in phi
        call farray_register_pde('phi_up_re',iphi_up_re)
        call farray_register_pde('phi_up_im',iphi_up_im)
        call farray_register_pde('phi_down_re',iphi_down_re)
        call farray_register_pde('phi_down_im',iphi_down_im)
        call farray_register_pde('dphi_up_re',idphi_up_re)
        call farray_register_pde('dphi_up_im',idphi_up_im)
        call farray_register_pde('dphi_down_re',idphi_down_re)
        call farray_register_pde('dphi_down_im',idphi_down_im)
        iphi=iphi_up_re
        idphi=idphi_up_re
        lwaterfall=.false.
      else
        call farray_register_pde('phi',iphi)
        call farray_register_pde('dphi',idphi)
        lphi_weakcharge=.false.
        lphi_hypercharge=.false.
        if (lwaterfall) then
          call farray_register_pde('psi',ipsi)
          call farray_register_pde('dpsi',idpsi)
        endif
      endif
!
      if (lflrw) call farray_register_ode('lna',ilna)
      if (lrho_chi) call farray_register_ode('infl_rho_chi',iinfl_rho_chi)
!
!  for power spectra, it is convenient to use ispecialvar and
!
      ispecialvar=iphi
      ispecialvar2=idphi
!
      call put_shared_variable('ddotam',ddotam_all,caller='register_klein_gordon')
      call put_shared_variable('Hscript',Hscript)
      call put_shared_variable('e2m_all',e2m_all)
      call put_shared_variable('b2m_all',b2m_all)
      call put_shared_variable('sigEm_all',sigEm_all,caller='register_klein_gordon')
      call put_shared_variable('sigBm_all',sigBm_all,caller='register_klein_gordon')
      call put_shared_variable('echarge',echarge,caller='register_klein_gordon')
      call put_shared_variable('lrho_chi',lrho_chi)
      call put_shared_variable('lphi_doublet',lphi_doublet)
      call put_shared_variable('lphi_weakcharge',lphi_weakcharge)
      call put_shared_variable('lphi_hypercharge',lphi_hypercharge)
      call put_shared_variable('lwaterfall',lwaterfall)
      call put_shared_variable('lflrw',lflrw)
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!
      use SharedVariables, only: get_shared_variable, put_shared_variable
      use FArrayManager, only: farray_index_by_name_ode, farray_index_by_name
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: iLCDM_lna
!
      if (lflrw) then
        iLCDM_lna=farray_index_by_name_ode('iLCDM_lna')
        if (iLCDM_lna>0) call fatal_error('initialize_special', 'there is a conflict with iLCDM_lna')
      endif
!
!  set phimass**2
!
      phimass2=phimass**2
      if (lwaterfall) psimass2=psimass**2
!
      if (lmagnetic .and. lem_backreact) then
        call get_shared_variable('alpf',alpf,caller='initialize_klein_gordon')
        call get_shared_variable('alpfpsi',alpfpsi,caller='initialize_klein_gordon')
        call get_shared_variable('lphi_linear_regime',lphi_linear_regime)
        call get_shared_variable('sigE_prefactor',sigE_prefactor)
        call get_shared_variable('sigB_prefactor',sigB_prefactor)
        call get_shared_variable('lcollinear_EB',lcollinear_EB)
        call get_shared_variable('lcollinear_EB_aver',lcollinear_EB_aver)
        call get_shared_variable('lnoncollinear_EB',lnoncollinear_EB)
        call get_shared_variable('lnoncollinear_EB_aver',lnoncollinear_EB_aver)
        call get_shared_variable('lmass_suppression',lmass_suppression)
        call get_shared_variable('lphi_hom',lphi_hom)
        call get_shared_variable('lpsi_hom',lpsi_hom)
        call get_shared_variable('lallow_bprime_zero',lallow_bprime_zero)
        call get_shared_variable('mass_chi',mass_chi)
      else
        if (.not.associated(alpf)) allocate(alpf,lphi_linear_regime, &
          sigE_prefactor, sigB_prefactor, lcollinear_EB, lcollinear_EB_aver, &
          lnoncollinear_EB, lnoncollinear_EB_aver, lmass_suppression, &
          lphi_hom, lpsi_hom, lallow_bprime_zero, mass_chi, alpfpsi)
        alpf=0.
        alpfpsi=0.
        lphi_linear_regime=.false.
        sigE_prefactor=0.
        sigB_prefactor=0.
        lcollinear_EB=.false.
        lcollinear_EB_aver=.false.
        lnoncollinear_EB=.false.
        lnoncollinear_EB_aver=.false.
        lmass_suppression=.false.
        lphi_hom=.false.
        lpsi_hom=.false.
        lallow_bprime_zero=.false.
        mass_chi=0.
      endif

      ! if iee = 0 then disp_current module is not called
      if (iee /= 0 .and. lphi_hypercharge) then
        call get_shared_variable('coupl_gy',coupl_gy)
        ! call get_shared_variable('llongitudinalE',llongitudinalE)
      else
        if (.not.associated(coupl_gy)) then
          ! allocate(coupl_gy, llongitudinalE)
          allocate(coupl_gy)
          coupl_gy=0.
          ! llongitudinalE=.false.
        endif
        lphi_hypercharge=.false.
      endif
      !  if iWW = 0 then electroweaksu2 module is not called
      if (iWW /= 0 .and. lphi_weakcharge) then
        call get_shared_variable('coupl_gw',coupl_gw)
        ! call get_shared_variable('llongitudinalW',llongitudinalW)
      else
        if (.not.associated(coupl_gw)) then
          ! allocate(coupl_gw, llongitudinalW)
          allocate(coupl_gw)
          coupl_gw=0.
          ! llongitudinalW=.false.
        endif
        lphi_weakcharge=.false.
      endif
!
      call keep_compiler_quiet(f)
      if (lphi_hypercharge) then
        ia0 = farray_index_by_name('a0')
      endif
      if (lphi_weakcharge) then
        iW0 = farray_index_by_name('W0')
      endif
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
      real :: Vpotential, Hubble_ini, phi_gam, amplphi_BD, amplee_BD, deriv_prefactor
      integer :: j
      real :: lnascale
!
      intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      do j=1,ninit
        select case (initspecial(j))
          case ('nothing'); if (lroot) print*,'init_special: nothing'
          case ('constant')
            f(:,:,:,iphi)=f(:,:,:,iphi)+amplphi
            f(:,:,:,idphi)=f(:,:,:,idphi)+ampldphi
            if (lwaterfall) then
              f(:,:,:,ipsi)=f(:,:,:,ipsi)+amplpsi
              f(:,:,:,idpsi)=f(:,:,:,idpsi)+ampldpsi
            endif
          case ('phi=sinkx')
            f(:,:,:,iphi)=f(:,:,:,iphi) &
              +spread(spread(amplphi*sin(kx_phi*x),2,my),3,mz)
          case ('phi=tanhkx')
            f(:,:,:,iphi)=f(:,:,:,iphi) &
              +spread(spread(.5*amplphi*(1.+tanh(kx_phi*(x-offset))),2,my),3,mz)
          ! sine-Gordon solution
          case ('phi=atan_exp_kx')
            phi_gam=1./sqrt(1.-phi_v**2)
            f(:,:,:,iphi)=f(:,:,:,iphi) &
              +spread(spread(4.*amplphi*atan(exp(phi_gam*kx_phi*(x-offset))),2,my),3,mz)
            f(:,:,:,idphi)=f(:,:,:,idphi)+spread(spread( &
              -4.*amplphi*kx_phi*phi_gam*phi_v*exp(phi_gam*kx_phi*(x-offset)) &
              /(exp(2.*phi_gam*kx_phi*(x-offset))+1.) &
              ,2,my),3,mz)
          case ('nophi')
            Vpotential=.5*phimass2*phi0**2
            dphi0=0.
            tstart=-sqrt(3./(8.*pi))/(ascale_ini*sqrt(Vpotential))
            t=tstart
            Hubble_ini=sqrt(8.*pi/3.*(.5*dphi0**2+.5*phimass2*phi0**2*ascale_ini**2))
            lnascale=log(ascale_ini)
            if (lroot .and. lflrw) f_ode(ilna)=lnascale
  !
          case ('default')
            Vpotential=.5*phimass2*phi0**2
            Hubble_ini=sqrt(8.*pi/3.*(.5*phimass2*phi0**2*ascale_ini**2))
            if (lcompute_dphi0) dphi0=-sqrt(1/(12.*pi))*phimass*ascale_ini
            tstart=-1/(ascale_ini*Hubble_ini)
            t=tstart
            lnascale=log(ascale_ini)
            f(:,:,:,iphi)   = f(:,:,:,iphi)   + phi0
            f(:,:,:,idphi)  = f(:,:,:,idphi)  + dphi0
            if (lroot .and. lflrw) then
              f_ode(ilna)   =lnascale
              a2                 =exp(f_ode(ilna))**2
              Hscript            =Hubble_ini/exp(lnascale)
            endif
          case ('gaussian-noise')
            call gaunoise(amplphi,f,iphi)
          case ('sinwave-phase')
            call hat(amplphi,f,iphi,width_phi,kx_phi,ky_phi,kz_phi)
            f(:,:,:,iphi)=f(:,:,:,iphi)+offset
          case ('phi_power_randomphase')
            call power_randomphase_hel(amplphi,initpower_phi,initpower2_phi, &
              cutoff_phi, ncutoff_phi, kpeak_phi, f, iphi, iphi, &
              relhel_phi, kgaussian_phi, lskip_projection_phi, lvectorpotential, &
              lscale_tobox, lpower_profile_file=.false., lno_noise=lno_noise_phi)
          case ('dphi_power_randomphase')
            call power_randomphase_hel(ampldphi,initpower_dphi,initpower2_dphi, &
              cutoff_dphi,ncutoff_phi,kpeak_dphi,f,idphi,idphi, &
              relhel_phi,kgaussian_phi, lskip_projection_phi, lvectorpotential, &
              lscale_tobox, lpower_profile_file=.false., lno_noise=lno_noise_dphi)
  !
  !  For Bunch-Davies, the amplitude Hubble_ini is used.
  !  We apply this optionally here also to the gauge field.
  !
          case ('Bunch-Davies')
            if (lroot) print*,'Hubble_ini=',Hubble_ini
            ! alberto: where is Hubble_ini defined for this initial condition?
            amplphi_BD=amplphi*Hubble_ini
            deriv_prefactor=1.
            call bunch_davies(f,iphi,iphi,idphi,idphi, &
                              amplphi_BD,kpeak_phi,deriv_prefactor)
            if (amplee_BD_prefactor/=0.) then
              deriv_prefactor=deriv_prefactor_ee
              amplee_BD=amplee_BD_prefactor*Hubble_ini
              call bunch_davies(f,iax,iaz,iex,iez,amplee_BD,kpeak_phi,deriv_prefactor)
            endif
          case ('phi_doublet')
            if (.not.lphi_doublet) &
                call fatal_error("init_special: lphi_doublet=.false. but initspecial='phi_doublet'", &
                                 trim(initspecial(j)))
            ! alberto: example initialization of a Higgs doublet with a sine wave in the
            ! real part of the upper component and the same for the lower component
            f(:,:,:,iphi_up_re)=f(:,:,:,iphi_up_re) &
              +spread(spread(amplphi*sin(kx_phi*x),2,my),3,mz)
            f(:,:,:,iphi_down_re)=f(:,:,:,iphi_down_re) &
              +spread(spread(amplphi*sin(kx_phi*x),2,my),3,mz)
          case default
            call fatal_error("init_special: No such initspecial: ", trim(initspecial(j)))
        endselect
      enddo
!
!  initial condition for energy density of charged particles
!
      if (lroot .and. lrho_chi) then
        select case (init_rho_chi)
          case ('zero');  f_ode(iinfl_rho_chi)=0.
          case ('given'); f_ode(iinfl_rho_chi)=rho_chi_init
          case default
            call fatal_error("init_special: No such init_rho_chi: ", trim(init_rho_chi))
        endselect
      endif
!
      call mpibcast_real(a2)
      call mpibcast_real(Hscript)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special

    integer :: i
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
!
!  Call pencils phi and dphi
      lpenc_requested(i_phi)=.true.
      lpenc_requested(i_dphi)=.true.

      if (lphi_doublet) then
        lpenc_requested(i_phi_doublet)=.true.
        lpenc_requested(i_dphi_doublet)=.true.
        lpenc_requested(i_phi_doublet_mod)=.true.
        if (lphi_hypercharge .or. lphi_weakcharge) then
          lpenc_requested(i_cov_der)=.true.
          if (lphi_hypercharge) then
            lpenc_requested(i_GammaY)=.true.
            lpenc_requested(i_aa)=.true.
          endif
          if (lphi_weakcharge) then
            lpenc_requested(i_GammaW1)=.true.
            lpenc_requested(i_GammaW2)=.true.
            lpenc_requested(i_GammaW3)=.true.
            lpenc_requested(i_W1)=.true.
            lpenc_requested(i_W2)=.true.
            lpenc_requested(i_W3)=.true.
          endif
        endif
      endif

      if (lwaterfall) then
        lpenc_requested(i_psi)=.true.
        lpenc_requested(i_dpsi)=.true.
      endif
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
      use Sub, only: grad, div
      use Deriv, only: der
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
      integer ::  i, j

! phi
      if (lpencil(i_phi)) p%phi = f(l1:l2,m,n,iphi)
! dphi
      if (lpencil(i_dphi)) p%dphi=f(l1:l2,m,n,idphi)
! psi
      if (lpencil(i_psi)) p%psi = f(l1:l2,m,n,ipsi)
! dpsi
      if (lpencil(i_dpsi)) p%dpsi=f(l1:l2,m,n,idpsi)
! phi_doublet (only computes the 3 remaining components, the first is in p%phi)
      if (lpencil(i_phi_doublet_mod)) then
        do i=1,3
          p%phi_doublet(:,i)=f(l1:l2,m,n,iphi+i)
          p%dphi_doublet(:,i)=f(l1:l2,m,n,idphi+i)
        enddo
        p%phi_doublet_mod=sqrt(f(l1:l2,m,n,iphi_up_re)**2 + &
                               f(l1:l2,m,n,iphi_up_im)**2 + &
                               f(l1:l2,m,n,iphi_down_re)**2 + &
                               f(l1:l2,m,n,iphi_down_im)**2)
      endif
! gphi
      if (lpencil(i_gphi)) call grad(f,iphi,p%gphi)
! gpsi
      if (lpencil(i_gpsi)) call grad(f,ipsi,p%gpsi)
!
! cov_der computation for Higgs doublet with weak and/or hypercharge
! covariant derivative is D_mu = partial_mu - i g W_mu^a tau^a/2 - i g' Y B_mu/2
! where tau^a are the Pauli matrices and Y is the hypercharge
!
      if (lpencil(i_cov_der)) then
        do i=0,3
          p%cov_der(:, 1, i+1)=f(l1:l2,m,n,idphi+i)
          p%cov_der(:, 1, i+1)=f(l1:l2,m,n,idphi+i)
          if (.not. lphi_hom) then
            do j=1,3
              call der(f, iphi+i, dfdxs(:, i+1, j), j)
              p%cov_der(:, j+1, i+1) = dfdxs(:, i, j)
            enddo
          endif
        enddo
        ! when lphi_hypercharge is chosen and disp_current.f90 is used
        ! add terms for hypercharge to covariant derivative
        if (lphi_hypercharge) then
          if (ia0 > 0) then
            ! U(1) term (first block), proportional to a0 (Y0)
            p%cov_der(:,1,1) = p%cov_der(:,1,1) + 0.5*coupl_gy*f(l1:l2,m,n,ia0)*f(l1:l2,m,n,iphi_up_im)
            p%cov_der(:,1,2) = p%cov_der(:,1,2) - 0.5*coupl_gy*f(l1:l2,m,n,ia0)*f(l1:l2,m,n,iphi_up_re)
            p%cov_der(:,1,3) = p%cov_der(:,1,3) + 0.5*coupl_gy*f(l1:l2,m,n,ia0)*f(l1:l2,m,n,iphi_down_im)
            p%cov_der(:,1,4) = p%cov_der(:,1,4) - 0.5*coupl_gy*f(l1:l2,m,n,ia0)*f(l1:l2,m,n,iphi_down_re)
          endif
          ! Remaining blocks (no Y0, multiply two field components)
          do i = 1, 3
            p%cov_der(:,i+1,1) = p%cov_der(:,i+1,1) + 0.5*coupl_gy*f(l1:l2,m,n,iaa+i-1)*f(l1:l2,m,n,iphi_up_im)
            p%cov_der(:,i+1,2) = p%cov_der(:,i+1,2) - 0.5*coupl_gy*f(l1:l2,m,n,iaa+i-1)*f(l1:l2,m,n,iphi_up_re)
            p%cov_der(:,i+1,3) = p%cov_der(:,i+1,3) + 0.5*coupl_gy*f(l1:l2,m,n,iaa+i-1)*f(l1:l2,m,n,iphi_down_im)
            p%cov_der(:,i+1,4) = p%cov_der(:,i+1,4) - 0.5*coupl_gy*f(l1:l2,m,n,iaa+i-1)*f(l1:l2,m,n,iphi_down_re)
          enddo
        endif

        ! when lphi_hypercharge is chosen and electroweak_su2.f90 is used
        ! add terms for weakcharge to covariant derivative
        if (lphi_weakcharge) then
          if (iW0 > 0) then
            p%cov_der(:,1,1) = p%cov_der(:,1,1) + &
                0.5*coupl_gw*(f(l1:l2,m,n,iW0)*f(l1:l2,m,n,iphi_down_im) - &
                f(l1:l2,m,n,iW0+1)*f(l1:l2,m,n,iphi_down_re) + f(l1:l2,m,n,iW0+2)*f(l1:l2,m,n,iphi_up_im))

            p%cov_der(:,1,2) = p%cov_der(:,1,2) - &
                  0.5*coupl_gw*(f(l1:l2,m,n,iW0)*f(l1:l2,m,n,iphi_down_re) + &
                  f(l1:l2,m,n,iW0+1)*f(l1:l2,m,n,iphi_down_im) + f(l1:l2,m,n,iW0+2)*f(l1:l2,m,n,iphi_up_re))

            p%cov_der(:,1,3) = p%cov_der(:,1,3) + &
                  0.5*coupl_gw*(f(l1:l2,m,n,iW0)*f(l1:l2,m,n,iphi_up_im) + &
                  f(l1:l2,m,n,iW0+1)*f(l1:l2,m,n,iphi_up_re) - f(l1:l2,m,n,iW0+2)*f(l1:l2,m,n,iphi_down_im))

            p%cov_der(:,1,4) = p%cov_der(:,1,4) - &
                  0.5*coupl_gw*(f(l1:l2,m,n,iW0)*f(l1:l2,m,n,iphi_up_re) - &
                  f(l1:l2,m,n,iW0+1)*f(l1:l2,m,n,iphi_up_im) - f(l1:l2,m,n,iW0+2)*f(l1:l2,m,n,iphi_down_re))
          endif

          ! iWW:iWW+2 -> W1x:W1z
          ! iWW+3:iWW+5 -> W2x:W2z
          ! iWW+6:iWW+8 -> W3x:W3z
          p%cov_der(:,2,1) = p%cov_der(:,2,1) + &
                0.5*coupl_gw*(f(l1:l2,m,n,iWW1)*f(l1:l2,m,n,iphi_down_im) - &
                f(l1:l2,m,n,iWW2)*f(l1:l2,m,n,iphi_down_re) + &
                f(l1:l2,m,n,iWW3)*f(l1:l2,m,n,iphi_up_im))

          p%cov_der(:,2,2) = p%cov_der(:,2,2) - &
                0.5*coupl_gw*(f(l1:l2,m,n,iWW1)*f(l1:l2,m,n,iphi_down_re) + &
                f(l1:l2,m,n,iWW2)*f(l1:l2,m,n,iphi_down_im) + &
                f(l1:l2,m,n,iWW3)*f(l1:l2,m,n,iphi_up_re))

          p%cov_der(:,2,3) = p%cov_der(:,2,3) + &
                0.5*coupl_gw*(f(l1:l2,m,n,iWW1)*f(l1:l2,m,n,iphi_up_im) + &
                f(l1:l2,m,n,iWW2)*f(l1:l2,m,n,iphi_up_re) - &
                f(l1:l2,m,n,iWW3)*f(l1:l2,m,n,iphi_down_im))

          p%cov_der(:,2,4) = p%cov_der(:,2,4) - &
                0.5*coupl_gw*(f(l1:l2,m,n,iWW1)*f(l1:l2,m,n,iphi_up_re) - &
                f(l1:l2,m,n,iWW2)*f(l1:l2,m,n,iphi_up_im) - &
                f(l1:l2,m,n,iWW3)*f(l1:l2,m,n,iphi_down_re))

          p%cov_der(:,3,1) = p%cov_der(:,3,1) + &
                0.5*coupl_gw*(f(l1:l2,m,n,iWW1+1)*f(l1:l2,m,n,iphi_down_im) - &
                f(l1:l2,m,n,iWW2+1)*f(l1:l2,m,n,iphi_down_re) + &
                f(l1:l2,m,n,iWW3+1)*f(l1:l2,m,n,iphi_up_im))

          p%cov_der(:,3,2) = p%cov_der(:,3,2) - &
                0.5*coupl_gw*(f(l1:l2,m,n,iWW1+1)*f(l1:l2,m,n,iphi_down_re) + &
                f(l1:l2,m,n,iWW2+1)*f(l1:l2,m,n,iphi_down_im) + &
                f(l1:l2,m,n,iWW3+1)*f(l1:l2,m,n,iphi_up_re))

          p%cov_der(:,3,3) = p%cov_der(:,3,3) + &
                0.5*coupl_gw*(f(l1:l2,m,n,iWW1+1)*f(l1:l2,m,n,iphi_up_im) + &
                f(l1:l2,m,n,iWW2+1)*f(l1:l2,m,n,iphi_up_re) - &
                f(l1:l2,m,n,iWW3+1)*f(l1:l2,m,n,iphi_down_im))

          p%cov_der(:,3,4) = p%cov_der(:,3,4) - &
                0.5*coupl_gw*(f(l1:l2,m,n,iWW1+1)*f(l1:l2,m,n,iphi_up_re) - &
                f(l1:l2,m,n,iWW2+1)*f(l1:l2,m,n,iphi_up_im) - &
                f(l1:l2,m,n,iWW3+1)*f(l1:l2,m,n,iphi_down_re))

          p%cov_der(:,4,1) = p%cov_der(:,4,1) + &
                0.5*coupl_gw*(f(l1:l2,m,n,iWW1+2)*f(l1:l2,m,n,iphi_down_im) - &
                f(l1:l2,m,n,iWW2+2)*f(l1:l2,m,n,iphi_down_re) + &
                f(l1:l2,m,n,iWW3+2)*f(l1:l2,m,n,iphi_up_im))

          p%cov_der(:,4,2) = p%cov_der(:,4,2) - &
                0.5*coupl_gw*(f(l1:l2,m,n,iWW1+2)*f(l1:l2,m,n,iphi_down_re) + &
                f(l1:l2,m,n,iWW2+2)*f(l1:l2,m,n,iphi_down_im) + &
                f(l1:l2,m,n,iWW3+2)*f(l1:l2,m,n,iphi_up_re))

          p%cov_der(:,4,3) = p%cov_der(:,4,3) + &
                0.5*coupl_gw*(f(l1:l2,m,n,iWW1+2)*f(l1:l2,m,n,iphi_up_im) + &
                f(l1:l2,m,n,iWW2+2)*f(l1:l2,m,n,iphi_up_re) - &
                f(l1:l2,m,n,iWW3+2)*f(l1:l2,m,n,iphi_down_im))

          p%cov_der(:,4,4) = p%cov_der(:,4,4) - &
                0.5*coupl_gw*(f(l1:l2,m,n,iWW1+2)*f(l1:l2,m,n,iphi_up_re) - &
                f(l1:l2,m,n,iWW2+2)*f(l1:l2,m,n,iphi_up_im) - &
                f(l1:l2,m,n,iWW3+2)*f(l1:l2,m,n,iphi_down_re))
        endif
        ! p%cov_der = cov_der
      endif
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine get_Hscript_and_a2(Hscript,a2rhom_all)

      real, intent(out) :: Hscript
      real, intent(in), optional  :: a2rhom_all
!
!  Choice of prescription for Hscript
!  alberto: changed default to 'set' with Hscript0=0, removed lzeroHubble
!           as it trivially corresponds to new default choice
!
      select case (Hscript_choice)
        case ('set')
          Hscript=Hscript0
          a2=1.
          a21=1./a2
        case ('friedmann')
          Hscript=sqrt((8.*pi/3.)*a2rhom_all)
          if (lgpu) call get_a2
        case default
          call fatal_error("dspecial_dt: No such Hscript_choice: ", trim(Hscript_choice))
      endselect

    endsubroutine get_Hscript_and_a2
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
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
!   4-sep-25/alberto: adapted from backreact_infl
!   6-sep-25/alberto: added Higgs doublet case
!   14-sep-25/alberto: added second scalar field psi for waterfall potential
!
      use Diagnostics, only: sum_mn_name, max_mn_name, save_name
      use Sub, only: dot_mn, del2, div
      use Deriv, only: der
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: Vprime, Vprimepsi, Vprime_aux, total_fric
      real, dimension (nx, 4) :: del2phi_doublet=0.
      real, dimension (nx) :: tmp, del2phi, del2psi
      real :: pref_Vprime=1., pref_Hubble=2., pref_del2=1., pref_alpf
      type (pencil_case) :: p
      integer :: i
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!  Choice of different potentials.
!  For the 1-cos profile, -Vprime (on the rhs) enters with -sin().
!
      select case (Vprime_choice)
        case ('quadratic'); Vprime=phimass2*p%phi
        case ('quartic'); Vprime=phimass2*p%phi+(lambda_phi/6.)*p%phi**3
        case ('cos-profile'); Vprime=phimass2*lambda_phi*sin(lambda_phi*p%phi)
        ! for doublet case, Vprime = (dV/d|Phi|)/|Phi|
        case ('doublet')
          if (.not.lphi_doublet) &
              call fatal_error("dspecial_dt: lphi_doublet=.false. but Vprime_choice='doublet'", &
                               trim(Vprime_choice))
          Vprime=2*lambda_phi*(p%phi_doublet_mod**2 - eta_phi**2)*p%phi_doublet_mod
        case ('waterfall')
          if (.not.lwaterfall) &
              call fatal_error("dspecial_dt: lwaterfall=.false. but Vprime_choice='waterfall'", &
                               trim(Vprime_choice))
          Vprime=phimass2*p%phi+coupl_phipsi**2*p%phi*p%psi**2
          Vprimepsi=lambda_psi*p%psi**3-psimass2*p%psi+coupl_phipsi**2*p%psi*p%phi**2
        case default
          call fatal_error("dspecial_dt: No such Vprime_choice: ", trim(Vprime_choice))
      endselect
!
!  Update df.
!  dphi/dt = psi
!  dpsi/dt = - ...
!
! alberto: determine prefactors for the different terms beforehand
!
      if (lconf_time) then
        pref_Vprime=a2; pref_alpf=a21
      ! alberto: for cosmic time, should coefficient of Hscript be 3?
      else
        pref_Hubble=3.; pref_Vprime=1.; pref_del2=a21
        pref_alpf=a21**2
      endif
!
! alberto: right-hand-side for Klein-Gordon equation with Higgs doublet
!           in presence of U(1) and/or SU(2) gauge fields
!
      if (lphi_doublet) then

        if (lhiggs_friction) then
          total_fric=p%phi*p%dphi
          do i=1,3
            total_fric=total_fric+p%phi_doublet(:,i)*p%dphi_doublet(:,i)
          enddo
          total_fric=total_fric/max(p%phi_doublet_mod**2, 1e-30)
        endif
        do i=0,3
          ! dphi/dt = dphi
          df(l1:l2,m,n,iphi+i)=df(l1:l2,m,n,iphi+i)+f(l1:l2,m,n,idphi+i)
          ! laplacian of the 4 components of the Higgs doublet
          if (c_phi/=0. .and. .not. lphi_hom) then
            call del2(f, iphi+i, del2phi_doublet(:, i+1))
          endif
          ! alberto: added Higgs friction (see eq. 3 in 1902.02751)
          if (lhiggs_friction) then
            df(l1:l2,m,n,iphi+i)=df(l1:l2,m,n,iphi+i)- & 
                higgs_friction*f(l1:l2,m,n,iphi+i)*total_fric
          endif
        enddo

        if (c_phi/=0. .and. lphi_hypercharge) then
          ! terms of the covariant Laplacian from U(1) gauge fields
          ! note that p%phi = f(l1:l2,m,n,iphi) = phi_up_re,
          ! p%phi_doublet(:,1) = phi_up_im,
          ! p%phi_doublet(:,2) = phi_down_re,
          ! p%phi_doublet(:,3) = phi_down_im
          ! del2phi_up_re
          del2phi_doublet(:,1) = del2phi_doublet(:,1) - &
            0.5*coupl_gy*(-p%aa(:,1)*dfdxs(:,2,1) - p%aa(:,2)*dfdxs(:,2,2) - &
            p%aa(:,3)*dfdxs(:,2,3) - p%aa(:,1)*p%cov_der(:,2,2) - &
            p%aa(:,2)*p%cov_der(:,3,2) - p%aa(:,3)*p%cov_der(:,4,2) - &
            p%GammaY*p%phi_doublet(:,1))

          ! del2phi_up_im
          del2phi_doublet(:,2) = del2phi_doublet(:,2) + &
            0.5*coupl_gy*(-p%aa(:,1)*dfdxs(:,1,1) - p%aa(:,2)*dfdxs(:,1,2) - &
            p%aa(:,3)*dfdxs(:,1,3) - p%aa(:,1)*p%cov_der(:,2,1) - &
            p%aa(:,2)*p%cov_der(:,3,1) - p%aa(:,3)*p%cov_der(:,4,1) - &
            p%GammaY*p%phi)

          ! del2phi_down_re
          del2phi_doublet(:,3) = del2phi_doublet(:,3) - &
            0.5*coupl_gy*(-p%aa(:,1)*dfdxs(:,4,1) - p%aa(:,2)*dfdxs(:,4,2) - &
            p%aa(:,3)*dfdxs(:,4,3) - p%aa(:,1)*p%cov_der(:,2,4) - &
            p%aa(:,2)*p%cov_der(:,3,4) - p%aa(:,3)*p%cov_der(:,4,4) - &
            p%GammaY*p%phi_doublet(:,3))

          ! del2phi_down_im
          del2phi_doublet(:,4) = del2phi_doublet(:,4) + &
            0.5*coupl_gy*(-p%aa(:,1)*dfdxs(:,3,1) - p%aa(:,2)*dfdxs(:,3,2) - &
            p%aa(:,3)*dfdxs(:,3,3) - p%aa(:,1)*p%cov_der(:,2,3) - &
            p%aa(:,2)*p%cov_der(:,3,3) - p%aa(:,3)*p%cov_der(:,4,3) - &
            p%GammaY*p%phi_doublet(:,2))
        endif
        if (c_phi/=0. .and. lphi_weakcharge) then
          ! terms of the covariant Laplacian from SU(2) gauge fields
          ! del2phi_up_re
          del2phi_doublet(:,1) = del2phi_doublet(:,1) - &
            0.5*coupl_gw*(-p%W1(:,1)*dfdxs(:,4,1) - &
            p%W1(:,2)*dfdxs(:,4,2) - p%W1(:,3)*dfdxs(:,4,3) + &
            p%W2(:,1)*dfdxs(:,3,1) + p%W2(:,2)*dfdxs(:,3,2) + &
            p%W2(:,3)*dfdxs(:,3,3) - p%W3(:,1)*dfdxs(:,2,1) + &
            p%W3(:,2)*dfdxs(:,2,2) + p%W3(:,3)*dfdxs(:,2,3) - &
            p%W1(:,1)*p%cov_der(:,2,4) - p%W1(:,2)*p%cov_der(:,3,4) - &
            p%W1(:,3)*p%cov_der(:,4,4) + p%W2(:,1)*p%cov_der(:,2,3) + &
            p%W2(:,2)*p%cov_der(:,3,3) + p%W2(:,3)*p%cov_der(:,4,3) - &
            p%W3(:,1)*p%cov_der(:,2,2) - p%W3(:,2)*p%cov_der(:,3,2) - &
            p%W3(:,3)*p%cov_der(:,4,2) - p%GammaW3*p%phi_doublet(:,1) + &
            p%GammaW2*p%phi_doublet(:,2) - p%GammaW1*p%phi_doublet(:,3))

          ! del2phi_up_im
          del2phi_doublet(:,2) = del2phi_doublet(:,2) + &
            0.5*coupl_gw*(-p%W1(:,1)*dfdxs(:,3,1) - &
            p%W1(:,2)*dfdxs(:,3,2) - p%W1(:,3)*dfdxs(:,3,3) - &
            p%W2(:,1)*dfdxs(:,4,1) - p%W2(:,2)*dfdxs(:,4,2) - &
            p%W2(:,3)*dfdxs(:,4,3) - p%W3(:,1)*dfdxs(:,1,1) - &
            p%W3(:,2)*dfdxs(:,1,2) - p%W3(:,3)*dfdxs(:,1,3) - &
            p%W1(:,1)*p%cov_der(:,2,3) - p%W1(:,2)*p%cov_der(:,3,3) - &
            p%W1(:,3)*p%cov_der(:,4,3) - p%W2(:,1)*p%cov_der(:,2,4) - &
            p%W2(:,2)*p%cov_der(:,3,4) - p%W2(:,3)*p%cov_der(:,4,4) - &
            p%W3(:,1)*p%cov_der(:,2,1) - p%W3(:,2)*p%cov_der(:,3,1) - &
            p%W3(:,3)*p%cov_der(:,4,1) - p%gammaW3*p%phi - &
            p%gammaW1*p%phi_doublet(:,2) - p%gammaW2*p%phi_doublet(:,3))

          ! del2phi_down_re
          del2phi_doublet(:,3) = del2phi_doublet(:,3) - &
            0.5*coupl_gw*(-p%W1(:,1)*dfdxs(:,2,1) - &
            p%W1(:,2)*dfdxs(:,2,2) - p%W1(:,3)*dfdxs(:,2,3) - &
            p%W2(:,1)*dfdxs(:,1,1) - p%W2(:,2)*dfdxs(:,1,2) - &
            p%W2(:,3)*dfdxs(:,1,3) + p%W3(:,1)*dfdxs(:,4,1) + &
            p%W3(:,2)*dfdxs(:,4,2) + p%W3(:,3)*dfdxs(:,4,3) - &
            p%W1(:,1)*p%cov_der(:,2,2) - p%W1(:,2)*p%cov_der(:,3,2) - &
            p%W1(:,3)*p%cov_der(:,4,2) - p%W2(:,1)*p%cov_der(:,2,1) - &
            p%W2(:,2)*p%cov_der(:,3,1) - p%W2(:,3)*p%cov_der(:,4,1) + &
            p%W3(:,1)*p%cov_der(:,2,4) + p%W3(:,2)*p%cov_der(:,3,4) + &
            p%W3(:,3)*p%cov_der(:,4,4) + p%gammaW3*p%phi_doublet(:,2) - &
            p%gammaW2*p%phi - p%gammaW1*p%phi_doublet(:,1))

          ! del2phi_down_im
          del2phi_doublet(:,4) = del2phi_doublet(:,4) + &
            0.5*coupl_gw*(-p%W1(:,1)*dfdxs(:,1,1) - &
            p%W1(:,2)*dfdxs(:,1,2) - p%W1(:,3)*dfdxs(:,1,3) + &
            p%W2(:,1)*dfdxs(:,2,1) + p%W2(:,2)*dfdxs(:,2,2) + &
            p%W2(:,3)*dfdxs(:,2,3) + p%W3(:,1)*dfdxs(:,3,1) + &
            p%W3(:,2)*dfdxs(:,3,2) + p%W3(:,3)*dfdxs(:,3,3) - &
            p%W1(:,1)*p%cov_der(:,2,1) - p%W1(:,2)*p%cov_der(:,3,1) - &
            p%W1(:,3)*p%cov_der(:,4,1) + p%W2(:,1)*p%cov_der(:,2,2) + &
            p%W2(:,2)*p%cov_der(:,3,2) + p%W2(:,3)*p%cov_der(:,4,2) + &
            p%W3(:,1)*p%cov_der(:,2,3) + p%W3(:,2)*p%cov_der(:,3,3) + &
            p%W3(:,3)*p%cov_der(:,4,3) + p%gammaW3*p%phi_doublet(:,2) - &
            p%gammaW1*p%phi + p%gammaW2*p%phi_doublet(:,1))

        endif
        do i=0,3
          Vprime_aux=0.
          if (Vprime_choice=='doublet') then
            Vprime_aux=Vprime*f(l1:l2,m,n,iphi+i)
          else
            if (i == 0) Vprime_aux=Vprime
          endif
          df(l1:l2,m,n,idphi+i)=df(l1:l2,m,n,idphi+i) - &
            pref_Hubble*Hscript*f(l1:l2,m,n,idphi+i) - &
            pref_Vprime*Vprime_aux
          if (c_phi/=0) &
            df(l1:l2,m,n,idphi+i)=df(l1:l2,m,n,idphi+i) + &
                c_phi**2*pref_del2*del2phi_doublet(:,i+1)
        enddo

      ! end if lphi_doublet
      else
        ! dphi/dt = dphi
        df(l1:l2,m,n,iphi)=df(l1:l2,m,n,iphi)+p%dphi
        df(l1:l2,m,n,idphi)=df(l1:l2,m,n,idphi) - &
              pref_Hubble*Hscript*p%dphi-pref_Vprime*Vprime
        if (c_phi/=0 .and. .not. lphi_hom) then
          call del2(f, iphi, del2phi)
          df(l1:l2,m,n,idphi)=df(l1:l2,m,n,idphi) + c_phi**2*pref_del2*del2phi
        endif
        ! second scalar if lwaterfall
        if (lwaterfall) then
          df(l1:l2,m,n,ipsi)=df(l1:l2,m,n,ipsi)+p%dpsi
          df(l1:l2,m,n,idpsi)=df(l1:l2,m,n,idpsi) - &
                pref_Hubble*Hscript*p%dpsi-pref_Vprime*Vprimepsi
          if (c_psi/=0 .and. .not. lpsi_hom) then
            call del2(f, ipsi, del2psi)
            df(l1:l2,m,n,idpsi)=df(l1:l2,m,n,idpsi) + c_psi**2*pref_del2*del2psi
          endif
        endif
!
!  magnetic terms, add (alpf/a^2)*(E.B) to dphi'/dt equation
!  only if no lphi_doublet
!
        if (lmagnetic .and. lem_backreact) then
          if (lphi_hom .and. .not. lphi_linear_regime) then
            ! if (lconf_time) then
            !   df(l1:l2,m,n,idphi)=df(l1:l2,m,n,idphi)+alpf*edotbm_all*a21
            ! else
            !   df(l1:l2,m,n,idphi)=df(l1:l2,m,n,idphi)+alpf*edotbm_all*a21**2
            ! endif
            df(l1:l2,m,n,idphi)=df(l1:l2,m,n,idphi)+pref_alpf*alpf*edotbm_all
          endif
          if (.not. lphi_hom) then
            call dot_mn(p%el,p%bb,tmp)
            df(l1:l2,m,n,idphi)=df(l1:l2,m,n,idphi)+pref_alpf*alpf*tmp
          endif
          ! backreaction in second scalar psi if lwaterfall
          if (lwaterfall) then
            if (lpsi_hom) then
              df(l1:l2,m,n,idpsi)=df(l1:l2,m,n,idpsi)+pref_alpf*alpfpsi*edotbm_all
            else
              call dot_mn(p%el,p%bb,tmp)
              df(l1:l2,m,n,idpsi)=df(l1:l2,m,n,idpsi)+pref_alpf*alpfpsi*tmp
            endif
          endif
        endif
      endif
!
!  Total contribution to the timestep.
!  If Ndiv=0 is set, we compute instead an advective timestep based on the Alfven speed.
!  vA=B/sqrt(rho_chi), so dt=C_M*dx/vA. In practice, C_M (=cdt_rho_chi) can be 20.
!
      if (lfirst.and.ldt.and.ldt_klein_gordon) then
        if (Ndiv==0.) then
          if (lrho_chi) then
            advec2=advec2+(b2m_all/f_ode(iinfl_rho_chi))*dxyz_2/cdt_rho_chi**2
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
    subroutine read_sums_from_device

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

    endsubroutine read_sums_from_device
!***********************************************************************
    subroutine dspecial_dt_ode
!
      use SharedVariables, only: get_shared_variable
!
      if(lgpu) call read_sums_from_device
      call get_Hscript_and_a2(Hscript,a2rhom_all)
      if (lflrw) df_ode(ilna)=df_ode(ilna)+Hscript
!
!  Energy density of the charged particles.
!  This is currently only done for <sigE>*<E^2>, and not for <sigE*E^2>.
!
      if (lrho_chi) then
        if (lnoncollinear_EB .or. lnoncollinear_EB_aver .or. &
          lcollinear_EB .or. lcollinear_EB_aver) then
          df_ode(iinfl_rho_chi)=df_ode(iinfl_rho_chi)-4.*Hscript*f_ode(iinfl_rho_chi) &
            +(sigEm_all*e2m_all+sigBm_all*edotbm_all)/ascale**3
        else
          df_ode(iinfl_rho_chi)=df_ode(iinfl_rho_chi)-4.*Hscript*f_ode(iinfl_rho_chi)
        endif
      endif
!
!  Diagnostics
!
      if (.not. lmultithread) then
        sigEm_all_diagnos = sigEm_all
        sigBm_all_diagnos = sigBm_all
        call calc_ode_diagnostics_special(f_ode)
      endif

    endsubroutine dspecial_dt_ode
!***********************************************************************
    subroutine calc_ode_diagnostics_special(f_ode)
  
      use Diagnostics 
      
      real, dimension(max_n_odevars), intent(in) :: f_ode
      real :: rho_chi, lnascale
      real :: Hscript_diagnos

      if (lrho_chi) then
        rho_chi=f_ode(iinfl_rho_chi)
      else
        rho_chi=0.
      endif

      if (ldiagnos) then
        call get_Hscript_and_a2(Hscript_diagnos,a2rhom_all_diagnos)
        if(lflrw) lnascale=f_ode(ilna)
        call save_name(Hscript_diagnos,idiag_Hscriptm)
        call save_name(lnascale,idiag_lnam)
        call save_name(ddotam_all_diagnos,idiag_ddotam)
        call save_name(a2rhopm_all_diagnos,idiag_a2rhopm)
        call save_name(a2rhom_all_diagnos,idiag_a2rhom)
        call save_name(a2rhophim_all_diagnos,idiag_a2rhophim)
        call save_name(a2rhogphim_all_diagnos,idiag_a2rhogphim)
        call save_name(a2rhopsim_all_diagnos,idiag_a2rhopsim)
        call save_name(a2rhogpsim_all_diagnos,idiag_a2rhogpsim)
        call save_name(rho_chi,idiag_rho_chi)
        call save_name(sigEm_all_diagnos,idiag_sigEma)
        call save_name(sigBm_all_diagnos,idiag_sigBma)
        if (lnoncollinear_EB_aver .or. lcollinear_EB_aver) &
          call save_name(count_eb0_all,idiag_count_eb0a)
      endif

    endsubroutine calc_ode_diagnostics_special
!***********************************************************************
    subroutine calc_diagnostics_special(f,p)

      use Diagnostics
      real, dimension(mx,my,mz,mfarray) :: f
      type(pencil_case) :: p

      if (ldiagnos) then
        call sum_mn_name(p%phi,idiag_phim)
        if (idiag_phi2m/=0) call sum_mn_name(p%phi**2,idiag_phi2m)
        if (idiag_phirms/=0) call sum_mn_name(p%phi**2,idiag_phirms,lsqrt=.true.)
        call sum_mn_name(p%dphi,idiag_dphim)
        if (idiag_dphi2m/=0) call sum_mn_name(p%dphi**2,idiag_dphi2m)
        if (idiag_dphirms/=0) call sum_mn_name(p%dphi**2,idiag_dphirms,lsqrt=.true.)
        if (lwaterfall) then
          call sum_mn_name(p%psi,idiag_psim)
          if (idiag_psi2m/=0) call sum_mn_name(p%psi**2,idiag_psi2m)
          if (idiag_psirms/=0) call sum_mn_name(p%psi**2,idiag_psirms,lsqrt=.true.)
          call sum_mn_name(p%dpsi,idiag_dpsim)
          if (idiag_dpsi2m/=0) call sum_mn_name(p%dpsi**2,idiag_dpsi2m)
          if (idiag_dpsirms/=0) call sum_mn_name(p%dpsi**2,idiag_dpsirms,lsqrt=.true.)
        endif
      endif

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
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_phim=0; idiag_phi2m=0; idiag_phirms=0
        idiag_dphim=0; idiag_dphi2m=0; idiag_dphirms=0
        idiag_Hscriptm=0; idiag_lnam=0; idiag_ddotam=0
        idiag_a2rhopm=0; idiag_a2rhom=0; idiag_a2rhophim=0
        idiag_a2rhogphim=0; idiag_rho_chi=0; idiag_sigEma=0
        idiag_sigBma=0; idiag_count_eb0a=0
        idiag_psim=0; idiag_psi2m=0; idiag_psirms=0
        idiag_dpsim=0; idiag_dpsi2m=0; idiag_dpsirms=0
        idiag_a2rhogpsim=0; idiag_a2rhopsim=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'phim',idiag_phim)
        call parse_name(iname,cname(iname),cform(iname),'phi2m',idiag_phi2m)
        call parse_name(iname,cname(iname),cform(iname),'phirms',idiag_phirms)
        call parse_name(iname,cname(iname),cform(iname),'dphim',idiag_dphim)
        call parse_name(iname,cname(iname),cform(iname),'dphi2m',idiag_dphi2m)
        call parse_name(iname,cname(iname),cform(iname),'dphirms',idiag_dphirms)
        call parse_name(iname,cname(iname),cform(iname),'Hscriptm',idiag_Hscriptm)
        call parse_name(iname,cname(iname),cform(iname),'lnam',idiag_lnam)
        call parse_name(iname,cname(iname),cform(iname),'ddotam',idiag_ddotam)
        call parse_name(iname,cname(iname),cform(iname),'a2rhopm',idiag_a2rhopm)
        call parse_name(iname,cname(iname),cform(iname),'a2rhom',idiag_a2rhom)
        call parse_name(iname,cname(iname),cform(iname),'a2rhophim',idiag_a2rhophim)
        call parse_name(iname,cname(iname),cform(iname),'a2rhogphim',idiag_a2rhogphim)
        call parse_name(iname,cname(iname),cform(iname),'rho_chi',idiag_rho_chi)
        call parse_name(iname,cname(iname),cform(iname),'sigEma',idiag_sigEma)
        call parse_name(iname,cname(iname),cform(iname),'sigBma',idiag_sigBma)
        call parse_name(iname,cname(iname),cform(iname),'count_eb0a',idiag_count_eb0a)
        call parse_name(iname,cname(iname),cform(iname),'psim',idiag_psim)
        call parse_name(iname,cname(iname),cform(iname),'psi2m',idiag_psi2m)
        call parse_name(iname,cname(iname),cform(iname),'psirms',idiag_psirms)
        call parse_name(iname,cname(iname),cform(iname),'dpsim',idiag_dpsim)
        call parse_name(iname,cname(iname),cform(iname),'dpsi2m',idiag_dpsi2m)
        call parse_name(iname,cname(iname),cform(iname),'dpsirms',idiag_dpsirms)
        call parse_name(iname,cname(iname),cform(iname),'a2rhopsim',idiag_a2rhopsim)
        call parse_name(iname,cname(iname),cform(iname),'a2rhogpsim',idiag_a2rhogpsim)
      enddo
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_echarge

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
            echarge=1./sqrt(1./.35**2+41./(48.*pi**2)*alog(mass_zboson/energy_scale))
        endselect
      else
        echarge=echarge_const
      endif

    endsubroutine get_echarge
!***********************************************************************
    subroutine get_sigE_and_B

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

    endsubroutine get_sigE_and_B
!***********************************************************************
    subroutine prep_rhs_special
!
!  1-aug-25/TP: coded
!
      call get_Hscript_and_a2(Hscript,a2rhom_all)
      call get_echarge
      call get_sigE_and_B

    endsubroutine prep_rhs_special
!***********************************************************************
    subroutine get_a2
  
      real :: lnascale
      if(lflrw) then
        lnascale=f_ode(ilna)
        ascale=exp(lnascale)
      endif
      a2=ascale**2
      a21=1./a2

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
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real :: sigE1m,sigB1m
!
!  If requested, calculate here <dphi**2+gphi**2+(4./3.)*(E^2+B^2)/a^2>.
!  This needs to be done on all processors, because otherwise ascale
!  is not known on all processors.
!
      call get_a2
      call mpibcast_real(a2)
      call mpibcast_real(a21)
!
!  In the following loop, go through all penciles and add up results to get e2m, etc.
!
      ddotam=0.; a2rhopm=0.; a2rhom=0.; e2m=0; b2m=0; edotbm=0; a2rhophim=0.; a2rhogphim=0.
      sigE1m=0.; sigB1m=0.; a2rhogpsim=0.; a2rhopsim=0.
!
!  In the following, sum over all mn pencils.
!
      do n=n1,n2
      do m=m1,m2
        call prep_ode_right(f,sigE1m,sigB1m)
      enddo
      enddo
!
      a2rhopm=a2rhopm/nwgrid
      a2rhom=a2rhom/nwgrid
      a2rhophim=a2rhophim/nwgrid
      a2rhogphim=a2rhogphim/nwgrid
      a2rhopsim=a2rhopsim/nwgrid
      a2rhogpsim=a2rhogpsim/nwgrid
      ddotam=(four_pi_over_three/nwgrid)*ddotam
      if (lphi_hom .or. lrho_chi .or. lnoncollinear_EB .or. &
          lnoncollinear_EB_aver .or. lpsi_hom) then
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
      call mpireduce_sum(a2rhopm,a2rhopm_all)
      call mpiallreduce_sum(a2rhom,a2rhom_all)
      call mpireduce_sum(a2rhophim,a2rhophim_all)
      call mpireduce_sum(a2rhogphim,a2rhogphim_all)
      call mpireduce_sum(a2rhopsim,a2rhopsim_all)
      call mpireduce_sum(a2rhogpsim,a2rhogpsim_all)
      call mpiallreduce_sum(ddotam,ddotam_all)
      a2rhom_all_diagnos     = a2rhom_all
      a2rhopm_all_diagnos    = a2rhopm_all
      a2rhophim_all_diagnos  = a2rhophim_all
      a2rhopsim_all_diagnos  = a2rhopsim_all
      a2rhogphim_all_diagnos = a2rhogphim_all
      a2rhogpsim_all_diagnos = a2rhogpsim_all
      ddotam_all_diagnos     = ddotam_all

      if (lroot .and. lflrw) then
        call get_Hscript_and_a2(Hscript,a2rhom_all)
      endif
!
!  Broadcast to other processors, and each processor uses put_shared_variable
!  to get the values to other subroutines.
!
      call mpibcast_real(Hscript)
      call mpibcast_real(e2m_all)
      call mpibcast_real(b2m_all)

    endsubroutine special_after_boundary
!***********************************************************************
    subroutine prep_ode_right(f,sigE1m,sigB1m)
!
      use Sub, only: dot2_mn, grad, curl, dot_mn
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, intent(inout) :: sigE1m,sigB1m
      real, dimension (nx,3) :: el, bb, gphi, gpsi
      real, dimension (nx) :: e2, b2, gphi2, dphi, a2rhop, a2rho
      real, dimension (nx) :: ddota, phi, Vpotential, edotb, sigE1, sigB1
      real, dimension (nx) :: boost, gam_EB, eprime, bprime, jprime1
      real, dimension (nx) :: psi, dpsi, gpsi2, a2rhopsi_tmp
!
!  if requested, calculate here <dphi**2+gphi**2+(4./3.)*(E^2+B^2)/a^2>
!  rhop is purely an output quantity
!  It is called a2rhom because rhom=a2rhom/a^2.
!
      phi=f(l1:l2,m,n,iphi)
      dphi=f(l1:l2,m,n,idphi)
      a2rho=0.5*dphi**2
      a2rhop=dphi**2
      ! if (lphi_hom) then
      !   a2rhop=dphi**2
      !   !a2rho=0.5*dphi**2
      !   !a2rhophim=a2rhophim+sum(a2rho)
      ! else
      if (.not. lphi_hom) then
        call grad(f,iphi,gphi)    !MR: the ghost zones are not necessarily updated!!!
        ! alberto: this function is called from special_after_boundary so shouldn't have the ghost zones updated?
        call dot2_mn(gphi,gphi2)
        a2rhogphim=a2rhogphim+sum(0.5*gphi2)
        a2rhop=a2rhop+onethird*gphi2
        a2rho=a2rho+0.5*gphi2
        !a2rhophim=a2rhophim+sum(a2rho)
      endif
      a2rhophim=a2rhophim+sum(a2rho)

      if (lwaterfall) then
        psi=f(l1:l2,m,n,ipsi)
        dpsi=f(l1:l2,m,n,idpsi)
        a2rhopsi_tmp=0.5*dpsi**2
        a2rhop=a2rhop+dpsi**2

        if (.not. lpsi_hom) then
          call grad(f,ipsi,gpsi)
          call dot2_mn(gpsi,gpsi2)
          a2rhogpsim=a2rhogpsim+sum(0.5*gpsi2)
          a2rhop=a2rhop+onethird*gpsi2
          a2rhopsi_tmp=a2rhopsi_tmp+0.5*gpsi2
          !a2rhophim=a2rhophim+sum(a2rho)
        endif
        a2rhopsim=a2rhopsim+sum(a2rhopsi_tmp)
        a2rho=a2rho+a2rhopsi_tmp
      endif
!
!  Note the .5*fourthird factor in front of (e2+b2)*a21, but that is
!  just for rhop, which is output quantity.
!
      if (iex/=0 .and. lem_backreact) then
        el=f(l1:l2,m,n,iex:iez)
        call curl(f,iaa,bb)          !MR: the ghost zones are not necessarily updated!!!
        call dot2_mn(bb,b2)
        call dot2_mn(el,e2)
        a2rhop=a2rhop+(.5*fourthird)*(e2+b2)*a21
        if (.not. lphi_linear_regime) a2rho=a2rho+.5*(e2+b2)*a21
        if (lrho_chi) then
          a2rho=a2rho+scale_rho_chi_Heqn*a2*f_ode(iinfl_rho_chi)
        endif
      endif
!
      a2rhopm=a2rhopm+sum(a2rhop)
!
!  Choice of different potentials
!
      select case (Vprime_choice)
        case ('quadratic')  ; Vpotential=.5*phimass2*phi**2
        case ('quartic')    ; Vpotential=phimass2*phi+(lambda_phi/6.)*phi**3  !(to be corrected)
        case ('cos-profile'); Vpotential=phimass2*lambda_phi*sin(lambda_phi*phi)  !(to be corrected)
        case ('waterfall')
          Vpotential=0.5*phimass2*phi**2 + .25*lambda_psi*psi**4
          if (lambda_psi /= 0.) then
            Vpotential=Vpotential+.25*psimass2**2/lambda_psi
          endif
          Vpotential=Vpotential-0.5*psimass2*psi**2+.5*coupl_phipsi**2*phi**2*psi**2
        case default
          call fatal_error("special_after_boundary: No such Vprime_choice: ",trim(Vprime_choice))
      endselect
!
!  compute ddotam = a"/a (needed for GW module)
!
      ddota=-dphi**2+4.*a2*Vpotential
      !if (lphi_hom) then
      !  ddota=-dphi**2+4.*a2*Vpotential
      !else
      if (.not. lphi_hom) ddota=ddota-gphi2
      !endif
      ! contribution from second scalar field
      if (lwaterfall) then
        ddota=ddota-dpsi**2
        if (.not. lpsi_hom) ddota=ddota-gpsi2
      endif
      ddotam=ddotam+sum(ddota)
      a2rho=a2rho+a2*Vpotential
      a2rhom=a2rhom+sum(a2rho)
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
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr
    use General , only: string_to_enum

    integer, parameter :: n_pars=100
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call string_to_enum(enum_hscript_choice,hscript_choice)
    call string_to_enum(enum_vprime_choice,vprime_choice)
    call copy_addr(enum_hscript_choice,p_par(1)) ! int
    call copy_addr(enum_vprime_choice,p_par(2)) ! int
    call copy_addr(iphi,p_par(3)) ! int
    call copy_addr(idphi,p_par(4)) ! int
    call copy_addr(lnoncollinear_eb,p_par(5)) ! bool
    call copy_addr(lnoncollinear_eb_aver,p_par(6)) ! bool
    call copy_addr(lcollinear_eb,p_par(7)) ! bool
    call copy_addr(lcollinear_eb_aver,p_par(8)) ! bool
    call string_to_enum(enum_echarge_type,echarge_type)
    call copy_addr(enum_echarge_type,p_par(9)) ! int
    call copy_addr(echarge_const,p_par(10))
    call copy_addr(hscript0,p_par(11))
    call copy_addr(phimass2,p_par(12))
    call copy_addr(lambda_phi,p_par(13))
    call copy_addr(lconf_time,p_par(14)) ! bool
    call copy_addr(c_phi,p_par(15)) 
    call copy_addr(lem_backreact,p_par(16))
    call copy_addr(ldt_klein_gordon,p_par(17)) ! bool
    call copy_addr(ndiv,p_par(18)) ! int
    call copy_addr(lrho_chi,p_par(19)) ! bool
    call copy_addr(iinfl_rho_chi,p_par(20)) ! int
    call copy_addr(cdt_rho_chi,p_par(21))
    call copy_addr(lflrw,p_par(22)) ! bool
    call copy_addr(ilna,p_par(23)) ! int
    call copy_addr(scale_rho_chi_heqn,p_par(24))
    call copy_addr(ipsi,p_par(25)) ! int
    call copy_addr(idpsi,p_par(26)) ! int
    call copy_addr(iphi_up_re,p_par(27)) ! int
    call copy_addr(iphi_up_im,p_par(28)) ! int
    call copy_addr(iphi_down_re,p_par(29)) ! int
    call copy_addr(iphi_down_im,p_par(30)) ! int
    call copy_addr(idphi_up_re,p_par(31)) ! int
    call copy_addr(idphi_up_im,p_par(32)) ! int
    call copy_addr(idphi_down_re,p_par(33)) ! int
    call copy_addr(idphi_down_im,p_par(34)) ! int
    call copy_addr(psimass2,p_par(35))
    call copy_addr(lambda_psi,p_par(36))
    call copy_addr(coupl_phipsi,p_par(37))
    call copy_addr(c_psi,p_par(38))
    call copy_addr(eta_phi,p_par(45))
    call copy_addr(lhiggs_friction,p_par(46)) ! bool
    call copy_addr(lwaterfall,p_par(47)) ! bool
    call copy_addr(higgs_friction,p_par(48))
    call copy_addr(lphi_doublet,p_par(49)) ! bool
    call copy_addr(lphi_weakcharge,p_par(50)) ! bool
    call copy_addr(lphi_hypercharge,p_par(51)) ! bool
    call copy_addr(ia0,p_par(52)) ! int
    call copy_addr(iw0,p_par(53)) ! int


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
