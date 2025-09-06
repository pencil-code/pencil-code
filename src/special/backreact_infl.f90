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
  integer :: iinfl_rho_chi=0
  real :: ncutoff_phi=1., infl_v=.1
  real :: axionmass=1.06e-6, axionmass2, ascale_ini=1.
  real :: phi0=.44, dphi0=-1.69e-7, c_light_axion=1., lambda_axion=0., eps=.01
  real :: amplphi=.1, ampldphi=.0, kx_phi=1., ky_phi=0., kz_phi=0., phase_phi=0., width=.1, offset=0.
  real :: initpower_phi=0.,  cutoff_phi=0.,  initpower2_phi=0.
  real :: initpower_dphi=0., cutoff_dphi=0., initpower2_dphi=0.
  real :: kgaussian_phi=0.,kpeak_phi=0., kgaussian_dphi=0., kpeak_dphi=0.
  real :: relhel_phi=0.
  real :: ddotam, a2rhopm, a2rhopm_all, a2rhom, a2rhom_all
  real :: edotbm, edotbm_all, e2m, e2m_all, b2m, b2m_all, a2rhophim, a2rhophim_all
  real :: sigE1m_all_nonaver, sigB1m_all_nonaver,sigEm_all,sigBm_all,sigEm_all_diagnos,sigBm_all_diagnos
  real :: a2rhogphim, a2rhogphim_all
  real :: a2, a21, Hscript
  real :: Hscript0=0., scale_rho_chi_Heqn=1., rho_chi_init=0., cdt_rho_chi=1.
  real :: amplee_BD_prefactor=0., deriv_prefactor_ee=-1.
  real :: echarge=.0, echarge_const=.303
  real :: count_eb0_all=0.
  real, target :: ddotam_all
  real, pointer :: alpf
  real, pointer :: sigE_prefactor, sigB_prefactor, mass_chi
  real, dimension (nx) :: dt1_special
  logical :: lcompute_dphi0=.true.
  logical :: lbackreact_infl=.true., lem_backreact=.true., lzeroHubble=.false.
  logical :: lscale_tobox=.true.,ldt_backreact_infl=.true., lconf_time=.true.
  logical :: lskip_projection_phi=.false., lvectorpotential=.false., lflrw=.false.
  logical :: lrho_chi=.false., lno_noise_phi=.false., lno_noise_dphi=.false.
  logical, pointer :: lphi_hom, lphi_linear_regime, lnoncollinear_EB, lnoncollinear_EB_aver
  logical, pointer :: lcollinear_EB, lcollinear_EB_aver, lmass_suppression
  logical, pointer :: lallow_bprime_zero
  character (len=labellen) :: Vprime_choice='quadratic', Hscript_choice='default'
  character (len=labellen), dimension(ninit) :: initspecial='nothing'
  character (len=50) :: echarge_type='const', init_rho_chi='zero'
!
  namelist /special_init_pars/ &
      initspecial, phi0, dphi0, axionmass, eps, ascale_ini, &
      lcompute_dphi0, lem_backreact, &
      c_light_axion, lambda_axion, amplphi, ampldphi, lno_noise_phi, lno_noise_dphi, &
      kx_phi, ky_phi, kz_phi, phase_phi, width, offset, &
      initpower_phi, initpower2_phi, cutoff_phi, kgaussian_phi, kpeak_phi, &
      initpower_dphi, initpower2_dphi, cutoff_dphi, kpeak_dphi, &
      ncutoff_phi, lscale_tobox, Hscript0, Hscript_choice, infl_v, lflrw, &
      lrho_chi, scale_rho_chi_Heqn, amplee_BD_prefactor, deriv_prefactor_ee, &
      echarge_type, init_rho_chi, rho_chi_init
!
  namelist /special_run_pars/ &
      initspecial, phi0, dphi0, axionmass, eps, ascale_ini, &
      lbackreact_infl, lem_backreact, c_light_axion, lambda_axion, Vprime_choice, &
      !lem_backreact, c_light_axion, lambda_axion, Vprime_choice, &
      lzeroHubble, ldt_backreact_infl, Ndiv, Hscript0, Hscript_choice, infl_v, &
      lflrw, lrho_chi, scale_rho_chi_Heqn, echarge_type, cdt_rho_chi
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_phim=0      ! DIAG_DOC: $\left<\phi\right>$
  integer :: idiag_phi2m=0     ! DIAG_DOC: $\left<\phi^2\right>$
  integer :: idiag_phirms=0    ! DIAG_DOC: $\left<\phi^2\right>^{1/2}$
  integer :: idiag_dphim=0     ! DIAG_DOC: $\left<\phi'\right>$
  integer :: idiag_dphi2m=0    ! DIAG_DOC: $\left<(\phi')^2\right>$
  integer :: idiag_dphirms=0   ! DIAG_DOC: $\left<(\phi')^2\right>^{1/2}$
  integer :: idiag_Hscriptm=0   ! DIAG_DOC: $\left<{\cal a*H}\right>$
  integer :: idiag_lnam=0      ! DIAG_DOC: $\left<\ln a\right>$
  integer :: idiag_ddotam=0    ! DIAG_DOC: $a''/a$
  integer :: idiag_a2rhopm=0   ! DIAG_DOC: $a^2 (rho+p)$
  integer :: idiag_a2rhom=0    ! DIAG_DOC: $a^2 rho$
  integer :: idiag_a2rhophim=0  ! DIAG_DOC: $a^2 rho$
  integer :: idiag_a2rhogphim=0 ! DIAG_DOC: $0.5 <grad phi^2>$
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
     if (lflrw) then
!     call farray_register_ode('infl_hubble',iinfl_hubble)
       call farray_register_ode('infl_lna',iinfl_lna)
     endif
!
     if (lrho_chi) then
       call farray_register_ode('infl_rho_chi',iinfl_rho_chi)
     endif
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
!
      if (lmagnetic .and. lem_backreact) then
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
        call get_shared_variable('mass_chi',mass_chi)
      else
        if (.not.associated(alpf)) allocate(alpf,lphi_hom,lphi_linear_regime, &
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
!  SAMPLE IMPLEMENTATION
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
            if (lroot .and. lflrw) then
              f_ode(iinfl_lna) =lnascale
!              f(iinfl_hubble) =Hubble_ini
            endif
!
          case ('default')
            Vpotential=.5*axionmass2*phi0**2
            Hubble_ini=sqrt(8.*pi/3.*(.5*axionmass2*phi0**2*ascale_ini**2))
!            dphi0=-ascale_ini*sqrt(2*eps/3.*Vpotential)
            if (lcompute_dphi0) dphi0=-sqrt(1/(12.*pi))*axionmass*ascale_ini
           ! dphi0=-sqrt(1/(12.*pi))*axionmass*ascale_ini
           ! dphi0=-sqrt(16*pi/3)*axionmass*ascale_ini
            tstart=-1/(ascale_ini*Hubble_ini)
            t=tstart
            lnascale=log(ascale_ini)
            f(:,:,:,iinfl_phi)   =f(:,:,:,iinfl_phi)   +phi0
            f(:,:,:,iinfl_dphi)  =f(:,:,:,iinfl_dphi)  +dphi0
            if (lroot .and. lflrw) then
              f_ode(iinfl_lna)   =lnascale
              a2                 =exp(f_ode(iinfl_lna))**2
              Hscript            =Hubble_ini/exp(lnascale)
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
!
          case ('Bunch-Davies')
            if (lroot) print*,'Hubble_ini=',Hubble_ini
            amplphi_BD=amplphi*Hubble_ini
            deriv_prefactor=1.
            call bunch_davies(f,iinfl_phi,iinfl_phi,iinfl_dphi,iinfl_dphi,amplphi_BD,kpeak_phi,deriv_prefactor)
            if (amplee_BD_prefactor/=0.) then
              deriv_prefactor=deriv_prefactor_ee
              amplee_BD=amplee_BD_prefactor*Hubble_ini
              call bunch_davies(f,iax,iaz,iex,iez,amplee_BD,kpeak_phi,deriv_prefactor)
            endif
          case default
            call fatal_error("init_special: No such initspecial: ", trim(initspecial(j)))
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
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
!      if (lmagnetic .and. lbackreact_infl) lpenc_requested(i_infl_a21)=.true.
!
!  pencil for gradient of phi
!
      lpenc_requested(i_gphi)=.true.
      if (lmagnetic .and. lem_backreact) lpenc_requested(i_infl_dphi)=.true.
!
!  Magnetic field needed for Maxwell stress
!
      if (lmagnetic) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_el)=.true.
        if (lrho_chi .or. lnoncollinear_EB .or. lnoncollinear_EB_aver .or. &
          lcollinear_EB .or. lcollinear_EB_aver) lpenc_requested(i_e2)=.true.
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
! infl_gphi
      if (lpencil(i_gphi)) call grad(f,iinfl_phi,p%gphi)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine get_Hscript_and_a2(Hscript,a2rhom_all)
      real, intent(OUT) :: Hscript
      real, intent(IN)  :: a2rhom_all
!
!  Choice of prescription for Hscript
!
      select case (Hscript_choice)
        case ('default')
          Hscript=sqrt((8.*pi/3.)*a2rhom_all)
          if (lgpu) call get_a2
        case ('set')
          Hscript=Hscript0
          a2=1.
          a21=1./a2
        case default
          call fatal_error("dspecial_dt: No such Hscript_choice: ", trim(Hscript_choice))
      endselect

!  Possibility of turning off evolution of scale factor and Hubble parameter
!  By default, lzeroHubble=F, so we use the calculation from above.
!
      if (lzeroHubble) then
        a2=1.
        a21=1./a2
        Hscript=0.
      endif
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
      use Sub, only: dot_mn, del2
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: phi, dphi, Vprime
      real, dimension (nx) :: tmp, del2phi
!      real :: tmp2
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!
      phi=f(l1:l2,m,n,iinfl_phi)
      dphi=f(l1:l2,m,n,iinfl_dphi)
!
!  Choice of different potentials.
!  For the 1-cos profile, -Vprime (on the rhs) enters with -sin().
!
      select case (Vprime_choice)
        case ('quadratic'); Vprime=axionmass2*phi
        case ('quartic'); Vprime=axionmass2*phi+(lambda_axion/6.)*phi**3
        case ('cos-profile'); Vprime=axionmass2*lambda_axion*sin(lambda_axion*phi)
        case default
          call fatal_error("dspecial_dt: No such Vprime_choice: ", trim(Vprime_choice))
      endselect
!
!  Update df.
!  dphi/dt = psi
!  dpsi/dt = - ...
!
        df(l1:l2,m,n,iinfl_phi)=df(l1:l2,m,n,iinfl_phi)+f(l1:l2,m,n,iinfl_dphi)
        if (lconf_time) then
          df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)-2.*Hscript*dphi-a2*Vprime
        else
          df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)-2.*Hscript*dphi-Vprime
        endif
!
!  speed of light term
!
        if (c_light_axion/=0. .and. .not. lphi_hom) then
          call del2(f,iinfl_phi,del2phi)
          if (lconf_time) then
            df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)+c_light_axion**2*del2phi
          else
            df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)+c_light_axion**2*a21*del2phi
          endif
        endif
!
!  magnetic terms, add (alpf/a^2)*(E.B) to dphi'/dt equation
!
      if (lmagnetic .and. lem_backreact) then
        if (lconf_time) then
          if (lphi_hom) then
            if (.not. lphi_linear_regime) &
              df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)+alpf*edotbm_all*a21
          else
            call dot_mn(p%el,p%bb,tmp)
            df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)+alpf*tmp*a21
          endif
        else
          if (lphi_hom) then
            if (.not. lphi_linear_regime) &
              df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)+alpf*edotbm_all*a21**2
          else
            call dot_mn(p%el,p%bb,tmp)
            df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)+alpf*tmp*a21**2
          endif
        endif
      endif
!
!  Total contribution to the timestep.
!  If Ndiv=0 is set, we compute instead an advective timestep based on the Alfven speed.
!  vA=B/sqrt(rho_chi), so dt=C_M*dx/vA. In practice, C_M (=cdt_rho_chi) can be 20.
!
      if (lfirst.and.ldt.and.ldt_backreact_infl) then
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
!      if (lfirst.and.ldt.and.ldt_backreact_infl) then
!        tmp2 = axionmass*sqrt(a2)
!        if (tmp2 > Hscript) then
!          dt1_special = Ndiv*abs(tmp2)
!        else
!          dt1_special = Ndiv*abs(Hscript)
!        endif
!        dt1_max=max(dt1_max,dt1_special)
!      endif
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
!     use Magnetic, only: eta_xtdep
!
!
      if(lgpu) call read_sums_from_device
      call get_Hscript_and_a2(Hscript,a2rhom_all)
      if (lflrw) then
        df_ode(iinfl_lna)=df_ode(iinfl_lna)+Hscript
      endif
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
      
      real, dimension(max_n_odevars), intent(IN) :: f_ode
      real :: rho_chi, lnascale
      real :: Hscript_diagnos

      if (lrho_chi) then
        rho_chi=f_ode(iinfl_rho_chi)
      else
        rho_chi=0.
      endif

      if (ldiagnos) then
        call get_Hscript_and_a2(Hscript_diagnos,a2rhom_all_diagnos)
        if(lflrw) lnascale=f_ode(iinfl_lna)
        call save_name(Hscript_diagnos,idiag_Hscriptm)
        call save_name(lnascale,idiag_lnam)
        call save_name(ddotam_all_diagnos,idiag_ddotam)
        call save_name(a2rhopm_all_diagnos,idiag_a2rhopm)
        call save_name(a2rhom_all_diagnos,idiag_a2rhom)
        call save_name(a2rhophim_all_diagnos,idiag_a2rhophim)
        call save_name(a2rhogphim_all_diagnos,idiag_a2rhogphim)
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
      real, dimension(nx) :: dphi,phi

      if (ldiagnos) then
        dphi=f(l1:l2,m,n,iinfl_dphi)
        phi=f(l1:l2,m,n,iinfl_phi)
        call sum_mn_name(phi,idiag_phim)
        if (idiag_phi2m/=0) call sum_mn_name(phi**2,idiag_phi2m)
        if (idiag_phirms/=0) call sum_mn_name(phi**2,idiag_phirms,lsqrt=.true.)
        call sum_mn_name(dphi,idiag_dphim)
        if (idiag_dphi2m/=0) call sum_mn_name(dphi**2,idiag_dphi2m)
        if (idiag_dphirms/=0) call sum_mn_name(dphi**2,idiag_dphirms,lsqrt=.true.)
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
      enddo
!!
!!!  write column where which magnetic variable is stored
!!      if (lwr) then
!!        call farray_index_append('idiag_SPECIAL_DIAGNOSTIC',idiag_SPECIAL_DIAGNOSTIC)
!!      endif
!!
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
        lnascale=f_ode(iinfl_lna)
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
      sigE1m=0.; sigB1m=0.
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
      call mpireduce_sum(a2rhopm,a2rhopm_all)
      call mpiallreduce_sum(a2rhom,a2rhom_all)
      call mpireduce_sum(a2rhophim,a2rhophim_all)
      call mpireduce_sum(a2rhogphim,a2rhogphim_all)
      call mpiallreduce_sum(ddotam,ddotam_all)
      a2rhom_all_diagnos     = a2rhom_all
      a2rhopm_all_diagnos    = a2rhopm_all
      a2rhophim_all_diagnos  = a2rhophim_all
      a2rhogphim_all_diagnos = a2rhogphim_all
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
      real, intent(INOUT) :: sigE1m,sigB1m
      real, dimension (nx,3) :: el, bb, gphi
      real, dimension (nx) :: e2, b2, gphi2, dphi, a2rhop, a2rho
      real, dimension (nx) :: ddota, phi, Vpotential, edotb, sigE1, sigB1
      real, dimension (nx) :: boost, gam_EB, eprime, bprime, jprime1
!
!  if requested, calculate here <dphi**2+gphi**2+(4./3.)*(E^2+B^2)/a^2>
!  rhop is purely an output quantity
!  It is called a2rhom because rhom=a2rhom/a^2.
!
      phi=f(l1:l2,m,n,iinfl_phi)
      dphi=f(l1:l2,m,n,iinfl_dphi)
      if (lphi_hom) then
        a2rhop=dphi**2
        a2rho=0.5*dphi**2
        a2rhophim=a2rhophim+sum(a2rho)
      else
        call grad(f,iinfl_phi,gphi)    !MR: the ghost zones are not necessarily updated!!!
        call dot2_mn(gphi,gphi2)
        a2rhogphim=a2rhogphim+sum(0.5*gphi2)
        a2rhop=dphi**2+onethird*gphi2
        a2rho=0.5*(dphi**2+gphi2)
        a2rhophim=a2rhophim+sum(a2rho)
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
        case ('quadratic')  ; Vpotential=.5*axionmass2*phi**2
        case ('quartic')    ; Vpotential=axionmass2*phi+(lambda_axion/6.)*phi**3  !(to be corrected)
        case ('cos-profile'); Vpotential=axionmass2*lambda_axion*sin(lambda_axion*phi)  !(to be corrected)
        case default
          call fatal_error("special_after_boundary: No such Vprime_choice: ",trim(Vprime_choice))
      endselect
!
!  compute ddotam = a"/a (needed for GW module)
!
      if (lphi_hom) then
        ddota=-dphi**2+4.*a2*Vpotential
      else
        ddota=-dphi**2-gphi2+4.*a2*Vpotential
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
    call copy_addr(lem_backreact,p_par(17))
    call copy_addr(ldt_backreact_infl,p_par(18)) ! bool
    call copy_addr(ndiv,p_par(19)) ! int
    call copy_addr(lrho_chi,p_par(20)) ! bool
    call copy_addr(iinfl_rho_chi,p_par(21)) ! int
    call copy_addr(cdt_rho_chi,p_par(22))
    call copy_addr(lflrw,p_par(23)) ! bool
    call copy_addr(iinfl_lna,p_par(24)) ! int
    call copy_addr(scale_rho_chi_heqn,p_par(25))

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
