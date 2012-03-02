! $Id$
!
!  This module solves the Poisson equation for pressure for 
!  the linearized density
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ldensity = .false.
! CPARAM logical, parameter :: lanelastic = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 6
! COMMUNICATED AUXILIARIES 5
!
! PENCILS PROVIDED glnrho(3); grho(3); gpp(3); 
! PENCILS PROVIDED uglnrho; ugrho
! PENCILS PROVIDED glnrho2; del2lnrho; del2rho; del6lnrho; del6rho
! PENCILS PROVIDED hlnrho(3,3); sglnrho(3); uij5glnrho(3),transprho
! PENCILS PROVIDED ekin
! PENCILS PROVIDED rho; rho1; lnrho
!
!***************************************************************
module Density
!
  use Cparam
  use Cdata
  use Messages
  use EquationOfState
  use Sub
  use Diagnostics
!
  use Special
!
  implicit none
!
  include '../density.h'
!
  real, dimension (ninit) :: ampllnrho=0.0, widthlnrho=0.1
  real, dimension (ninit) :: rho_left=1.0, rho_right=1.0
  real, dimension (ninit) :: amplrho=0.0, phase_lnrho=0.0, radius_lnrho=0.5
  real, dimension (ninit) :: kx_lnrho=1.0, ky_lnrho=1.0, kz_lnrho=1.0
  real, dimension (ninit) :: kxx_lnrho=0.0, kyy_lnrho=0.0, kzz_lnrho=0.0
  real :: lnrho_const=0.0, rho_const=1.0
  real :: cdiffrho=0.0, diffrho=0.0, diffrho_hyper3=0.0, diffrho_shock=0.0
  real :: eps_planet=0.5, q_ell=5.0, hh0=0.0
  real :: xblob=0., yblob=0., zblob=0.
  real :: co1_ss=0.,co2_ss=0.,Sigma1=150.
  real :: lnrho_int=0.,lnrho_ext=0.,damplnrho_int=0.,damplnrho_ext=0.
  real :: wdamp=0.,density_floor=-1.0
  real :: radial_percent_smooth=10.,rshift=0.0
  real, dimension(3) :: diffrho_hyper3_aniso=0.
  real, dimension(mz) :: lnrho_init_z=0.0,del2lnrho_init_z=0.0
  real, dimension(mz) :: dlnrhodz_init_z=0.0, glnrho2_init_z=0.0
  real, target :: plaw=0.0
  real :: lnrho_z_shift=0.0
  real, dimension (nz,3) :: glnrhomz
  real :: powerlr=3.0, zoverh=1.5, hoverr=0.05
  real :: init_average_density
  integer, parameter :: ndiff_max=4
  logical :: lcontinuity_gas=.true.
  logical :: lupw_lnrho=.false.,lupw_rho=.false.
  logical :: ldiff_normal=.false.,ldiff_hyper3=.false.,ldiff_shock=.false.
  logical :: ldiff_hyper3lnrho=.false.,ldiff_hyper3_aniso=.false.
  logical :: ldiff_hyper3_polar=.false.,lanti_shockdiffusion=.false.
  logical :: lfreeze_lnrhoint=.false.,lfreeze_lnrhoext=.false.
  logical :: lfreeze_lnrhosqu=.false.,lexponential_smooth=.false.
  logical :: lrho_as_aux=.false., ldiffusion_nolog=.false.
  logical :: lshare_plaw=.false.,lmassdiff_fix=.false.
  logical :: lcheck_negative_density=.false.
  logical :: lcalc_glnrhomean=.false.
  logical, pointer :: lanelastic_lin

!
  character (len=labellen), dimension(ninit) :: initlnrho='nothing'
  character (len=labellen) :: strati_type='lnrho_ss'
  character (len=labellen), dimension(ndiff_max) :: idiff=''
  character (len=intlen) :: iinit_str
  complex :: coeflnrho=0.
!
  integer :: iglobal_gg=0
  integer :: niter=1
!
  namelist /density_init_pars/ &
      ampllnrho,initlnrho,widthlnrho,                    &
      rho_left,rho_right,lnrho_const,rho_const,cs2bot,cs2top,       &
      radius_lnrho,eps_planet,xblob,yblob,zblob,                    &
      b_ell,q_ell,hh0,rbound,lwrite_stratification,                 &
      mpoly,strati_type,beta_glnrho_global,radial_percent_smooth,   &
      kx_lnrho,ky_lnrho,kz_lnrho,amplrho,phase_lnrho,coeflnrho,     &
      kxx_lnrho, kyy_lnrho, kzz_lnrho,                              &
      co1_ss,co2_ss,Sigma1,idiff,ldensity_nolog,lexponential_smooth,&
      wdamp,plaw,lcontinuity_gas,density_floor,lanti_shockdiffusion,&
      rshift,lrho_as_aux,ldiffusion_nolog,lnrho_z_shift,            &
      lshare_plaw, powerlr, zoverh, hoverr
!
  namelist /density_run_pars/ &
      cdiffrho,diffrho,diffrho_hyper3,diffrho_shock,                &
      cs2bot,cs2top,lupw_lnrho,lupw_rho,idiff,                      &
      lnrho_int,lnrho_ext,damplnrho_int,damplnrho_ext,              &
      wdamp,lfreeze_lnrhoint,lfreeze_lnrhoext,                      &
      lnrho_const,plaw,lcontinuity_gas,                             &
      diffrho_hyper3_aniso,lfreeze_lnrhosqu,density_floor,          &
      lanti_shockdiffusion,lrho_as_aux,ldiffusion_nolog,            &
      lcheck_negative_density,lmassdiff_fix,niter
! diagnostic variables (need to be consistent with reset list below)
  integer :: idiag_rhom=0       ! DIAG_DOC: $\left<\varrho\right>$
                                ! DIAG_DOC:   \quad(mean density)
  integer :: idiag_rho2m=0      ! DIAG_DOC:
  integer :: idiag_lnrho2m=0    ! DIAG_DOC:
  integer :: idiag_drho2m=0     ! DIAG_DOC:
  integer :: idiag_drhom=0      ! DIAG_DOC:
  integer :: idiag_rhomin=0     ! DIAG_DOC:
  integer :: idiag_rhomax=0     ! DIAG_DOC:
  integer :: idiag_ugrhom=0     ! DIAG_DOC: $\left<\uv\cdot\nabla\varrho\right>$
  integer :: idiag_uglnrhom=0   ! DIAG_DOC:
  integer :: idiag_lnrhomphi=0  ! PHIAVG_DOC: $\left<\ln\varrho\right>_\varphi$
  integer :: idiag_rhomphi=0    ! PHIAVG_DOC: $\left<\varrho\right>_\varphi$
  integer :: idiag_dtd=0        ! DIAG_DOC:
  integer :: idiag_rhomz=0      ! DIAG_DOC:
  integer :: idiag_rhomy=0      ! DIAG_DOC:
  integer :: idiag_rhomx=0      ! DIAG_DOC:
  integer :: idiag_rhomxy=0     ! DIAG_DOC:
  integer :: idiag_rhomxz=0     ! DIAG_DOC:
  integer :: idiag_rhomr=0      ! DIAG_DOC:
  integer :: idiag_totmass=0    ! DIAG_DOC:
  integer :: idiag_mass=0       ! DIAG_DOC: $\int\varrho\,dV$
  integer :: idiag_divrhoum=0   ! DIAG_DOC: $\left<\nabla\cdot(\varrho\uv)\right>$
  integer :: idiag_divrhourms=0 ! DIAG_DOC: $\left|\nabla\cdot(\varrho\uv)\right|_{\rm rms}$
  integer :: idiag_divrhoumax=0 ! DIAG_DOC: $\left|\nabla\cdot(\varrho\uv)\right|_{\rm max}$
!
  contains
!***********************************************************************
    subroutine register_density()
!
!  Initialise variables which should know that we solve the
!  compressible hydro equations: ilnrho; increase nvar accordingly.
!
!   4-jun-02/axel: adapted from hydro
!
      use FArrayManager
!
      call farray_register_auxiliary('pp',ipp,communicated=.true.)
      call farray_register_auxiliary('rhs',irhs,vector=3,communicated=.true.)
      call farray_register_auxiliary('divu',idivu,communicated=.false.)
!
!  Identify version number (generated automatically by CVS).
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_density
!***********************************************************************
    subroutine initialize_density(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  For compatibility with other applications, we keep the possibility
!  of giving diffrho units of dxmin*cs0, but cs0 is not well defined general
!
!  24-nov-02/tony: coded
!  31-aug-03/axel: normally, diffrho should be given in absolute units
!
      use Deriv, only: der_pencil,der2_pencil
      use FArrayManager, only: farray_register_auxiliary, farray_register_global
      use Gravity, only: lnumerical_equilibrium
      use Mpicomm, only: stop_it
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
      integer :: i,ierr
      logical :: lnothing
!
      call get_shared_variable('lanelastic_lin',lanelastic_lin,ierr)
      if (ierr/=0) call stop_it("lanelastic_lin: "//&
        "there was a problem when sharing lanelastic_lin")
!
      if (lanelastic_lin) then
        call farray_register_auxiliary('rho_b',irho_b,communicated=.true.)
      else
        call farray_register_auxiliary('rho',irho,communicated=.true.)
      endif
!
!  initialize cs2cool to cs20
!  (currently disabled, because it causes problems with mdwarf auto-test)
!     cs2cool=cs20
!
      if (diffrho==0.) then
!
!  Made to work by adding diffrho + cdiffrho to the rprint reset list.
!
        diffrho=cdiffrho*dxmin*cs0
      endif
!
!  Turn off continuity equation. 
!
      lcontinuity_gas=.false.
      print*, 'initialize_density: density_anelastic, turned off continuity equation'
!
!  Initialize mass diffusion
!
      ldiff_normal=.false.
      ldiff_shock=.false.
      ldiff_hyper3=.false.
      ldiff_hyper3lnrho=.false.
      ldiff_hyper3_aniso=.false.
      ldiff_hyper3_polar=.false.
!
      lnothing=.false.
!
      do i=1,ndiff_max
        select case (idiff(i))
        case ('','none')
          if (lroot .and. (.not. lnothing)) print*,'no mass diffusion'
        case default
          write(unit=errormsg,fmt=*) 'initialize_density: ', &
              'You cannot have mass diffusion in anelastic approximation', trim(idiff(i))
          call fatal_error('initialize_density',errormsg)
        endselect
        lnothing=.true.
      enddo
!
! Tell the equation of state that we're here and what f variable we use
! DM+PC
! For anelastic case use pressure
        call select_eos_variable('pp',ipp)
!
! Do not allow inconsistency between rho0 (from eos) and rho_const 
! or lnrho0 and lnrho_const. 
!
      if (rho0/=rho_const) then
        if (lroot) then 
          print*,"WARNING!"
          print*,"inconsistency between the density constants from eos  "
          print*,"(rho0 or lnrho0) and the ones from the density module "
          print*,"(rho_const or lnrho_const). It may damage your        "
          print*,"simulation if you are using them in different places. "
          call warning("initialize_density","")
        endif
      endif
!
      if (lnumerical_equilibrium) then
         if (lroot) print*,'initializing global gravity in density'
         call farray_register_global('gg',iglobal_gg,vector=3)
      endif
!
    endsubroutine initialize_density
!***********************************************************************
    subroutine init_lnrho(f)
!
!  Initialise logarithmic or non-logarithmic density.
!
!   7-nov-01/wolf: coded
!  28-jun-02/axel: added isothermal
!  15-oct-03/dave: added spherical shell (kws)
!
      use General, only: itoa,complex_phase
      use Gravity, only: zref,z1,z2,gravz,nu_epicycle,potential, &
                         lnumerical_equilibrium
      use Initcond
      use IO
      use Mpicomm
!      use Selfgravity, only: rhs_poisson_const
      use InitialCondition, only: initial_condition_lnrho
      use SharedVariables, only: get_shared_variable
      use Poisson, only: inverse_laplacian
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: pot,prof
      real, dimension (nx,ny,nz) :: psi
      real, dimension (ninit) :: lnrho_left,lnrho_right
      real :: lnrhoint,cs2int,pot0
      real :: pot_ext,lnrho_ext,cs2_ext,tmp1,k_j2
      real :: zbot,ztop,haut
      real, dimension(1) :: mass_per_proc,pres_per_proc
      real :: average_density
      real, dimension (nx) :: r_mn,lnrho,TT,ss
      real, pointer :: gravx
      complex :: omega_jeans
      integer :: j, ierr
      logical :: lnothing
!
      intent(inout) :: f
      type (pencil_case) :: p
!
!  Define bottom and top height.
!
      zbot=xyz0(3)
      ztop=xyz0(3)+Lxyz(3)
!
!  Set default values for sound speed at top and bottom.
!  These may be updated in one of the following initialization routines.
!
      cs2top=cs20; cs2bot=cs20
!
!  Different initializations of lnrho (called from start).
!
      lnrho0      = log(rho0)
      lnrho_left  = log(rho_left)
      lnrho_right = log(rho_right)
!
      lnothing=.true.
!
      do j=1,ninit
!
        if (initlnrho(j)=='nothing') cycle
!
        lnothing=.false.
!
        iinit_str=itoa(j)
!
        select case (initlnrho(j))
!
        case ('const')
          if (lroot) print*,'initialize anelastic: const'
          do n=1,mz; do m=1,my
            if (lanelastic_lin) then
              f(1:mx,m,n,irho_b)=rho0
            else
              f(1:mx,m,n,irho)=rho0
            endif
          enddo; enddo
          f(:,:,:,ipp)=0.
!
        case ('test-poisson')
          if (lroot) print*,'initialize anelastic: test-poisson'
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ipp)=-2.*sin(x(l1:l2))*cos(z(n))
          enddo; enddo
!          call inverse_laplacian(f,f(l1:l2,m1:m2,n1:n2,ipp))
          call inverse_laplacian_z(f,f(l1:l2,m1:m2,n1:n2,ipp))
          print*, 'print results in binary file'
          write(11) f(l1:l2,4,n1:n2,ipp)
!
        case ('-ln(1+u2/2cs02)')
          f(:,:,:,ilnrho) = -alog(1. &
            +(f(:,:,:,iux)**2+f(:,:,:,iuy)**2+f(:,:,:,iuz)**2)/(2.*cs0**2))
!
        case ('anelastic')
          do m=1,my; do n=1,mz
            if (lanelastic_lin) then
              f(1:mx,m,n,ipp)=0.0
              f(1:mx,m,n,irho_b)=rho0*exp(gamma*gravz*z(n)/cs20) ! Define the base state density
            else
              f(1:mx,m,n,irho)=rho0*exp(gamma*gravz*z(n)/cs20)
              f(1:mx,m,n,ipp)=f(1:mx,m,n,irho)*cs20
            endif
          enddo; enddo
!
        case ('polytropic_simple')
          if (lanelastic_lin) then
            call polytropic_simple(f)
            f(:,:,:,ipp)=0.0
          else
            call fatal_error('init_lnrho','Not coded yet')
          endif
!
        case default
!
!  Catch unknown values
!
          write(unit=errormsg,fmt=*) 'No such value for initlnrho(' &
                            //trim(iinit_str)//'): ',trim(initlnrho(j))
          call fatal_error('init_lnrho',errormsg)
! 
        endselect
!
!  if the ipp f-array exists (e.g. in anelastic problems), set it
!  (for now corresponding to an isothermal eos)
!
!        if (ipp/=0.and.leos) f(:,:,:,ipp) = f(:,:,:,irho)*cs20
!
        if (lroot) print*,'init_lnrho: initlnrho('//trim(iinit_str)//') = ', &
            trim(initlnrho(j))
!
      enddo  ! End loop over initial conditions.
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_lnrho(f)
!
      if (lnothing.and.lroot) print*,'init_lnrho: nothing'
!
!  check that cs2bot,cs2top are ok
!  for runs with ionization or fixed ionization, don't print them
!
      if (leos_ionization .or. leos_fixed_ionization) then
        cs2top=impossible
        cs2bot=impossible
      else
        if (lroot) print*,'init_lnrho: cs2bot,cs2top=',cs2bot,cs2top
      endif
!
!  If unlogarithmic density considered, take exp of lnrho resulting from
!  initlnrho
!
!      if (ldensity_nolog) f(:,:,:,irho)=exp(f(:,:,:,ilnrho))
!
!  sanity check
!
!      if (notanumber(f(l1:l2,m1:m2,n1:n2,irho))) then
!        call error('init_rho', 'Infinit density values')
!      endif
!
    endsubroutine init_lnrho
!**********************************************************************
    subroutine calc_ldensity_pars(f)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
   endsubroutine calc_ldensity_pars
!***********************************************************************
    subroutine numerical_equilibrium(f)
!
!  sets gravity gg in order to achieve an numerical exact equilbrium
!  at t=0. This is only valid for the polytropic case, i.e.
!
!    (1/rho) grad(P) = cs20 (rho/rho0)^(gamma-2) grad(rho)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: lnrho,cs2
      real, dimension (nx,3) :: glnrho
      real, dimension (nx,3) :: gg_mn
      integer :: i,j,ilnrho
!
      do m=m1,m2
      do n=n1,n2
        lnrho=f(l1:l2,m,n,ilnrho)
        cs2=cs20*exp(gamma_m1*(lnrho-lnrho0))
        call grad(f,ilnrho,glnrho)
        do j=1,3
          gg_mn(:,j)=cs2*glnrho(:,j)
        enddo
        f(l1:l2,m,n,iglobal_gg:iglobal_gg+2)=gg_mn
      enddo
      enddo
!
    endsubroutine numerical_equilibrium
!***********************************************************************
    subroutine pencil_criteria_density()
!
!  All pencils that the Density module depends on are specified here.
!
!  19-11-04/anders: coded
!
      lpenc_requested(i_rho)=.true.
!      lpenc_requested(i_pp)=.true.
!      lpenc_requested(i_lnrho)=.true.
!      lpenc_requested(i_ugrho)=.true.
!
!      lpenc_diagnos2d(i_lnrho)=.true.
      lpenc_diagnos2d(i_rho)=.true.
!
!      if (idiag_lnrho2m/=0) lpenc_diagnos(i_lnrho)=.true.
!      if (idiag_ugrhom/=0) lpenc_diagnos(i_ugrho)=.true.
!      if (idiag_uglnrhom/=0) lpenc_diagnos(i_uglnrho)=.true.
!      if (idiag_divrhoum/=0.or.idiag_divrhourms/=0..or.idiag_divrhoumax/=0.) then
!         lpenc_diagnos(i_rho)=.true.
!!        lpenc_diagnos(i_uglnrho)=.true.
!         lpenc_diagnos(i_ugrho)=.true.
!         lpenc_diagnos(i_divu)=.true.
!      endif
!
    endsubroutine pencil_criteria_density
!***********************************************************************
    subroutine pencil_interdep_density(lpencil_in)
!
!  Interdependency among pencils from the Density module is specified here.
!
!  19-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
!      if (ldensity_nolog) then
!        if (lpencil_in(i_rho1)) lpencil_in(i_rho)=.true.
!      else
!        if (lpencil_in(i_rho)) lpencil_in(i_rho1)=.true.
!      endif
!      if (lpencil_in(i_uglnrho)) then
!        lpencil_in(i_uu)=.true.
!        lpencil_in(i_glnrho)=.true.
!      endif
!      if (lpencil_in(i_ugrho)) then
!        lpencil_in(i_uu)=.true.
!        lpencil_in(i_grho)=.true.
!      endif
!      if (lpencil_in(i_glnrho2)) lpencil_in(i_glnrho)=.true.
!      if (lpencil_in(i_sglnrho)) then
!        lpencil_in(i_sij)=.true.
!        lpencil_in(i_glnrho)=.true.
!      endif
!      if (lpencil_in(i_uij5glnrho)) then
!        lpencil_in(i_uij5)=.true.
!        lpencil_in(i_glnrho)=.true.
!      endif
!  The pencils glnrho and grho come in a bundle.
!      if (lpencil_in(i_glnrho) .and. lpencil_in(i_grho)) then
!        if (ldensity_nolog) then
!          lpencil_in(i_grho)=.false.
!        else
!          lpencil_in(i_glnrho)=.false.
!        endif
!      endif
!
    endsubroutine pencil_interdep_density
!***********************************************************************
    subroutine calc_pencils_density(f,p)
!
!  Dummy routine copied from nodensity.f90 
!
!  20-11-04/anders: coded
!
      use EquationOfState, only: lnrho0, rho0
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
      integer :: i, mm, nn, ierr,l, irhoxx
!
      if (ldensity_nolog) call fatal_error('density_anelastic','working with lnrho')
!
      if (lanelastic_lin) then
        irhoxx=irho_b
      else 
        irhoxx=irho
      endif
!
      p%rho=f(l1:l2,m,n,irhoxx)
      p%rho1=1./p%rho
      p%lnrho = log(p%rho)
! glnrho and grho
      if (lpencil(i_grho)) call grad(f, irhoxx, p%grho)
      if (lpencil(i_glnrho)) then 
        p%glnrho(:,1)=p%grho(:,1)*p%rho1
        p%glnrho(:,2)=p%grho(:,2)*p%rho1
        p%glnrho(:,3)=p%grho(:,3)*p%rho1
      endif
! uglnrho
      if (lpencil(i_uglnrho)) call dot(p%uu,p%glnrho,p%uglnrho)
! ugrho
      if (lpencil(i_ugrho)) then
        call grad(f,irhoxx,p%grho)
        call dot(p%uu,p%grho,p%ugrho)
      endif
! sglnrho
      if (lpencil(i_sglnrho)) call multmv(p%sij,p%glnrho,p%sglnrho)
!
    endsubroutine calc_pencils_density
!***********************************************************************
    subroutine density_before_boundary(f)
!
!  Actions to take before boundary conditions are set.
!
!   2-apr-08/anders: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!      if ( (.not.ldensity_nolog) .and. (irho/=0) ) &
!          f(l1:l2,m1:m2,n1:n2,irho)=exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
!
    endsubroutine density_before_boundary
!***********************************************************************
    subroutine dlnrho_dt(f,df,p)
!
! Dummy routine (taken from nodensity.f90 ) 
!  14-oct-09/dhruba: coded
!
      use sub
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,df,p
!!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dlnrho_dt: not SOLVING  dlnrho_dt in anelastic'
!      if (headtt) call identify_bcs('lnrho',ilnrho)
      if (headtt) call identify_bcs('pp',ipp)
      if (headtt) call identify_bcs('rhs',irhs)
      if (headtt) call identify_bcs('rhs',irhs+1)
      if (headtt) call identify_bcs('rhs',irhs+2)
      if (headtt) call identify_bcs('idivu',idivu)
!
!  Calculate density diagnostics
!
      if (ldiagnos) then
        if (idiag_rhom/=0)     call sum_mn_name(p%rho,idiag_rhom)
        if (idiag_totmass/=0)  call sum_mn_name(p%rho,idiag_totmass,lint=.true.)
        if (idiag_mass/=0)     call integrate_mn_name(p%rho,idiag_mass)
        if (idiag_rhomin/=0) &
            call max_mn_name(-p%rho,idiag_rhomin,lneg=.true.)
        if (idiag_rhomax/=0)   call max_mn_name(p%rho,idiag_rhomax)
        if (idiag_rho2m/=0)    call sum_mn_name(p%rho**2,idiag_rho2m)
        if (idiag_lnrho2m/=0)  call sum_mn_name(p%lnrho**2,idiag_lnrho2m)
        if (idiag_drho2m/=0)   call sum_mn_name((p%rho-rho0)**2,idiag_drho2m)
        if (idiag_drhom/=0)    call sum_mn_name(p%rho-rho0,idiag_drhom)
        if (idiag_ugrhom/=0)   call sum_mn_name(p%ugrho,idiag_ugrhom)
        if (idiag_divrhoum/=0) &
           call sum_mn_name(p%rho*p%divu+p%ugrho,idiag_divrhoum)
        if (idiag_divrhourms/=0) call sum_mn_name((p%rho*p%divu+p%ugrho)**2,idiag_divrhourms,lsqrt=.true.)
        if (idiag_divrhoumax/=0) call max_mn_name(p%rho*p%divu+p%ugrho,idiag_divrhoumax)
        if (idiag_dtd/=0) &
            call max_mn_name(diffus_diffrho/cdtv,idiag_dtd,l_dt=.true.)
      endif
!
    endsubroutine dlnrho_dt
!***********************************************************************
!  Here comes a collection of different density stratification routines
!***********************************************************************
    subroutine isothermal_density(f)
!
!  Isothermal stratification (for lnrho and ss)
!  This routine should be independent of the gravity module used.
!  When entropy is present, this module also initializes entropy.
!
!  Sound speed (and hence Temperature), and density (at infinity) are
!  initialised to their respective reference values:
!           sound speed: cs^2_0            from start.in
!           density: rho0 = exp(lnrho0)
!
!   8-jul-02/axel: incorporated/adapted from init_lnrho
!  11-jul-02/axel: fixed sign; should be tmp=gamma*pot/cs20
!  02-apr-03/tony: made entropy explicit rather than using tmp/-gamma
!  11-jun-03/tony: moved entropy initialisation to separate routine
!                  to allow isothermal condition for arbitrary density
!
      use Gravity
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: pot,tmp
      real :: cp1
!
!  Stratification depends on the gravity potential
!
      if (lroot) print*,'isothermal_density: isothermal stratification'
      if (gamma/=1.0) then
        if ((.not. lentropy) .and. (.not. ltemperature)) & 
          call fatal_error('isothermal_density','for gamma/=1.0, you need entropy or temperature!');
      endif
!
      call get_cp1(cp1)
      do n=n1,n2
        do m=m1,m2
          call potential(x(l1:l2),y(m),z(n),pot=pot)
          tmp=-gamma*pot/cs20
          f(l1:l2,m,n,ilnrho) = f(l1:l2,m,n,ilnrho) + lnrho0 + tmp
          if (lentropy) f(l1:l2,m,n,iss) = f(l1:l2,m,n,iss) &
               -gamma_m1*(f(l1:l2,m,n,ilnrho)-lnrho0)/gamma
          if (ltemperature) f(l1:l2,m,n,ilnTT)=log(cs20*cp1/gamma_m1)
        enddo
      enddo
!
!  cs2 values at top and bottom may be needed to boundary conditions.
!  The values calculated here may be revised in the entropy module.
!
      cs2bot=cs20
      cs2top=cs20
!
    endsubroutine isothermal_density
!***********************************************************************
    subroutine polytropic_simple(f)
!
!  Polytropic stratification (for lnrho and ss)
!  This routine should be independent of the gravity module used.
!
!  To maintain continuity with respect to the isothermal case,
!  one may want to specify cs20 (=1), and so zinfty is calculated from that.
!  On the other hand, for polytropic atmospheres it may be more
!  physical to specify zinfty (=1), ie the top of the polytropic atmosphere.
!  This is done if zinfty is different from 0.
!
!   8-jul-02/axel: incorporated/adapted from init_lnrho
!
      use Gravity, only: gravz_profile,gravz,zinfty,zref,zgrav,  &
                             potential,nu_epicycle
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: pot,dlncs2,r_mn
      real :: ggamma,ztop,zbot,zref2,pot_ext,lnrho_ref,ptop,pbot
!
!  identifier
!
      if (lroot) print*,'polytropic_simple: mpoly=',mpoly
!
!  The following is specific only to cases with gravity in the z direction
!  zref is calculated such that rho=rho0 and cs2=cs20 at z=zref.
!  Note: gravz is normally negative!
!
      if (lgravz) then
        if (gravz_profile=='const') then
          if (lroot.and.gravz==0.) print*,'polytropic_simple: divide by gravz=0'
          zref=zinfty-(mpoly+1.)*cs20/(-gamma*gravz)
        elseif (gravz_profile=='const_zero') then
          if (lroot.and.gravz==0.) print*,'polytropic_simple: divide by gravz=0'
          zref=zinfty-(mpoly+1.)*cs20/(-gamma*gravz)
        elseif (gravz_profile=='linear') then
          if (lroot.and.gravz==0.) print*,'polytropic_simple: divide by gravz=0'
          zref2=zinfty**2-(mpoly+1.)*cs20/(0.5*gamma*nu_epicycle**2)
          if (zref2<0) then
            if (lroot) print*,'polytropic_simple: zref**2<0 is not ok'
            zref2=0. !(and see what happens)
          endif
          zref=sqrt(zref2)
        else
          if (lroot) print*,'polytropic_simple: zref not prepared!'
        endif
        if (lroot) print*,'polytropic_simple: zref=',zref
!
!  check whether zinfty lies outside the domain (otherwise density
!  would vanish within the domain). At the moment we are not properly
!  testing the lower boundary on the case of a disc (commented out).
!
        ztop=xyz0(3)+Lxyz(3)
        zbot=xyz0(3)
        !-- if (zinfty<min(ztop,zgrav) .or. (-zinfty)>min(zbot,zgrav)) then
        if (zinfty<min(ztop,zgrav)) then
          if (lroot) print*,'polytropic_simple: domain too big; zinfty=',zinfty
          !call stop_it( &
          !         'polytropic_simply: rho and cs2 will vanish within domain')
        endif
      endif
!
!  stratification Gamma (upper case in the manual)
!
      ggamma=1.+1./mpoly
!
!  polytropic sphere with isothermal exterior
!  calculate potential at the stellar surface, pot_ext
!
      if (lgravr) then
        call potential(R=r_ext,POT=pot_ext)
        cs2top=-gamma/(mpoly+1.)*pot_ext
        lnrho_ref=mpoly*log(cs2top)-(mpoly+1.)
        print*,'polytropic_simple: pot_ext=',pot_ext
        do n=n1,n2; do m=m1,m2
          r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
          call potential(x(l1:l2),y(m),z(n),pot=pot)
!
!  density
!  these formulae assume lnrho0=0 and cs0=1
!
          where (r_mn > r_ext)
            f(l1:l2,m,n,irho_b)=exp(lnrho_ref-gamma*pot/cs2top)
          elsewhere
            dlncs2=log(-gamma*pot/((mpoly+1.)*cs20))
            f(l1:l2,m,n,irho_b)=exp(lnrho0+mpoly*dlncs2)
          endwhere
!
!  entropy
!
          if (lentropy) then
            where (r_mn > r_ext)
              f(l1:l2,m,n,iss_b)=-(1.-1./gamma)*log(f(l1:l2,m,n,irho_b))+log(cs2top)/gamma
            elsewhere
              dlncs2=log(-gamma*pot/((mpoly+1.)*cs20))
              f(l1:l2,m,n,iss_b)=mpoly*(ggamma/gamma-1.)*dlncs2
            endwhere
          endif
        enddo; enddo
      else
!
!  cartesian case with gravity in the z direction
!
        do n=n1,n2; do m=m1,m2
          call potential(x(l1:l2),y(m),z(n),pot=pot)
          dlncs2=log(-gamma*pot/((mpoly+1.)*cs20))
          f(l1:l2,m,n,irho_b)=rho0*exp(mpoly*dlncs2)
!          if (lentropy) f(l1:l2,m,n,iss_b)=mpoly*(ggamma/gamma-1.)*dlncs2
          if (lentropy) f(l1:l2,m,n,iss_b)=(ggamma/gamma-1.)*log(f(l1:l2,m,n,irho_b)/rho0)
        enddo; enddo
              do n=n1,n2
              write(*,*) f(l1,m1,n,irho_b)
              end do
!
!  cs2 values at top and bottom may be needed to boundary conditions.
!  In spherical geometry, ztop is z at the outer edge of the box,
!  so this calculation still makes sense.
!
        call potential(xyz0(1),xyz0(2),ztop,pot=ptop)
        cs2top=-gamma*ptop/(mpoly+1.)
!
!  In spherical geometry ztop should never be used.
!  Even in slab geometry ztop is not normally used.
!
        call potential(xyz0(1),xyz0(2),zbot,pot=pbot)
        cs2bot=-gamma*pbot/(mpoly+1.)
      endif
!
    endsubroutine polytropic_simple
!***********************************************************************
    subroutine read_density_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=density_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=density_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_density_init_pars
!***********************************************************************
    subroutine write_density_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=density_init_pars)
!
    endsubroutine write_density_init_pars
!***********************************************************************
    subroutine read_density_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=density_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=density_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_density_run_pars
!***********************************************************************
    subroutine write_density_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=density_run_pars)
!
    endsubroutine write_density_run_pars
!***********************************************************************
    subroutine rprint_density(lreset,lwrite)
!
!  reads and registers print parameters relevant for compressible part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname, inamex, inamey, inamez, inamexy, inamexz, irz, inamer
      logical :: lwr
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_rhom=0; idiag_rho2m=0; idiag_lnrho2m=0
        idiag_drho2m=0; idiag_drhom=0
        idiag_ugrhom=0; idiag_uglnrhom=0
        idiag_rhomin=0; idiag_rhomax=0; idiag_dtd=0
        idiag_lnrhomphi=0; idiag_rhomphi=0
        idiag_rhomz=0; idiag_rhomy=0; idiag_rhomx=0 
        idiag_rhomxy=0; idiag_rhomr=0; idiag_totmass=0
        idiag_rhomxz=0; idiag_divrhoum=0; idiag_divrhourms=0; idiag_divrhoumax=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) print*,'rprint_density: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'rhom',idiag_rhom)
        call parse_name(iname,cname(iname),cform(iname),'rho2m',idiag_rho2m)
        call parse_name(iname,cname(iname),cform(iname),'drho2m',idiag_drho2m)
        call parse_name(iname,cname(iname),cform(iname),'drhom',idiag_drhom)
        call parse_name(iname,cname(iname),cform(iname),'rhomin',idiag_rhomin)
        call parse_name(iname,cname(iname),cform(iname),'rhomax',idiag_rhomax)
        call parse_name(iname,cname(iname),cform(iname),'lnrho2m',idiag_lnrho2m)
        call parse_name(iname,cname(iname),cform(iname),'ugrhom',idiag_ugrhom)
        call parse_name(iname,cname(iname),cform(iname),'uglnrhom',idiag_uglnrhom)
        call parse_name(iname,cname(iname),cform(iname),'dtd',idiag_dtd)
        call parse_name(iname,cname(iname),cform(iname),'totmass',idiag_totmass)
        call parse_name(iname,cname(iname),cform(iname),'mass',idiag_mass)
        call parse_name(iname,cname(iname),cform(iname),'divrhoum',idiag_divrhoum)
        call parse_name(iname,cname(iname),cform(iname),'divrhourms',idiag_divrhourms)
        call parse_name(iname,cname(iname),cform(iname),'divrhoumax',idiag_divrhoumax)
!
!  alternatively, use these shorter names: drurms and drumax,
!  instead of divrhourms and divrhoumax
!
        call parse_name(iname,cname(iname),cform(iname),'drurms',idiag_divrhourms)
        call parse_name(iname,cname(iname),cform(iname),'drumax',idiag_divrhoumax)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rhomz',idiag_rhomz)
      enddo
!
!  check for those quantities for which we want xz-averages
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'rhomy',idiag_rhomy)
      enddo
!
!  check for those quantities for which we want yz-averages
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'rhomx',idiag_rhomx)
      enddo
!
!  check for those quantities for which we want phiz-averages
!
      do inamer=1,nnamer
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'rhomr',idiag_rhomr)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do inamexz=1,nnamexz
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'rhomxz',idiag_rhomxz)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'rhomxy',idiag_rhomxy)
      enddo
!
!  check for those quantities for which we want phi-averages
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),&
            'lnrhomphi',idiag_lnrhomphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'rhomphi',idiag_rhomphi)
      enddo
!
!  write column where which density variable is stored
!
      if (lwr) then
        write(3,*) 'i_rhom=',idiag_rhom
        write(3,*) 'i_rho2m=',idiag_rho2m
        write(3,*) 'i_drho2m=',idiag_drho2m
        write(3,*) 'i_drhom=',idiag_drhom
        write(3,*) 'i_rhomin=',idiag_rhomin
        write(3,*) 'i_rhomax=',idiag_rhomax
        write(3,*) 'i_lnrho2m=',idiag_lnrho2m
        write(3,*) 'i_ugrhom=',idiag_ugrhom
        write(3,*) 'i_uglnrhom=',idiag_uglnrhom
        write(3,*) 'i_rhomz=',idiag_rhomz
        write(3,*) 'i_rhomy=',idiag_rhomy
        write(3,*) 'i_rhomx=',idiag_rhomx
        write(3,*) 'i_rhomxy=',idiag_rhomxy
        write(3,*) 'i_rhomxz=',idiag_rhomxz
        write(3,*) 'nname=',nname
        write(3,*) 'ilnrho=',ilnrho
        write(3,*) 'irho=',irho
        write(3,*) 'i_lnrhomphi=',idiag_lnrhomphi
        write(3,*) 'i_rhomphi=',idiag_rhomphi
        write(3,*) 'i_rhomr=',idiag_rhomr
        write(3,*) 'i_dtd=',idiag_dtd
        write(3,*) 'i_totmass=',idiag_totmass
        write(3,*) 'i_mass=',idiag_mass
        write(3,*) 'i_divrhoum=',idiag_divrhoum
        write(3,*) 'i_divrhourms=',idiag_divrhourms
        write(3,*) 'i_divrhoumax=',idiag_divrhoumax
      endif
!
    endsubroutine rprint_density
!***********************************************************************
    subroutine get_init_average_density(f,init_average_density)
!
!  10-dec-09/piyali: added to pass initial average density 
!  equ.f90 
!
    use Diagnostics, only: integrate_mn,get_average_density
!
    real, dimension (mx,my,mz,mfarray):: f
    real, dimension(1) :: mass_per_proc
    real :: init_average_density
    intent(in):: f
!
    do m=m1,m2
    do n=n1,n2
      call integrate_mn(exp(f(l1:l2,m,n,ilnrho)),mass_per_proc(1))
    enddo
    enddo
    call get_average_density(mass_per_proc(1),init_average_density)
!
    endsubroutine get_init_average_density
!***********************************************************************
    subroutine get_slices_density(f,slices)
!
!  Write slices for animation of Density variables.
!
!  26-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      if (.not. lanelastic_lin) then
!  Loop over slices
      
        select case (trim(slices%name))

!  Density.

          case ('rho')
              slices%yz =f(ix_loc,m1:m2,n1:n2,irho)
              slices%xz =f(l1:l2,iy_loc,n1:n2,irho)
              slices%xy =f(l1:l2,m1:m2,iz_loc,irho)
              slices%xy2=f(l1:l2,m1:m2,iz2_loc,irho)
              if (lwrite_slice_xy3) slices%xy3=f(l1:l2,m1:m2,iz3_loc,irho)
              if (lwrite_slice_xy4) slices%xy4=f(l1:l2,m1:m2,iz4_loc,irho)
              slices%ready=.true.
!
!  Logarithmic density.
!
          case ('lnrho')
            call fatal_error('get_slices_density','Not working with lnrho anymore')
        endselect
      endif
!
    endsubroutine get_slices_density
!***********************************************************************
    subroutine get_slices_pressure(f,slices)
!
!  Write slices for animation of Pressure variables.
!
!  26-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Pressure.
!
        case ('pp')
            slices%yz =f(ix_loc,m1:m2,n1:n2,ipp)
            slices%xz =f(l1:l2,iy_loc,n1:n2,ipp)
            slices%xy =f(l1:l2,m1:m2,iz_loc,ipp)
            slices%xy2=f(l1:l2,m1:m2,iz2_loc,ipp)
            if (lwrite_slice_xy3) slices%xy3=f(l1:l2,m1:m2,iz3_loc,ipp)
            if (lwrite_slice_xy4) slices%xy4=f(l1:l2,m1:m2,iz4_loc,ipp)
            slices%ready=.true.
!
      endselect
!
    endsubroutine get_slices_pressure
!***********************************************************************
    subroutine anelastic_after_mn(f, p, df, mass_per_proc)
!
      use Poisson, only: inverse_laplacian
      use Mpicomm, only: initiate_isendrcv_bdry, finalize_isendrcv_bdry
      use Boundcond, only: boundconds
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
!      real, dimension (nx,ny,nz) :: pold
      type (pencil_case) :: p
      real, dimension (nx,3) :: gpp
      real, dimension (nx) :: phi_rhs_pencil
      real, dimension (1)  :: mass_per_proc
      real    :: average_density,average_pressure,init_average_density
      integer :: j, ju, l, i
!
!  Set first the boundary conditions on rhs
!
      call initiate_isendrcv_bdry(f,irhs,irhs+2)
      call finalize_isendrcv_bdry(f,irhs,irhs+2)
      call boundconds(f,irhs,irhs+2)
!
!  Find the divergence of rhs
!
!      pold(1:nx,1:ny,1:nz)=f(l1:l2,m1:m2,n1:n2,ipp)
      do n=n1,n2
      do m=m1,m2
        call div(f,irhs,phi_rhs_pencil)
        f(l1:l2,m,n,ipp)=phi_rhs_pencil
      enddo
      enddo
!
!  get pressure from inverting the Laplacian
!
      if (lperi(3)) then
        call inverse_laplacian(f,f(l1:l2,m1:m2,n1:n2,ipp))
!  refresh the ghost zones: periodic pressure
        call initiate_isendrcv_bdry(f,ipp)
        call finalize_isendrcv_bdry(f,ipp)
        call boundconds(f,ipp,ipp)
        if (.not. lanelastic_lin) then 
          call get_average_density(mass_per_proc(1),average_density)
          call get_average_pressure(init_average_density,average_density,average_pressure)
          f(:,:,:,ipp)=f(:,:,:,ipp)+average_pressure
        endif
      else
        call inverse_laplacian_z(f,f(l1:l2,m1:m2,n1:n2,ipp))       
!
! fill the ghost zones: 
! In the vertical direction: dP/dz=0 (symmetric)
!
        m=4; n=n1
        do i=1,nghost
           f(l1:l2,4,n1-i,ipp)= f(l1:l2,4,n1+i,ipp)
        enddo
        n=n2
        do i=1,nghost
          f(l1:l2,4,n2+i,ipp)=f(l1:l2,4,n2-i,ipp)
        enddo
! Bc in the horizontal direction: periodic
        f(1:l1-1,:,:,ipp) = f(l2i:l2,:,:,ipp)
        f(l2+1:mx,:,:,ipp) = f(l1:l1i,:,:,ipp)
      endif
!
!  Add the pressure gradient term to the NS equation
!
      do n=n1,n2; do m=m1,m2
        call grad(f,ipp,gpp)
        do j=1,3
          ju=j+iuu-1
          if (lanelastic_lin) then
!            df(l1:l2,m,n,ju)=df(l1:l2,m,n,ju)-gpp(:,j)/f(l1:l2,m,n,irho_b) &
!                + gamma*f(l1:l2,m,n,ipp)/(f(l1:l2,m,n,irho_b)*p%cs2)*p%gg(:,j)
            df(l1:l2,m,n,ju)=df(l1:l2,m,n,ju)-gpp(:,j)/f(l1:l2,m,n,irho_b)
          else
            df(l1:l2,m,n,ju)=df(l1:l2,m,n,ju)-gpp(:,j)/f(l1:l2,m,n,irho)
          endif 
        enddo
!
        if (.not. lanelastic_lin) then
          f(l1:l2,m,n,irho)=f(l1:l2,m,n,ipp)*cs20
!         call initiate_isendrcv_bdry(f,ilnrho)
!         call finalize_isendrcv_bdry(f,ilnrho)
!         call boundconds(f,ilnrho,ilnrho)
        endif
      enddo; enddo
!
    endsubroutine anelastic_after_mn
!***********************************************************************
    subroutine inverse_laplacian_z(f,phi)
!
!  Solve the Poisson equation by Fourier transforming in the xy-plane and
!  solving the discrete matrix equation in the z-direction.
!
!  02-mar-2012/dintrans+dubuffet: coded
!
      use Fourier, only: fourier_transform_xy
      use General, only: tridag
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny,nz) :: phi, b1, d2uz, b2
      real, dimension (nx) :: tmpx
      real, dimension (nz) :: a_tri, b_tri, c_tri, r_tri, u_tri
      real :: k2
      integer :: ikx, iky, j, k
      logical :: err
!
!  Identify version.
!
      if (lroot .and. ip<10) call svn_id( &
          '$Id$')
!
! test with the pressure
!
!      do n=n1,n2
!      do m=m1,m2
!!        call der2(f,iuz,tmpx,3)
!        call der(f,iuz,tmpx,3)
!        d2uz(:,m-nghost,n-nghost)=1.e-3*tmpx
!      enddo
!      enddo
!
!  The right-hand-side of the Poisson equation is purely real.
!
      b1 = 0.0
      b2 = 0.0
!
!  Forward transform (to k-space).
!
      call fourier_transform_xy(phi,b1)
!      call fourier_transform_xy(d2uz,b2)
!
!  Solve for discrete z-direction with zero density above and below z-boundary.
!
      do iky=1,ny
        do ikx=1,nx
          if ((kx_fft(ikx)==0.0) .and. (ky_fft(iky)==0.0)) then
            phi(ikx,iky,:) = 0.0
          else
            k2=kx_fft(ikx)**2+ky_fft(iky)**2
            a_tri=1.0/dz**2
            b_tri=-2.0/dz**2-k2
            c_tri=1.0/dz**2
            r_tri=phi(ikx,iky,:)
! P = 0 (useful to test the Poisson solver)
!            b_tri(1)=1.  ; c_tri(1)=0.  ; r_tri(1)=0.
!            b_tri(nz)=1. ; a_tri(nz)=0. ; r_tri(nz)=0.
! dP/dz = 0
            c_tri(1)=c_tri(1)+a_tri(1)
            a_tri(nz)=a_tri(nz)+c_tri(nz)
! delta P = 0 outside
!            b_tri(1)=b_tri(1)+2./dz**2/(2.+2*sqrt(k2)*dz+k2*dz**2)
!            b_tri(nz)=b_tri(nz)+2./dz**2/(2.+2*sqrt(k2)*dz+k2*dz**2)
! dP/dz = mu d2 uz /dz2
!            c_tri(1)=c_tri(1)+a_tri(1)
!            r_tri(1)=r_tri(1)+2./dz*d2uz(ikx,iky,1)
!            a_tri(nz)=a_tri(nz)+c_tri(nz)
!            r_tri(nz)=r_tri(nz)-2./dz*d2uz(ikx,iky,nz)
! P = mu d uz /dz
!            b_tri(1)=1.
!            c_tri(1)=0.
!            r_tri(1)=d2uz(ikx,iky,1)
!            b_tri(nz)=1.
!            a_tri(nz)=0.
!            r_tri(nz)=d2uz(ikx,iky,nz)
!
            call tridag(a_tri,b_tri,c_tri,r_tri,u_tri,err)
            phi(ikx,iky,:)=u_tri
          endif
        enddo
      enddo
!
      do iky=1,ny
        do ikx=1,nx
          if ((kx_fft(ikx)==0.0) .and. (ky_fft(iky)==0.0)) then
            b1(ikx,iky,:) = 0.0
          else
            k2=kx_fft(ikx)**2+ky_fft(iky)**2
            a_tri=1.0/dz**2
            b_tri=-2.0/dz**2-k2
            c_tri=1.0/dz**2
            r_tri=b1(ikx,iky,:)
! P = 0 (useful to test the Poisson solver)
!            b_tri(1)=1.  ; c_tri(1)=0.  ; r_tri(1)=0.
!            b_tri(nz)=1. ; a_tri(nz)=0. ; r_tri(nz)=0.
! dP/dz = 0
            c_tri(1)=c_tri(1)+a_tri(1)
            a_tri(nz)=a_tri(nz)+c_tri(nz)
! delta P = 0 outside
!            b_tri(1)=b_tri(1)+2./dz**2/(2.+2*sqrt(k2)*dz+k2*dz**2)
!            b_tri(nz)=b_tri(nz)+2./dz**2/(2.+2*sqrt(k2)*dz+k2*dz**2)
! dP/dz = mu d2 uz /dz2
!            c_tri(1)=c_tri(1)+a_tri(1)
!            r_tri(1)=r_tri(1)+2./dz*b2(ikx,iky,1)
!            a_tri(nz)=a_tri(nz)+c_tri(nz)
!            r_tri(nz)=r_tri(nz)-2./dz*b2(ikx,iky,nz)
! P = mu d uz /dz
!            b_tri(1)=1.
!            c_tri(1)=0.
!            r_tri(1)=b2(ikx,iky,1)
!            b_tri(nz)=1.
!            a_tri(nz)=0.
!            r_tri(nz)=b2(ikx,iky,nz)
!
            call tridag(a_tri,b_tri,c_tri,r_tri,u_tri,err)
            b1(ikx,iky,:)=u_tri
          endif
        enddo
      enddo
!
!  Inverse transform (to real space).
!
      call fourier_transform_xy(phi,b1,linv=.true.)
!
    endsubroutine inverse_laplacian_z
!***********************************************************************
    subroutine dynamical_diffusion(umax)
!
!  Dynamically set mass diffusion coefficient given fixed mesh Reynolds number.
!
!  27-jul-11/nils: coded
!
      real, intent(in) :: umax
!
      call keep_compiler_quiet(umax)
      call fatal_error('dynamical_diffusion',&
          'This subroutine is not yet implemented for the anelastic module.')
!
    endsubroutine dynamical_diffusion
!***********************************************************************
    subroutine impose_density_floor(f)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine impose_density_floor
!***********************************************************************
endmodule Density
