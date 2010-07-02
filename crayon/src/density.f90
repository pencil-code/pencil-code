! $Id: density.f90 13579 2010-03-31 11:58:14Z AxelBrandenburg $
!
!  This module takes care of the continuity equation.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ldensity = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED lnrho; rho; rho1; glnrho(3); grho(3); uglnrho; ugrho
! PENCILS PROVIDED glnrho2; del2lnrho; del2rho 
! PENCILS PROVIDED hlnrho(3,3); sglnrho(3) 
!
!***************************************************************
module Density
!
  use Cparam
  use Cdata
  use Messages
  use EquationOfState
  use Sub, only : keep_compiler_quiet
!
  implicit none
!
  include 'density.h'
!
  real, dimension (ninit) :: ampllnrho=0.0, widthlnrho=0.1
  real, dimension (ninit) :: rho_left=1.0, rho_right=1.0
  real, dimension (ninit) :: amplrho=0.0, phase_lnrho=0.0, radius_lnrho=0.5
  real, dimension (ninit) :: kx_lnrho=1.0, ky_lnrho=1.0, kz_lnrho=1.0
  real, dimension (ninit) :: kxx_lnrho=0.0, kyy_lnrho=0.0, kzz_lnrho=0.0
  real, dimension (nz,3) :: glnrhomz
  real, dimension (mz) :: lnrho_init_z=0.0, del2lnrho_init_z=0.0
  real, dimension (mz) :: dlnrhodz_init_z=0.0, glnrho2_init_z=0.0
  real, dimension (3) :: diffrho_hyper3_aniso=0.0
  real :: lnrho_const=0.0, rho_const=1.0
  real :: cdiffrho=0.0, diffrho=0.0
  real :: diffrho_hyper3=0.0, diffrho_hyper3_mesh=5.0, diffrho_shock=0.0
  real :: eps_planet=0.5, q_ell=5.0, hh0=0.0
  real :: xblob=0.0, yblob=0.0, zblob=0.0
  real :: co1_ss=0.0, co2_ss=0.0, Sigma1=150.0
  real :: lnrho_int=0.0, lnrho_ext=0.0, damplnrho_int=0.0, damplnrho_ext=0.0
  real :: wdamp=0.0, density_floor=-1.0
  real :: mass_source_Mdot=0.0, mass_source_sigma=0.0
  real :: radial_percent_smooth=10.0, rshift=0.0
  real, target :: plaw=0.0
  real :: lnrho_z_shift=0.0
  real :: powerlr=3.0, zoverh=1.5, hoverr=0.05
  complex :: coeflnrho=0.0
  integer, parameter :: ndiff_max=4
  integer :: iglobal_gg=0
  logical :: lmass_source=.false.,lcontinuity_gas=.true.
  logical :: ldiff_normal=.false.,ldiff_hyper3=.false.,ldiff_shock=.false.
  logical :: ldiff_hyper3lnrho=.false.,ldiff_hyper3_aniso=.false.
  logical :: ldiff_hyper3_polar=.false.,lanti_shockdiffusion=.false.
  logical :: ldiff_hyper3_mesh=.false.
  logical :: lfreeze_lnrhoint=.false.,lfreeze_lnrhoext=.false.
  logical :: lfreeze_lnrhosqu=.false.,lexponential_smooth=.false.
  logical :: lrho_as_aux=.false., ldiffusion_nolog=.false.
  logical :: lshare_plaw=.false.,lmassdiff_fix=.false.
  logical :: lcheck_negative_density=.false.
  logical :: lcalc_glnrhomean=.false.
  logical :: ldensity_profile_masscons=.false.
  character (len=labellen), dimension(ninit) :: initlnrho='nothing'
  character (len=labellen) :: strati_type='lnrho_ss'
  character (len=labellen), dimension(ndiff_max) :: idiff=''
  character (len=labellen) :: borderlnrho='nothing'
  character (len=labellen) :: mass_source_profile='cylindric'
  character (len=5) :: iinit_str
!
  namelist /density_init_pars/ &
      ampllnrho, initlnrho, widthlnrho, rho_left, rho_right, lnrho_const, &
      rho_const, cs2bot, cs2top, radius_lnrho, eps_planet, xblob, yblob, &
      zblob, b_ell, q_ell, hh0, rbound, mpoly, &
      strati_type, beta_glnrho_global, radial_percent_smooth, kx_lnrho, &
      ky_lnrho, kz_lnrho, amplrho, phase_lnrho, coeflnrho, kxx_lnrho, &
      kyy_lnrho,  kzz_lnrho, co1_ss, co2_ss, Sigma1, idiff, ldensity_nolog, &
      lexponential_smooth, wdamp, lcontinuity_gas, density_floor, &
      rshift, lrho_as_aux, ldiffusion_nolog, &
      lnrho_z_shift, powerlr,  zoverh,  hoverr
!
  namelist /density_run_pars/ &
      cdiffrho, diffrho, diffrho_shock, &
      cs2bot, cs2top, &
      idiff, lnrho_int, lnrho_ext, &
      damplnrho_int, damplnrho_ext, wdamp, &
      lnrho_const, plaw, lcontinuity_gas, &
      lfreeze_lnrhosqu, density_floor, lrho_as_aux, &
      ldiffusion_nolog, lcheck_negative_density, lmassdiff_fix, &
      lcalc_glnrhomean, ldensity_profile_masscons
!
!  Diagnostic variables (need to be consistent with reset list below).
!
  integer :: idiag_rhom=0       ! DIAG_DOC: $\left<\varrho\right>$
                                ! DIAG_DOC:   \quad(mean density)
  integer :: idiag_rho2m=0      ! DIAG_DOC:
  integer :: idiag_lnrho2m=0    ! DIAG_DOC:
  integer :: idiag_drho2m=0     ! DIAG_DOC:
  integer :: idiag_drhom=0      ! DIAG_DOC:
  integer :: idiag_rhomin=0     ! DIAG_DOC: $\min(\rho)$
  integer :: idiag_rhomax=0     ! DIAG_DOC: $\max(\rho)$
  integer :: idiag_ugrhom=0     ! DIAG_DOC: $\left<\uv\cdot\nabla\varrho\right>$
  integer :: idiag_uglnrhom=0   ! DIAG_DOC:
  integer :: idiag_lnrhomphi=0  ! PHIAVG_DOC: $\left<\ln\varrho\right>_\varphi$
  integer :: idiag_rhomphi=0    ! PHIAVG_DOC: $\left<\varrho\right>_\varphi$
  integer :: idiag_dtd=0        ! DIAG_DOC:
  integer :: idiag_rhomz=0      ! DIAG_DOC: $\left<\varrho\right>_{xy}$
  integer :: idiag_rhomy=0      ! DIAG_DOC:
  integer :: idiag_rhomx=0      ! DIAG_DOC:
  integer :: idiag_rhomxy=0     ! DIAG_DOC:
  integer :: idiag_rhomxz=0     ! DIAG_DOC:
  integer :: idiag_rhomr=0      ! DIAG_DOC:
  integer :: idiag_totmass=0    ! DIAG_DOC: $\int\varrho\,dV$
  integer :: idiag_mass=0       ! DIAG_DOC: $\int\varrho\,dV$
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
      use Sub
      use FArrayManager
!
      call farray_register_pde('lnrho',ilnrho)
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
          "$Id: density.f90 13579 2010-03-31 11:58:14Z AxelBrandenburg $")
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
      use FArrayManager
      use Mpicomm
      use SharedVariables
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      integer :: i
      logical :: lnothing
!
!  Set irho equal to ilnrho if we are considering non-logarithmic density.
!  Also reset the corresponding slot in varname.
!
      if (ldensity_nolog) then
        irho=ilnrho
        varname(irho)='rho'
      endif
!
!  initialize cs2cool to cs20
!  (currently disabled, because it causes problems with mdarf auto-test)
!     cs2cool=cs20
!
!
      if (diffrho==0.) then
!
!  Made to work by adding diffrho + cdiffrho to the rprint reset list.
!
        diffrho=cdiffrho*dxmin*cs0
      endif
!
!  Turn off continuity equation term for 0-D runs.
!
      if (nxgrid*nygrid*nzgrid==1) then
        lcontinuity_gas=.false.
        print*, 'initialize_density: 0-D run, turned off continuity equation'
      endif
!
!  Initialize mass diffusion
!
      ldiff_normal=.false.
      ldiff_shock=.false.
!
      lnothing=.false.
!
!  different choices of mass diffusion (if any)
!
      do i=1,ndiff_max
        select case (idiff(i))
        case ('normal')
          if (lroot) print*,'diffusion: div(D*grad(rho))'
          ldiff_normal=.true.
        case ('shock','diff-shock','diffrho-shock')
          if (lroot) print*,'diffusion: shock diffusion'
          ldiff_shock=.true.
        case ('','none')
          if (lroot .and. (.not. lnothing)) print*,'diffusion: nothing'
        case default
          write(unit=errormsg,fmt=*) 'initialize_density: ', &
              'No such value for idiff(',i,'): ', trim(idiff(i))
          call fatal_error('initialize_density',errormsg)
        endselect
        lnothing=.true.
      enddo
!
!  If we're timestepping, die or warn if the the diffusion coefficient that
!  corresponds to the chosen diffusion type is not set.
!
      if (lrun) then
        if (ldiff_normal.and.diffrho==0.0) &
            call warning('initialize_density', &
            'Diffusion coefficient diffrho is zero!')
        if (ldiff_shock.and.diffrho_shock==0.0) &
            call fatal_error('initialize_density', &
            'diffusion coefficient diffrho_shock is zero!')
      endif
!
! Tell the equation of state that we're here and what f variable we use
!
      if (ldensity_nolog) then
        call select_eos_variable('rho',irho)
      else
        call select_eos_variable('lnrho',ilnrho)
      endif
!
! Do not allow inconsistency between rho0 (from eos) and rho_const
! or lnrho0 and lnrho_const.
!
      if (rho0.ne.rho_const) then
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
!  Possible to store non log rho as auxiliary variable.
!
      if (lrho_as_aux) then
        if (ldensity_nolog) then
          if (lroot) print*, 'initialize_density: makes no sense to have '// &
              'lrho_as_aux=T if already evolving non log rho'
          call fatal_error('initialize_density','')
        else
          call farray_register_auxiliary('rho',irho,communicated=.true.)
        endif
      endif
!
!  For diffusion term with non-logarithmic density we need to save rho
!  as an auxiliary variable.
!
      if (ldiffusion_nolog .and. .not. lrho_as_aux) then
        if (lroot) then
          print*, 'initialize_density: must have lrho_as_aux=T '// &
              'for non-logarithmic diffusion'
          print*, '  (consider setting lrho_as_aux=T and'
          print*, '   !  MAUX CONTRIBUTION 1'
          print*, '   !  COMMUNICATED AUXILIARIES 1'
          print*, '   in cparam.local)'
        endif
        call fatal_error('initialize_density','')
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
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
      use General, only: chn,complex_phase
      use Gravity, only: zref,z1,z2,gravz,nu_epicycle,potential
      use Initcond
      use IO
      use Mpicomm
      use Sub
      use InitialCondition, only: initial_condition_lnrho
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: prof
      real, dimension (ninit) :: lnrho_left,lnrho_right
      real :: zbot,ztop
      integer :: j
      logical :: lnothing
!
      intent(inout) :: f
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
        call chn(j,iinit_str)
!
        select case (initlnrho(j))
!
        case ('zero', '0'); f(:,:,:,ilnrho)=0.
        case ('const_lnrho'); f(:,:,:,ilnrho)=lnrho_const
        case ('const_rho'); f(:,:,:,ilnrho)=log(rho_const)
        case ('constant'); f(:,:,:,ilnrho)=log(rho_left(j))
        case ('mode')
          call modes(ampllnrho(j),coeflnrho,f,ilnrho,kx_lnrho(j), &
              ky_lnrho(j),kz_lnrho(j))
        case ('blob')
          call blob(ampllnrho(j),f,ilnrho,radius_lnrho(j),xblob,yblob,zblob)
        case ('blob_hs')
          print*, 'init_lnrho: put a blob in hydrostatic equilibrium:'// &
          'radius_lnrho, ampllnrho, position=',radius_lnrho(j), &
          ampllnrho(j), xblob, yblob, zblob
          call blob(ampllnrho(j),f,ilnrho,radius_lnrho(j),xblob,yblob,zblob)
          call blob(-ampllnrho(j),f,iss,radius_lnrho(j),xblob,yblob,zblob)
        case ('hydrostatic-z', '1')
          print*, 'init_lnrho: use polytropic_simple instead!'
        case ('xjump')
          call jump(f,ilnrho,lnrho_left(j),lnrho_right(j),widthlnrho(j),'x')
        case ('yjump')
          call jump(f,ilnrho,lnrho_left(j),lnrho_right(j),widthlnrho(j),'y')
        case ('zjump')
          call jump(f,ilnrho,lnrho_left(j),lnrho_right(j),widthlnrho(j),'z')
        case ('soundwave-x')
          call soundwave(ampllnrho(j),f,ilnrho,kx=kx_lnrho(j))
        case ('soundwave-y')
          call soundwave(ampllnrho(j),f,ilnrho,ky=ky_lnrho(j))
        case ('soundwave-z')
          call soundwave(ampllnrho(j),f,ilnrho,kz=kz_lnrho(j))
        case ('sinwave-phase')
          call sinwave_phase(f,ilnrho,ampllnrho(j),kx_lnrho(j), &
              ky_lnrho(j),kz_lnrho(j),phase_lnrho(j))
        case ('sinwave-phase-nolog')
          do m=m1,m2; do n=n1,n2
            f(l1:l2,m,n,ilnrho) = f(l1:l2,m,n,ilnrho) + &
                alog(1+amplrho(j)*sin(kx_lnrho(j)*x(l1:l2)+ &
                ky_lnrho(j)*y(m)+kz_lnrho(j)*z(n)+phase_lnrho(j)))
          enddo; enddo
        case ('coswave-phase')
          call coswave_phase(f,ilnrho,ampllnrho(j),kx_lnrho(j), &
              ky_lnrho(j),kz_lnrho(j),phase_lnrho(j))
        case ('sinwave-x')
          call sinwave(ampllnrho(j),f,ilnrho,kx=kx_lnrho(j))
        case ('sinwave-y')
          call sinwave(ampllnrho(j),f,ilnrho,ky=ky_lnrho(j))
        case ('sinwave-z')
          call sinwave(ampllnrho(j),f,ilnrho,kz=kz_lnrho(j))
        case ('coswave-x')
          call coswave(ampllnrho(j),f,ilnrho,kx=kx_lnrho(j))
        case ('coswave-y')
          call coswave(ampllnrho(j),f,ilnrho,ky=ky_lnrho(j))
        case ('coswave-z')
          call coswave(ampllnrho(j),f,ilnrho,kz=kz_lnrho(j))
        case ('triquad')
          call triquad(ampllnrho(j),f,ilnrho,kx_lnrho(j), &
              ky_lnrho(j),kz_lnrho(j), kxx_lnrho(j), kyy_lnrho(j), &
              kzz_lnrho(j))
        case ('sinx_siny_sinz')
          call sinx_siny_sinz(ampllnrho(j),f,ilnrho, &
              kx_lnrho(j),ky_lnrho(j),kz_lnrho(j))
        case ('gaussian3d')
          call gaussian3d(ampllnrho(j),f,ilnrho,radius_lnrho(j))
        case ('gaussian-z')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnrho) = f(l1:l2,m,n,ilnrho) - &
                z(n)**2/(2*radius_lnrho(j)**2)
          enddo; enddo
        case ('gauss-z-offset')
          do n=n1,n2
             f(:,:,n,ilnrho) = f(:,:,n,ilnrho) + &
                alog(exp(f(:,:,n,ilnrho))+ &
                ampllnrho(j)*(exp(-(z(n)+lnrho_z_shift)**2/ &
                (2*radius_lnrho(j)**2))))
          enddo
        case ('gaussian-noise')
          If (lnrho_left(j) /= 0.) f(:,:,:,ilnrho)=lnrho_left(j)
          call gaunoise(ampllnrho(j),f,ilnrho,ilnrho)
        case ('gaussian-noise-x')
!
!  Noise, but just x-dependent.
!
          call gaunoise(ampllnrho(j),f,ilnrho,ilnrho)
          f(:,:,:,ilnrho)=spread(spread(f(:,4,4,ilnrho),2,my),3,mz) !(watch 1-D)
        case ('rho-jump-z', '2')
!
!  Density jump (for shocks).
!
          if (lroot) print*, 'init_lnrho: density jump; rho_left,right=', &
              rho_left(j), rho_right(j)
          if (lroot) print*, 'init_lnrho: density jump; widthlnrho=', &
              widthlnrho(j)
          do n=n1,n2; do m=m1,m2
            prof=0.5*(1.0+tanh(z(n)/widthlnrho(j)))
            f(l1:l2,m,n,ilnrho)=log(rho_left(j))+log(rho_left(j)/rho_right(j))*prof
          enddo; enddo
!
!  A*tanh(y/d) profile
!
        case ('tanhy')
          if (lroot) print*,'init_lnrho: tangential discontinuity'
          do m=m1,m2
            prof=ampllnrho(j)*tanh(y(m)/widthlnrho(j))
            do n=n1,n2
              f(l1:l2,m,n,ilnrho)=prof
            enddo
          enddo
        case ('sound-wave', '11')
!
!  Sound wave (should be consistent with hydro module).
!
          if (lroot) print*,'init_lnrho: x-wave in lnrho; ampllnrho=', &
              ampllnrho(j)
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnrho)=lnrho_const+ampllnrho(j)*sin(kx_lnrho(j)*x(l1:l2))
          enddo; enddo
        case ('sound-wave-exp')
!
!  Sound wave (should be consistent with hydro module).
!
          if (lroot) print*,'init_lnrho: x-wave in rho; ampllnrho=', &
              ampllnrho(j)
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnrho)=log(rho_const+amplrho(j)*sin(kx_lnrho(j)*x(l1:l2)))
          enddo; enddo
        case ('sound-wave2')
!
!  Sound wave (should be consistent with hydro module).
!
          if (lroot) print*,'init_lnrho: x-wave in lnrho; ampllnrho=', &
              ampllnrho(j)
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnrho)=lnrho_const+ampllnrho(j)*cos(kx_lnrho(j)*x(l1:l2))
          enddo; enddo
        case ('shock-tube', '13')
!
!  Shock tube test (should be consistent with hydro module).
!
          call information('init_lnrho','polytropic standing shock')
          do n=n1,n2; do m=m1,m2
            prof=0.5*(1.+tanh(x(l1:l2)/widthlnrho(j)))
            f(l1:l2,m,n,ilnrho)=log(rho_left(j))+ &
                (log(rho_right(j))-log(rho_left(j)))*prof
          enddo; enddo
        case ('sin-xy')
!
!  sin profile in x and y.
!
          call information('init_lnrho','lnrho=sin(x)*sin(y)')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnrho)=log(rho0) + &
                ampllnrho(j)*sin(kx_lnrho(j)*x(l1:l2))*sin(ky_lnrho(j)*y(m))
          enddo; enddo
        case ('sin-xy-rho')
!
!  sin profile in x and y, but in rho, not ln(rho).
!
          call information('init_lnrho','rho=sin(x)*sin(y)')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnrho)=log(rho0*(1+ &
                ampllnrho(j)*sin(kx_lnrho(j)*x(l1:l2))*sin(ky_lnrho(j)*y(m))))
          enddo; enddo
        case ('linear')
!
!  Linear profile in kk.xxx.
!
          call information('init_lnrho','linear profile')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnrho) = log(rho0) + &
                ampllnrho(j)*(kx_lnrho(j)*x(l1:l2)+ &
                ky_lnrho(j)*y(m)+kz_lnrho(j)*z(n))/ &
                sqrt(kx_lnrho(j)**2+ky_lnrho(j)**2+kz_lnrho(j)**2)
          enddo; enddo
        case ('planet')
!
!  Planet solution of Goodman, Narayan & Goldreich (1987).
!  (Simple 3-D)
!
          call planet(rbound,f,eps_planet,radius_lnrho(j), &
              gamma,cs20,rho0,widthlnrho(j),hh0)
        case ('planet_hc')
!
!  Planet solution of Goodman, Narayan & Goldreich (1987).
!  (3-D with hot corona)
!
          call planet_hc(amplrho(j),f,eps_planet, &
              radius_lnrho(j), gamma,cs20,rho0,widthlnrho(j))
        case ('Ferriere')
          call information('init_lnrho','Ferriere set in entropy')
        case ('Galactic-hs')
          call information('init_lnrho', &
              'Galactic hydrostatic equilibrium setup done in entropy')
        case ('compressive-shwave')
!  Should be consistent with density
          f(:,:,:,ilnrho) = log(rho_const + f(:,:,:,ilnrho))
!
        case default
!
!  Catch unknown values
!
          write(unit=errormsg,fmt=*) 'No such value for initlnrho(' &
                            //trim(iinit_str)//'): ',trim(initlnrho(j))
          call fatal_error('init_lnrho',errormsg)

        endselect
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
      if (lroot) print*,'init_lnrho: cs2bot,cs2top=',cs2bot,cs2top
!
!  If unlogarithmic density considered, take exp of lnrho resulting from
!  initlnrho
!
      if (ldensity_nolog) f(:,:,:,irho)=exp(f(:,:,:,ilnrho))
!
!  sanity check
!
      if (notanumber(f(l1:l2,m1:m2,n1:n2,ilnrho))) then
        call error('init_lnrho', 'Imaginary density values')
      endif
!
    endsubroutine init_lnrho
!***********************************************************************
    subroutine pencil_criteria_density()
!
!  All pencils that the Density module depends on are specified here.
!
!  19-11-04/anders: coded
!
      if (ldensity_nolog) lpenc_requested(i_rho)=.true.
      if (lcontinuity_gas) then
        lpenc_requested(i_divu)=.true.
        if (ldensity_nolog) then
          lpenc_requested(i_ugrho)=.true.
        else
          lpenc_requested(i_uglnrho)=.true.
        endif
      endif
      if (ldiff_shock) then
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_gshock)=.true.
        if (ldensity_nolog .or. ldiffusion_nolog) then
          lpenc_requested(i_grho)=.true.
          lpenc_requested(i_del2rho)=.true.
          if (ldiffusion_nolog) lpenc_requested(i_rho1)=.true.
        else
          lpenc_requested(i_glnrho)=.true.
          lpenc_requested(i_glnrho2)=.true.
          lpenc_requested(i_del2lnrho)=.true.
        endif
      endif
      if (ldiff_normal) then
        if (ldensity_nolog .or. ldiffusion_nolog) then
          lpenc_requested(i_del2rho)=.true.
          if (ldiffusion_nolog) lpenc_requested(i_rho1)=.true.
        else
          lpenc_requested(i_glnrho2)=.true.
          lpenc_requested(i_del2lnrho)=.true.
        endif
      endif
!
      lpenc_diagnos2d(i_lnrho)=.true.
      lpenc_diagnos2d(i_rho)=.true.
!
      if (idiag_rhom/=0 .or. idiag_rhomz/=0 .or. idiag_rhomy/=0 .or. &
           idiag_rhomx/=0 .or. idiag_rho2m/=0 .or. idiag_rhomin/=0 .or. &
           idiag_rhomax/=0 .or. idiag_rhomxy/=0 .or. idiag_rhomxz/=0 .or. &
           idiag_totmass/=0 .or. idiag_mass/=0 .or. idiag_drho2m/=0 .or. &
           idiag_drhom/=0) &
           lpenc_diagnos(i_rho)=.true.
      if (idiag_lnrho2m/=0) lpenc_diagnos(i_lnrho)=.true.
      if (idiag_ugrhom/=0) lpenc_diagnos(i_ugrho)=.true.
      if (idiag_uglnrhom/=0) lpenc_diagnos(i_uglnrho)=.true.
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
      if (ldensity_nolog) then
        if (lpencil_in(i_rho1)) lpencil_in(i_rho)=.true.
      else
        if (lpencil_in(i_rho)) lpencil_in(i_rho1)=.true.
      endif
      if (lpencil_in(i_grho)) then
        if (.not.ldensity_nolog) lpencil_in(i_rho)=.true.
      endif
      if (lpencil_in(i_glnrho)) then
        if (ldensity_nolog) lpencil_in(i_rho1)=.true.
      endif
      if (lpencil_in(i_uglnrho)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_glnrho)=.true.
      endif
      if (lpencil_in(i_ugrho)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_grho)=.true.
      endif
      if (lpencil_in(i_glnrho2)) lpencil_in(i_glnrho)=.true.
      if (lpencil_in(i_sglnrho)) then
        lpencil_in(i_sij)=.true.
        lpencil_in(i_glnrho)=.true.
      endif
      if (lpencil_in(i_del2lnrho)) then
        if (ldensity_nolog) then
          lpencil_in(i_rho1)=.true.
          lpencil_in(i_del2rho)=.true.
          lpencil_in(i_glnrho2)=.true.
        endif
      endif
!
    endsubroutine pencil_interdep_density
!***********************************************************************
    subroutine calc_pencils_density(f,p)
!
!  Calculate Density pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  19-11-04/anders: coded
!
      use Mpicomm, only: stop_it
      use Sub, only: grad,dot,dot2,u_dot_grad,del2,multmv,g2ij, dot_mn
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: tmp3
      integer :: i
!
      intent(inout) :: f,p
! lnrho
      if (lpencil(i_lnrho)) then
        if (ldensity_nolog) then
          p%lnrho=log(f(l1:l2,m,n,irho))
        else
          p%lnrho=f(l1:l2,m,n,ilnrho)
        endif
      endif
! rho1 and rho
      if (ldensity_nolog) then
        if (lpencil(i_rho)) then
          p%rho=f(l1:l2,m,n,irho)
          if (lcheck_negative_density .and. any(p%rho <= 0.)) &
            call fatal_error_local('calc_pencils_density', 'negative density detected')
        endif
        if (lpencil(i_rho1)) p%rho1=1.0/p%rho
      else
        if (lpencil(i_rho1)) p%rho1=exp(-f(l1:l2,m,n,ilnrho))
        if (lpencil(i_rho)) p%rho=1.0/p%rho1
      endif
! glnrho and grho
      if (lpencil(i_glnrho).or.lpencil(i_grho)) then
        call grad(f,ilnrho,tmp3)
        if (ldensity_nolog) then
          if (lpencil(i_grho)) p%grho=tmp3
          if (lpencil(i_glnrho)) then
            do i=1,3
              p%glnrho(:,i)=tmp3(:,i)*p%rho1
            enddo
          endif
        else
          if (lpencil(i_glnrho)) p%glnrho=tmp3
          if (lpencil(i_grho)) then
            if (irho/=0) then
              call grad(f,irho,p%grho)
            else
              do i=1,3
                p%grho(:,i)=p%rho*tmp3(:,i)
              enddo
            endif
          endif
        endif
      endif
! uglnrho
      if (lpencil(i_uglnrho)) then
        if (ldensity_nolog) then
          call dot(p%uu,p%glnrho,p%uglnrho)
        else
          call u_dot_grad(ilnrho,p%glnrho,p%uu,p%uglnrho)
        endif
      endif
! ugrho
      if (lpencil(i_ugrho)) then
        if (ldensity_nolog) then
          call u_dot_grad(ilnrho,p%grho,p%uu,p%ugrho)
        else
          call dot(p%uu,p%grho,p%ugrho)
        endif
      endif
! glnrho2
      if (lpencil(i_glnrho2)) call dot2(p%glnrho,p%glnrho2)
! del2rho
      if (lpencil(i_del2rho)) then
        if (ldensity_nolog) then
          call del2(f,ilnrho,p%del2rho)
        else
          if (irho/=0) then
            call del2(f,irho,p%del2rho)
          else
            if (headtt) &
                call fatal_error('calc_pencils_density',&
                'del2rho not available for logarithmic mass density')
          endif
        endif
      endif
! del2lnrho
      if (lpencil(i_del2lnrho)) then
        if (ldensity_nolog) then
          p%del2lnrho=p%rho1*p%del2rho-p%glnrho2
        else
          call del2(f,ilnrho,p%del2lnrho)
        endif
      endif
! hlnrho
      if (lpencil(i_hlnrho)) then
        if (ldensity_nolog) then
          if (headtt) then
            call fatal_error('calc_pencils_density', &
                             'hlnrho not available for non-logarithmic mass density')
          endif
        else
          call g2ij(f,ilnrho,p%hlnrho)
        endif
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
      if ( (.not.ldensity_nolog) .and. (irho/=0) ) &
          f(l1:l2,m1:m2,n1:n2,irho)=exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
!
    endsubroutine density_before_boundary
!***********************************************************************
    subroutine dlnrho_dt(df,p)
!
!  Continuity equation.
!  Calculate dlnrho/dt = - u.gradlnrho - divu
!
!   7-jun-02/axel: incoporated from subroutine pde
!
      use Diagnostics
      use Mpicomm, only: stop_it
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: fdiff, gshockglnrho, gshockgrho
!
      intent(in)  :: p
      intent(out) :: df
!
!  Identify module and boundary conditions.
!
      call timing('dlnrho_dt','entered',mnloop=.true.)
      if (headtt.or.ldebug) print*,'dlnrho_dt: SOLVE dlnrho_dt'
      if (headtt) call identify_bcs('lnrho',ilnrho)
!
!  Continuity equation.
!
      if (lcontinuity_gas) then
        if (ldensity_nolog) then
          df(l1:l2,m,n,irho)   = df(l1:l2,m,n,irho)   - p%ugrho - p%rho*p%divu
        else
          df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) - p%uglnrho - p%divu
        endif
      endif
!
!  Mass diffusion
!
      fdiff=0.0
!
      if (ldiff_normal) then  ! Normal diffusion operator
        if (ldensity_nolog) then
          fdiff = fdiff + diffrho*p%del2rho
        else
          if (ldiffusion_nolog) then
            fdiff = fdiff + p%rho1*diffrho*p%del2rho
          else
            fdiff = fdiff + diffrho*(p%del2lnrho+p%glnrho2)
          endif
        endif
        if (lfirst.and.ldt) diffus_diffrho=diffus_diffrho+diffrho
        if (headtt) print*,'dlnrho_dt: diffrho=', diffrho
      endif
!
!  Shock diffusion
!
      if (ldiff_shock) then
        if (ldensity_nolog) then
          call dot_mn(p%gshock,p%grho,gshockgrho)
          df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + &
               diffrho_shock*p%shock*p%del2rho + diffrho_shock*gshockgrho
        else
          if (ldiffusion_nolog) then
            call dot_mn(p%gshock,p%grho,gshockgrho)
            df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + p%rho1*( &
                 diffrho_shock*p%shock*p%del2rho + diffrho_shock*gshockgrho )
          else
            call dot_mn(p%gshock,p%glnrho,gshockglnrho)
            df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + &
                 diffrho_shock*p%shock*(p%del2lnrho+p%glnrho2) + &
                 diffrho_shock*gshockglnrho
          endif
        endif
        if (lfirst.and.ldt) diffus_diffrho=diffus_diffrho+diffrho_shock*p%shock
        if (headtt) print*,'dlnrho_dt: diffrho_shock=', diffrho_shock
      endif
!
!  Add diffusion term to continuity equation
!
      if (ldensity_nolog) then
        df(l1:l2,m,n,irho)   = df(l1:l2,m,n,irho)   + fdiff
      else
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + fdiff
      endif
!
!  Multiply diffusion coefficient by Nyquist scale.
!
      if (lfirst.and.ldt) then
        diffus_diffrho =diffus_diffrho *dxyz_2
        if (headtt.or.ldebug) then
          print*,'dlnrho_dt: max(diffus_diffrho ) =', maxval(diffus_diffrho)
        endif
      endif
!
!  2d-averages
!  Note that this does not necessarily happen with ldiagnos=.true.
!
      call timing('dlnrho_dt','before l2davgfirst',mnloop=.true.)
      if (l2davgfirst) then
        if (idiag_rhomxz/=0)    call ysum_mn_name_xz(p%rho,idiag_rhomxz)
        if (idiag_rhomxy/=0)    call zsum_mn_name_xy(p%rho,idiag_rhomxy)
      endif
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1
!
      if (l1davgfirst) then
        if (idiag_rhomz/=0)    call xysum_mn_name_z(p%rho,idiag_rhomz)
        if (idiag_rhomx/=0)    call yzsum_mn_name_x(p%rho,idiag_rhomx)
        if (idiag_rhomy/=0)    call xzsum_mn_name_y(p%rho,idiag_rhomy)
      endif
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
        if (idiag_uglnrhom/=0) call sum_mn_name(p%uglnrho,idiag_uglnrhom)
        if (idiag_dtd/=0) &
            call max_mn_name(diffus_diffrho/cdtv,idiag_dtd,l_dt=.true.)
      endif
      call timing('dlnrho_dt','finished',mnloop=.true.)
!
    endsubroutine dlnrho_dt
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
!  Reads and registers print parameters relevant for continuity equation.
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Diagnostics
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
!  (This needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_rhom=0; idiag_rho2m=0; idiag_lnrho2m=0
        idiag_drho2m=0; idiag_drhom=0
        idiag_ugrhom=0; idiag_uglnrhom=0
        idiag_rhomin=0; idiag_rhomax=0; idiag_dtd=0
        idiag_lnrhomphi=0; idiag_rhomphi=0
        idiag_rhomz=0; idiag_rhomy=0; idiag_rhomx=0
        idiag_rhomxy=0; idiag_rhomr=0; idiag_totmass=0; idiag_mass=0
        idiag_rhomxz=0
      endif
!
!  iname runs through all possible names that may be listed in print.in.
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
        call parse_name(iname,cname(iname),cform(iname),'uglnrhom', &
            idiag_uglnrhom)
        call parse_name(iname,cname(iname),cform(iname),'dtd',idiag_dtd)
        call parse_name(iname,cname(iname),cform(iname),'totmass',idiag_totmass)
        call parse_name(iname,cname(iname),cform(iname),'mass',idiag_mass)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rhomz', &
            idiag_rhomz)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'rhomy', &
            idiag_rhomy)
      enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'rhomx', &
            idiag_rhomx)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do inamexz=1,nnamexz
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'rhomxz', &
            idiag_rhomxz)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'rhomxy', &
            idiag_rhomxy)
      enddo
!
!  Write column where which density variable is stored.
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
      endif
!
    endsubroutine rprint_density
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
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Density.
!
        case ('rho')
          if (ldensity_nolog) then
            slices%yz =f(ix_loc,m1:m2,n1:n2,irho)
            slices%xz =f(l1:l2,iy_loc,n1:n2,irho)
            slices%xy =f(l1:l2,m1:m2,iz_loc,irho)
            slices%xy2=f(l1:l2,m1:m2,iz2_loc,irho)
            slices%ready=.true.
          else
            slices%yz =exp(f(ix_loc,m1:m2,n1:n2,ilnrho))
            slices%xz =exp(f(l1:l2,iy_loc,n1:n2,ilnrho))
            slices%xy =exp(f(l1:l2,m1:m2,iz_loc,ilnrho))
            slices%xy2=exp(f(l1:l2,m1:m2,iz2_loc,ilnrho))
            slices%ready=.true.
          endif
!
!  Logarithmic density.
!
        case ('lnrho')
          if (ldensity_nolog) then
            slices%yz =alog(f(ix_loc,m1:m2,n1:n2,irho))
            slices%xz =alog(f(l1:l2,iy_loc,n1:n2,irho))
            slices%xy =alog(f(l1:l2,m1:m2,iz_loc,irho))
            slices%xy2=alog(f(l1:l2,m1:m2,iz2_loc,irho))
            slices%ready=.true.
          else
            slices%yz =f(ix_loc,m1:m2,n1:n2,ilnrho)
            slices%xz =f(l1:l2,iy_loc,n1:n2,ilnrho)
            slices%xy =f(l1:l2,m1:m2,iz_loc,ilnrho)
            slices%xy2=f(l1:l2,m1:m2,iz2_loc,ilnrho)
            slices%ready=.true.
          endif
!
      endselect
!
    endsubroutine get_slices_density
!***********************************************************************
endmodule Density
