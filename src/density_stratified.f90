! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ldensity = .true.
! CPARAM logical, parameter :: lanelastic = .false.
! CPARAM logical, parameter :: lboussinesq = .false.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED rhos; rhos1
! PENCILS PROVIDED grhos(3); glnrhos(3); ugrhos
!
! PENCILS PROVIDED rho; rho1
! PENCILS PROVIDED grho(3); glnrho(3)
! PENCILS PROVIDED ekin
!
! PENCILS PROVIDED lnrho
! PENCILS PROVIDED ugrho; uglnrho
! PENCILS PROVIDED transprho
! PENCILS PROVIDED del2rho; del2lnrho; del6lnrho
! PENCILS PROVIDED hlnrho(3,3)
! PENCILS PROVIDED uij5glnrho(3); sglnrho(3)
!
!***************************************************************
module Density
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'density.h'
  integer :: pushpars2c, pushdiags2c  ! should be procedure pointer (F2003)
!
  character(len=labellen), dimension(ninit) :: initrho = 'nothing'
  real, dimension(ninit) :: amplrho = 0.0
  logical :: lconserve_mass = .false.
  logical :: ldiff_shock = .false.
  logical :: ldiff_hyper3_mesh = .false.
  logical :: lmassdiff_fix = .false.
  real :: density_floor = 0.0
  real :: diffrho_shock = 0.0
  real :: diffrho_hyper3_mesh = 0.0
  real, dimension(3) :: beta_glnrho_global=0.0, beta_glnrho_scaled=0.0
!
  namelist /density_init_pars/ initrho, amplrho, beta_glnrho_global, lconserve_mass, lmassdiff_fix
!
  namelist /density_run_pars/ density_floor, diffrho_hyper3_mesh, diffrho_shock, lconserve_mass, lmassdiff_fix
!
!  Diagnostic Variables
!
  integer :: idiag_mass = 0     ! DIAG_DOC: $\int\rho\,d^3x$
  integer :: idiag_rhomin = 0   ! DIAG_DOC: $\min\left|\rho\right|$
  integer :: idiag_rhomax = 0   ! DIAG_DOC: $\max\left|\rho\right|$
  integer :: idiag_drhom = 0    ! DIAG_DOC: $\langle\Delta\rho/\rho_0\rangle$
  integer :: idiag_drho2m = 0   ! DIAG_DOC: $\langle\left(\Delta\rho/\rho_0\right)^2\rangle$
  integer :: idiag_drhorms = 0  ! DIAG_DOC: $\langle\Delta\rho/\rho_0\rangle_{rms}$
  integer :: idiag_drhomax = 0  ! DIAG_DOC: $\max\left|\Delta\rho/\rho_0\right|$
!
!  xy-averages
!
  integer :: idiag_drhomz = 0    ! XYAVG_DOC: $\langle\Delta\rho/\rho_0\rangle_{xy}$
  integer :: idiag_drho2mz = 0   ! XYAVG_DOC: $\langle\left(\Delta\rho/\rho_0\right)^2\rangle_{xy}$
!
!  xz-averages
!
  integer :: idiag_drhomy = 0    ! XZAVG_DOC: $\langle\Delta\rho/\rho_0\rangle_{xz}$
  integer :: idiag_drho2my = 0   ! XZAVG_DOC: $\langle\left(\Delta\rho/\rho_0\right)^2\rangle_{xz}$
!
!  yz-averages
!
  integer :: idiag_drhomx = 0    ! YZAVG_DOC: $\langle\Delta\rho/\rho_0\rangle_{yz}$
  integer :: idiag_drho2mx = 0   ! YZAVG_DOC: $\langle\left(\Delta\rho/\rho_0\right)^2\rangle_{yz}$
!
!  y-averages
!
  integer :: idiag_drhomxz = 0   ! YAVG_DOC: $\langle\Delta\rho/\rho_0\rangle_y$
  integer :: idiag_drho2mxz = 0  ! YAVG_DOC: $\langle\left(\Delta\rho/\rho_0\right)^2\rangle_y$
!
!  z-averages
!
  integer :: idiag_drhomxy = 0   ! ZAVG_DOC: $\langle\Delta\rho/\rho_0\rangle_z$
  integer :: idiag_drho2mxy = 0  ! ZAVG_DOC: $\langle\left(\Delta\rho/\rho_0\right)^2\rangle_z$
  integer :: idiag_sigma = 0     ! ZAVG_DOC; $\Sigma\equiv\int\rho\,\mathrm{d}z$
!
!  Module Variables
!
  real, dimension(mz) :: rho0z = 0.0
  real, dimension(mz) :: dlnrho0dz = 0.0
  real :: mass0 = 0.0
!
!  Dummy Variables
!
  real, dimension(nz) :: glnrhomz = 0.0
  logical :: lcalc_glnrhomean = .false.
  logical :: lupw_lnrho = .false.
!
  contains
!***********************************************************************
    subroutine register_density()
!
!  Register Density module.
!
!  28-feb-13/ccyang: coded.
!
      use FArrayManager, only: farray_register_pde
!
!  Register relative density as dynamical variable rho.
!
      call farray_register_pde('rho', irho)
!
!  This module conflicts with Gravity module.
!
      if (lgrav) call fatal_error('register_density', 'conflicting with Gravity module. ')
!
!  This module does not consider logarithmic density.
!
      ldensity_nolog = .true.
!
!  Identify version number.
!
      if (lroot) call svn_id("$Id$")
!
    endsubroutine register_density
!***********************************************************************
    subroutine initialize_density(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  07-sep-14/ccyang: coded.
!
      use EquationOfState, only: select_eos_variable, get_stratz
      use SharedVariables, only: put_shared_variable
      use DensityMethods,  only: initialize_density_methods
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      integer :: ierr
!
!  Tell the equation of state that we're here and what f variable we use.
!
      call select_eos_variable('rho', irho)
      call initialize_density_methods
!
!  Get density stratification.
!
      if (lstratz) call get_stratz(z, rho0z, dlnrho0dz)
!
!  Disable the force-free considerations.
!
      call put_shared_variable('lffree', .false., ierr)
      if (ierr /= 0) call fatal_error('initialize_density', 'cannot share variable lffree. ')
!
!  Check the switches.
!
      shock: if (diffrho_shock > 0.0) then
        ldiff_shock = .true.
        if (lroot) print *, 'initialize_density: shock mass diffusion; diffrho_shock = ', diffrho_shock
      endif shock
!
      hyper3_mesh: if (diffrho_hyper3_mesh > 0.0) then
        ldiff_hyper3_mesh = .true.
        if (lroot) print *, 'initialize_density: mesh hyper-diffusion; diffrho_hyper3_mesh = ', diffrho_hyper3_mesh
      endif hyper3_mesh
!
!  Get the total mass.
!
      if (lconserve_mass) mass0 = total_mass(f)
!
    endsubroutine initialize_density
!***********************************************************************
    subroutine init_lnrho(f)
!
!  Initializes the density field.
!
!  27-feb-13/ccyang: coded.
!
      use Initcond, only: gaunoise
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      integer :: j
!
      init: do j = 1, ninit
        if (initrho(j) == 'nothing') cycle init
!
        init_cond: select case (initrho(j))
!
        case ('zero') init_cond
          f(:,:,:,irho) = 0.0
!
        case ('const') init_cond
          f(:,:,:,irho) = amplrho(j)
!
        case ('noise', 'gaussian') init_cond
          call gaunoise(amplrho(j), f, irho, irho)
!
        case default init_cond
          call fatal_error('init_lnrho', 'unknown initial condition ' // trim(initrho(j)))
!
        endselect init_cond
!
      enddo init
!
    endsubroutine init_lnrho
!***********************************************************************
    subroutine pencil_criteria_density()
!
!  All pencils that the Density module depends on are specified here.
!
!  02-nov-13/ccyang: coded.
!
      lpenc_requested(i_rhos) = .true.
      lpenc_requested(i_uu) = .true.
      lpenc_requested(i_divu) = .true.
      lpenc_requested(i_ugrhos) = .true.
!
      shock: if (ldiff_shock) then
        lpenc_requested(i_shock) = .true.
        lpenc_requested(i_gshock) = .true.
        lpenc_requested(i_grhos) = .true.
      endif shock
!
      massdiff: if (lmassdiff_fix) then
        hydro: if (lhydro) then
          lpenc_requested(i_rhos1) = .true.
          lpenc_requested(i_uu) = .true.
        endif hydro
        if (lthermal_energy) lpenc_requested(i_u2) = .true.
      endif massdiff
!
!  Diagnostic Pencils
!
      if (idiag_mass /= 0 .or. idiag_rhomin /= 0 .or. idiag_rhomax /= 0) lpenc_diagnos(i_rho) = .true.
      if (idiag_sigma /= 0) lpenc_diagnos(i_rho) = .true.
!
    endsubroutine pencil_criteria_density
!***********************************************************************
    subroutine pencil_interdep_density(lpencil_in)
!
!  Interdependency among pencils from the Density module is specified here.
!
!  02-nov-13/ccyang: coded.
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_rhos1)) lpencil_in(i_rhos) = .true.
!
      glnrhos: if (lpencil_in(i_glnrhos)) then
        lpencil_in(i_rhos1) = .true.
        lpencil_in(i_grhos) = .true.
      endif glnrhos
!
      ugrhos: if (lpencil_in(i_ugrhos)) then
        lpencil_in(i_grhos) = .true.
        lpencil_in(i_uu) = .true.
      endif ugrhos
!
      if (lpencil_in(i_rho)) lpencil_in(i_rhos) = .true.
!
      if (lpencil_in(i_rho1)) lpencil_in(i_rho) = .true.
!
      grho: if (lpencil_in(i_grho)) then
        lpencil_in(i_rho) = .true.
        lpencil_in(i_grhos) = .true.
      endif grho
!
      glnrho: if (lpencil_in(i_glnrho)) then
        lpencil_in(i_rho1) = .true.
        lpencil_in(i_grho) = .true.
      endif glnrho
!
      ekin: if (lpencil_in(i_ekin)) then
        lpencil_in(i_rho)=.true.
        lpencil_in(i_u2)=.true.
      endif ekin
!
    endsubroutine pencil_interdep_density
!***********************************************************************
    subroutine calc_pencils_density(f, p)
!
!  Calculate Density pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  02-nov-13/ccyang: coded.
!
      use Sub, only: grad, u_dot_grad
      use EquationOfState, only: cs20
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      type(pencil_case), intent(inout) :: p
!
! rhos: density scaled by stratification
!
      if (lpencil(i_rhos)) p%rhos = 1.0 + f(l1:l2,m,n,irho)
!
! rhos1
!
      if (lpencil(i_rhos1)) p%rhos1 = 1.0 / p%rhos
!
! grhos
!
      if (lpencil(i_grhos)) call grad(f, irho, p%grhos)
!
! glnrhos
!
      if (lpencil(i_glnrhos)) p%glnrhos = spread(p%rhos1,2,3) * p%grhos
!
! ugrhos
!
      if (lpencil(i_ugrhos)) call u_dot_grad(f, irho, p%grhos, p%uu, p%ugrhos)
!
! rho
!
      if (lpencil(i_rho)) p%rho = rho0z(n) * p%rhos
!
! rho1
!
      if (lpencil(i_rho1)) p%rho1 = 1.0 / p%rho
!
! grho
!
      grho: if (lpencil(i_grho)) then
        p%grho = rho0z(n) * p%grhos
        p%grho(:,3) = p%grho(:,3) + dlnrho0dz(n) * p%rho
      endif grho
!
! glnrho
!
      if (lpencil(i_glnrho)) p%glnrho = spread(p%rho1,2,3) * p%grho
!
! ekin
!
      if (lpencil(i_ekin)) p%ekin = 0.5 * p%rho * p%u2
!
!  Currently not required pencils.
!
      lnrho: if (lpencil(i_lnrho)) then
        call fatal_error('calc_pencils_density', 'lnrho is not available. ')
        p%lnrho = 0.0
      endif lnrho
!
      ugrho: if (lpencil(i_ugrho)) then
        call fatal_error('calc_pencils_density', 'ugrho is not available. ')
        p%ugrho = 0.0
      endif ugrho
!
      uglnrho: if (lpencil(i_uglnrho)) then
        call fatal_error('calc_pencils_density', 'uglnrho is not available. ')
        p%uglnrho = 0.0
      endif uglnrho
!
      transprho: if (lpencil(i_transprho)) then
        call fatal_error('calc_pencils_density', 'transprho is not available. ')
        p%transprho = 0.0
      endif transprho
!
      del2rho: if (lpencil(i_del2rho)) then
        call fatal_error('calc_pencils_density', 'del2rho is not available. ')
        p%del2rho = 0.0
      endif del2rho
!
      del2lnrho: if (lpencil(i_del2lnrho)) then
        call fatal_error('calc_pencils_density', 'del2lnrho is not available. ')
        p%del2lnrho = 0.0
      endif del2lnrho
!
      del6lnrho: if (lpencil(i_del6lnrho)) then
        call fatal_error('calc_pencils_density', 'del6lnrho is not available. ')
        p%del6lnrho = 0.0
      endif del6lnrho
!
      hlnrho: if (lpencil(i_hlnrho)) then
        call fatal_error('calc_pencils_density', 'hlnrho is not available. ')
        p%hlnrho = 0.0
      endif hlnrho
!
      uij5glnrho: if (lpencil(i_uij5glnrho)) then
        call fatal_error('calc_pencils_density', 'uij5glnrho is not available. ')
        p%uij5glnrho = 0.0
      endif uij5glnrho
!
      sglnrho: if (lpencil(i_sglnrho)) then
        call fatal_error('calc_pencils_density', 'sglnrho is not available. ')
        p%sglnrho = 0.0
      endif sglnrho
!
    endsubroutine calc_pencils_density
!***********************************************************************
    subroutine dlnrho_dt(f,df,p)
!
!  Continuity equation.
!
!  02-nov-13/ccyang: coded.
!  24-aug-14/ccyang: add mass diffusion correction.
!
      use Sub, only: identify_bcs, dot_mn, del2
      use Deriv, only: der6
      use Special, only: special_calc_density
      use Diagnostics
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type(pencil_case), intent(in) :: p
!
      real, dimension(nx) :: fdiff, del2rhos, penc, src_density
      integer :: j
!
!  Start the clock for this procedure.
!
      call timing('dlnrho_dt', 'entered', mnloop=.true.)
!
!  Identify module and boundary conditions.
!
      if (headtt .or. ldebug) print *, 'dlnrho_dt: SOLVE the continuity equation. '
      if (headtt) call identify_bcs('rho', irho)
!
!  Find the rate of change.
!
      penc = p%uu(:,3) * dlnrho0dz(n)
      df(l1:l2,m,n,irho) = df(l1:l2,m,n,irho) - p%ugrhos - p%rhos * (p%divu + penc)
      if (lfirst .and. ldt) then
        src_density = abs(penc)
        maxsrc = max(maxsrc, src_density)
      endif
!
      fdiff = 0.0
!
!  Shock mass diffusion
!
      shock: if (ldiff_shock) then
        call del2(f, irho, del2rhos)
        call dot_mn(p%gshock, p%grhos, penc)
        fdiff = fdiff + diffrho_shock * (p%shock * del2rhos + penc)
        if (lfirst .and. ldt) diffus_diffrho = diffus_diffrho + diffrho_shock * dxyz_2 * p%shock
      endif shock
!
!  Mesh hyper-diffusion
!
      hyper3: if (ldiff_hyper3_mesh) then
        dir: do j = 1, 3
          call der6(f, irho, penc, j, ignoredx=.true.)
          fdiff = fdiff + diffrho_hyper3_mesh * penc * dline_1(:,j)
        enddo dir
        if (lfirst .and. ldt) diffus_diffrho3 = diffus_diffrho3 &
                                              + diffrho_hyper3_mesh * (abs(dline_1(:,1)) + abs(dline_1(:,2)) + abs(dline_1(:,3)))
      endif hyper3
!
      df(l1:l2,m,n,irho) = df(l1:l2,m,n,irho) + fdiff
!
!  Mass diffusion corrections.
!
      massdiff: if (lmassdiff_fix) then
        hydro: if (lhydro) then
          penc = fdiff * p%rhos1
          forall(j = iux:iuz) df(l1:l2,m,n,j) = df(l1:l2,m,n,j) - penc * p%uu(:,j-iuu+1)
        endif hydro
        energy: if (lthermal_energy) then
          df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) - 0.5 * rho0z(n) * fdiff * p%u2
        elseif (lenergy) then energy
          call fatal_error('dlnrho_dt', 'mass diffusion correction for this energy equation is not implemented. ')
        endif energy
      endif massdiff
!
!  Call optional user-defined calculations.
!
      if (lspecial) call special_calc_density(f, df, p)
!
!  2D averages
!
      avg_2d: if (l2davgfirst) then
        penc = f(l1:l2,m,n,irho)
        if (idiag_drhomxz /= 0) call ysum_mn_name_xz(penc, idiag_drhomxz)
        if (idiag_drho2mxz /= 0) call ysum_mn_name_xz(penc**2, idiag_drho2mxz)
        if (idiag_drhomxy /= 0) call zsum_mn_name_xy(penc, idiag_drhomxy)
        if (idiag_drho2mxy /= 0) call zsum_mn_name_xy(penc**2, idiag_drho2mxy)
        if (idiag_sigma /= 0) call zsum_mn_name_xy(p%rho, idiag_sigma, lint=.true.)
      endif avg_2d
!
!  1D averages
!
      avg_1d: if (l1davgfirst) then
        penc = f(l1:l2,m,n,irho)
!       xy-averages
        if (idiag_drhomz /= 0) call xysum_mn_name_z(penc, idiag_drhomz)
        if (idiag_drho2mz /= 0) call xysum_mn_name_z(penc**2, idiag_drho2mz)
!       xz-averages
        if (idiag_drhomy /= 0) call xzsum_mn_name_y(penc, idiag_drhomy)
        if (idiag_drho2my /= 0) call xzsum_mn_name_y(penc**2, idiag_drho2my)
!       yz-averages
        if (idiag_drhomx /= 0) call yzsum_mn_name_x(penc, idiag_drhomx)
        if (idiag_drho2mx /= 0) call yzsum_mn_name_x(penc**2, idiag_drho2mx)
      endif avg_1d
!
!  Diagnostics
!
      diagnos: if (ldiagnos) then
!
        rho: if (idiag_mass /= 0 .or. idiag_rhomin /= 0 .or. idiag_rhomax /= 0) then
          penc = p%rho
          if (idiag_mass /= 0) call integrate_mn_name(penc, idiag_mass)
          if (idiag_rhomin /= 0) call max_mn_name(-penc, idiag_rhomin, lneg=.true.)
          if (idiag_rhomax /= 0) call max_mn_name(penc, idiag_rhomax)
        endif rho
!
        drho: if (idiag_drhom /= 0 .or. idiag_drho2m /= 0 .or. idiag_drhorms /= 0 .or. idiag_drhomax /= 0) then
          penc = f(l1:l2,m,n,irho)
          if (idiag_drhom /= 0) call sum_mn_name(penc, idiag_drhom)
          if (idiag_drho2m /= 0) call sum_mn_name(penc**2, idiag_drho2m)
          if (idiag_drhorms /= 0) call sum_mn_name(penc**2, idiag_drhorms, lsqrt=.true.)
          if (idiag_drhomax /= 0) call max_mn_name(abs(penc), idiag_drhomax)
        endif drho
!
      endif diagnos
!
!  Stop the clock.
!
      call timing('dlnrho_dt', 'finished', mnloop=.true.)
!
    endsubroutine dlnrho_dt
!***********************************************************************
    subroutine calc_diagnostics_density(f,p)

      real, dimension (mx,my,mz,mfarray) :: f
      type(pencil_case) :: p

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)

    endsubroutine calc_diagnostics_density
!***********************************************************************
    subroutine split_update_density(f)
!
!  Integrate operator split terms and/or perform post-time-step
!  processing.
!
!  18-jun-14/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      real :: a, a1
!
!  Conserve mass and momentum.
!
      consrv: if (lconserve_mass) then
        a = mass0 / total_mass(f)
        a1 = 1.0 / a
        zscan: do n = n1, n2
          yscan: do m = m1, m2
            f(l1:l2,m,n,irho) = a * f(l1:l2,m,n,irho) + (a - 1.0)
            if (lhydro) f(l1:l2,m,n,iux:iuz) = a1 * f(l1:l2,m,n,iux:iuz)
          enddo yscan
        enddo zscan
      endif consrv
!
    endsubroutine split_update_density
!***********************************************************************
    subroutine impose_density_floor(f)
!
!  Check negative density or impose a density floor.
!
!  27-feb-13/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      integer :: i, j, k
      real :: drhos_min
!
!  Impose density floor or
!
      dens_floor: if (density_floor > 0.0) then
        scan1: do k = n1, n2
          drhos_min = (density_floor - rho0z(k)) / rho0z(k)
          where(f(l1:l2,m1:m2,k,irho) < drhos_min) f(l1:l2,m1:m2,k,irho) = drhos_min
        enddo scan1
!
!  Trap any negative density.
!
      else dens_floor
        neg_dens: if (any(f(l1:l2,m1:m2,n1:n2,irho) <= -1.0)) then
          scan2z: do k = n1, n2
            scan2y: do j = m1, m2
              scan2x: do i = l1, l2
                if (f(i,j,k,irho) <= -1.0) print 10, rho0z(k) * (1.0 + f(i,j,k,irho)), x(i), y(j), z(k)
                10 format (1x, 'Negative density ', es13.6, ' is detected at x = ', es10.3, ', y = ', es10.3, ', and z = ', es10.3)
              enddo scan2x
            enddo scan2y
          enddo scan2z
          call fatal_error_local('impose_density_floor', 'detected negative density. ')
        endif neg_dens
!
      endif dens_floor
!
    endsubroutine impose_density_floor
!***********************************************************************
    subroutine read_density_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=density_init_pars, IOSTAT=iostat)
!
    endsubroutine read_density_init_pars
!***********************************************************************
    subroutine write_density_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=density_init_pars)
!
    endsubroutine write_density_init_pars
!***********************************************************************
    subroutine read_density_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=density_run_pars, IOSTAT=iostat)
!
    endsubroutine read_density_run_pars
!***********************************************************************
    subroutine write_density_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=density_run_pars)
!
    endsubroutine write_density_run_pars
!***********************************************************************
    subroutine rprint_density(lreset, lwrite)
!
!  Reads and registers print parameters relevant for continuity equation.
!
!  27-feb-13/ccyang: coded.
!
      use Diagnostics, only: parse_name
      use FArrayManager, only: farray_index_append
!
      logical, intent(in) :: lreset
      logical, intent(in), optional :: lwrite
!
      logical :: lwr
      integer :: iname
!
      lwr = .false.
      if (present(lwrite)) lwr = lwrite
!
!  Reset everything in case of reset.
!
      reset: if (lreset) then
!       Diagnostics
        idiag_mass = 0
        idiag_rhomin = 0
        idiag_rhomax = 0
        idiag_drhom = 0
        idiag_drho2m = 0
        idiag_drhorms = 0
        idiag_drhomax = 0
!       xy-averages
        idiag_drhomz = 0
        idiag_drho2mz = 0
!       xz-averages
        idiag_drhomy = 0
        idiag_drho2my = 0
!       yz-averages
        idiag_drhomx = 0
        idiag_drho2mx = 0
!       y-averages
        idiag_drhomxz = 0
        idiag_drho2mxz = 0
!       z-averages
        idiag_sigma = 0
        idiag_drhomxy = 0
        idiag_drho2mxy = 0
      endif reset
!
!  Check for diagnostic variables listed in print.in.
!
      whole: do iname = 1, nname
        call parse_name(iname, cname(iname), cform(iname), 'mass', idiag_mass)
        call parse_name(iname, cname(iname), cform(iname), 'rhomin', idiag_rhomin)
        call parse_name(iname, cname(iname), cform(iname), 'rhomax', idiag_rhomax)
        call parse_name(iname, cname(iname), cform(iname), 'drhom', idiag_drhom)
        call parse_name(iname, cname(iname), cform(iname), 'drho2m', idiag_drho2m)
        call parse_name(iname, cname(iname), cform(iname), 'drhorms', idiag_drhorms)
        call parse_name(iname, cname(iname), cform(iname), 'drhomax', idiag_drhomax)
      enddo whole
!
!  Check for xy-averages listed in xyaver.in.
!
      xyaver: do iname = 1, nnamez
        call parse_name(iname, cnamez(iname), cformz(iname), 'drhomz', idiag_drhomz)
        call parse_name(iname, cnamez(iname), cformz(iname), 'drho2mz', idiag_drho2mz)
      enddo xyaver
!
!  Check for xz-averages listed in xzaver.in.
!
      xzaver: do iname = 1, nnamey
        call parse_name(iname, cnamey(iname), cformy(iname), 'drhomy', idiag_drhomy)
        call parse_name(iname, cnamey(iname), cformy(iname), 'drho2my', idiag_drho2my)
      enddo xzaver
!
!  Check for yz-averages listed in yzaver.in.
!
      yzaver: do iname = 1, nnamex
        call parse_name(iname, cnamex(iname), cformx(iname), 'drhomx', idiag_drhomx)
        call parse_name(iname, cnamex(iname), cformx(iname), 'drho2mx', idiag_drho2mx)
      enddo yzaver
!
!  Check for y-averages listed in yaver.in.
!
      yaver: do iname = 1, nnamexz
        call parse_name(iname, cnamexz(iname), cformxz(iname), 'drhomxz', idiag_drhomxz)
        call parse_name(iname, cnamexz(iname), cformxz(iname), 'drho2mxz', idiag_drho2mxz)
      enddo yaver
!
!  Check for z-averages listed in zaver.in.
!
      zaver: do iname = 1, nnamexy
        call parse_name(iname, cnamexy(iname), cformxy(iname), 'drhomxy', idiag_drhomxy)
        call parse_name(iname, cnamexy(iname), cformxy(iname), 'drho2mxy', idiag_drho2mxy)
        call parse_name(iname, cnamexy(iname), cformxy(iname), 'sigma', idiag_sigma)
      enddo zaver
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then 
        where(cnamev=='rho'.or.cnamev=='drho') cformv='DEFINED'
      endif
!
!  Write column where which density variable is stored.
!
      if (lwr) then
        call farray_index_append('ilnrho',0)
        call farray_index_append('irho',irho)
      endif
!
    endsubroutine rprint_density
!***********************************************************************
    subroutine get_slices_density(f,slices)
!
!  Write slices for animation of Density variables.
!
!  27-feb-13/ccyang: coded.
!
      use Slices_methods, only: assign_slices_scal

      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      type(slice_data), intent(inout) :: slices
!
      var: select case (trim(slices%name))
!
      case ('drho');  call assign_slices_scal(slices,f,irho)
!
      case ('rho') var
        if (lwrite_slice_yz)  slices%yz  = (1.+f(ix_loc,m1:m2, n1:n2,  irho))*spread(rho0z(n1:n2),1,ny)
        if (lwrite_slice_xz)  slices%xz  = (1.+f(l1:l2, iy_loc,n1:n2,  irho))*spread(rho0z(n1:n2),1,nx)
        if (lwrite_slice_xz2) slices%xz2 = (1.+f(l1:l2, iy2_loc,n1:n2,  irho))*spread(rho0z(n1:n2),1,nx)
        if (lwrite_slice_xy)  slices%xy  = (1.+f(l1:l2, m1:m2, iz_loc, irho))*rho0z(iz_loc)
        if (lwrite_slice_xy2) slices%xy2 = (1.+f(l1:l2, m1:m2, iz2_loc,irho))*rho0z(iz2_loc)
        if (lwrite_slice_xy3) slices%xy3 = (1.+f(l1:l2, m1:m2, iz3_loc,irho))*rho0z(iz3_loc)
        if (lwrite_slice_xy4) slices%xy4 = (1.+f(l1:l2, m1:m2, iz4_loc,irho))*rho0z(iz4_loc)
        slices%ready = .true.
!
      endselect var
!
    endsubroutine get_slices_density
!***********************************************************************
    subroutine dynamical_diffusion(uc)
!
!  Dynamically set mass diffusion coefficient given fixed mesh Reynolds number.
!
!  28-feb-13/ccyang: coded
!
!  Input Argument
!      uc
!          Characteristic velocity of the system.
!
      real, intent(in) :: uc
!
!  Hyper-diffusion coefficient
!
      if (ldiff_hyper3_mesh) diffrho_hyper3_mesh = pi5_1 * uc / re_mesh / sqrt(3.0)
!
    endsubroutine dynamical_diffusion
!***********************************************************************
!***********************************************************************
!
!  LOCAL ROUTINES GO BELOW HERE.
!
!***********************************************************************
!***********************************************************************
    real function total_mass(f)
!
!  Gets the total mass.
!
!  18-jun-14/ccyang: coded
!
      use Mpicomm, only: mpiallreduce_sum
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      real :: mass_loc, a
!
      mass_loc = 0.0
      zscan: do n = n1, n2
        a = rho0z(n) * zprim(n)
        do m = m1, m2
          mass_loc = mass_loc + sum(a * yprim(m) * xprim(l1:l2) * (1.0 + f(l1:l2,m,n,irho)))
        enddo
      enddo zscan
      call mpiallreduce_sum(mass_loc, total_mass)
!
    endfunction total_mass
!***********************************************************************
!***********************************************************************
!
!  DUMMY BUT PUBLIC ROUTINES GO BELOW HERE.
!
!***********************************************************************
!***********************************************************************
    subroutine anelastic_after_mn(f, p, df, mass_per_proc)
!
!  Dummy
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(in) :: df
      real, dimension(1), intent(in) :: mass_per_proc
      type(pencil_case), intent(in) :: p
!
      call keep_compiler_quiet(f, df)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(mass_per_proc)
!
    endsubroutine anelastic_after_mn
!***********************************************************************
    subroutine boussinesq(f)
!
!  dummy routine for the Boussinesq approximation
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine boussinesq
!***********************************************************************
    subroutine density_after_boundary(f)
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
  endsubroutine density_after_boundary
!***********************************************************************
    subroutine density_before_boundary(f)
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine density_before_boundary
!***********************************************************************
    subroutine get_init_average_density(f, init_average_density)
!
!  10-dec-09/piyali: added to pass initial average density
!
    real, dimension(mx,my,mz,mfarray), intent(in) :: f
    real, intent(in) :: init_average_density
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(init_average_density)
!
    endsubroutine get_init_average_density
!***********************************************************************
    subroutine get_slices_pressure(f, slices)
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      type(slice_data), intent(in) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_pressure
!***********************************************************************
    real function mean_density(f)
!
!  Calculate mean density of the whole box.
!
!  06-mar-15/ccyang: adapted.
!
      use Mpicomm, only: mpiallreduce_sum
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      integer :: m, n
      real :: tmp
!
      mean_density = 0.0
!
      zscan: do n = n1, n2
        tmp = 0.0
        yscan: do m = m1, m2
          tmp = tmp + sum((1.0 + f(l1:l2,m,n,irho)) * dVol_x(l1:l2)) * dVol_y(m)
        enddo yscan
        mean_density = mean_density + tmp * rho0z(n) * dVol_z(n)
      enddo zscan
!
      mean_density = mean_density / box_volume
!
      comm: if (ncpus > 1) then
        call mpiallreduce_sum(mean_density, tmp)
        mean_density = tmp
      endif comm
!
    endfunction mean_density
!***********************************************************************
    subroutine update_char_vel_density(f)
!
!  Updates characteristic velocity for slope-limited diffusion.
!
!  21-oct-15/MR: coded
!
      real, dimension(mx,my,mz,mfarray), intent(INOUT) :: f
!
      if (lslope_limit_diff) call fatal_error('update_char_vel_density', 'not implemented')
!
    endsubroutine update_char_vel_density
!***********************************************************************
    subroutine impose_density_ceiling(f)
!
!  Dummy routine.
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f

      call keep_compiler_quiet(f)

    endsubroutine impose_density_ceiling
!***********************************************************************
endmodule Density
