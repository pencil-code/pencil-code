! $Id$
!
!  This module takes care of everything related to dust velocity
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ldustvelocity = .true.
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED divud(ndustspec); ood(3,ndustspec); od2(ndustspec)
! PENCILS PROVIDED oud(ndustspec); ud2(ndustspec); udij(3,3,ndustspec)
! PENCILS PROVIDED sdij(3,3,ndustspec); udgud(3,ndustspec); uud(3,ndustspec)
! PENCILS PROVIDED del2ud(3,ndustspec); del6ud(3,ndustspec)
! PENCILS PROVIDED graddivud(3,ndustspec)
!
!***************************************************************

module Dustvelocity
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'dustvelocity.h'
!ajwm SHOULDN'T REALLY BE SHARED
!ajwm but are used consistently with the Dustdensity module
!ajwm - not good but for reasons of dust density / velocity interaction
  public :: rhods, surfd, mdplus, mdminus
  public :: ad, scolld, ustcst, tausd1, tausd
  public :: unit_md, dust_chemistry, mumon, mmon, mi, md
!
  integer, parameter :: nvisc_max=4
  complex, dimension (7) :: coeff=0.
  real, dimension(ndustspec,ndustspec) :: scolld
  real, dimension(nx,ndustspec) :: tausd1
  real, dimension(ndustspec) :: md=1.0, mdplus, mdminus, ad=0.
  real, dimension(ndustspec) :: surfd, mi, rhodsad1
  real, dimension(ndustspec) :: tausd=1.0, betad=0.0
  real :: betad0=0.
  real, dimension(ndustspec) :: nud=0.0, nud_hyper3=0.0, nud_shock=0.0, nud_hyper3_mesh=5.0
  real :: uudx0=0.0, uudy0=0.0, uudz0=0.0
  real :: ampluud=0.0, ampl_udx=0.0, ampl_udy=0.0, ampl_udz=0.0
  real :: phase_udx=0.0, phase_udy=0.0, phase_udz=0.0
  real :: kx_uud=1.0, ky_uud=1.0, kz_uud=1.0
  real :: rhods=1.0, nd0=1.0, md0=1.0, rhod0=1.0, mu_ext=0.
  real :: ad0=0.0, ad1=0.0, dimd1=0.333333, deltamd=1.0
  real :: nud_all=0.0, betad_all=0.0, tausd_all=0.0
  real :: viscd_exponent=0.0, adref_nud=0.0
  real :: mmon, mumon, mumon1, surfmon, ustcst, unit_md=1.0
  real :: beta_dPdr_dust=0.0, beta_dPdr_dust_scaled=0.0,cdtd=0.2
  real :: gravx_dust=0.0
  real :: Omega_pseudo=0.0, u0_gas_pseudo=0.0, tausgmin=0.0, tausg1max=0.0
  real :: dust_pressure_factor=1.
  real :: shorttauslimit=0.0, shorttaus1limit=0.0
  real :: scaleHtaus=1.0, z0taus=0.0, widthtaus=1.0
  logical :: llin_radiusbins=.false., llog_massbins=.true.
  logical :: ladvection_dust=.true.,lcoriolisforce_dust=.true.
  logical :: ldragforce_dust=.true.,ldragforce_gas=.false.
  logical :: ldust_pressure=.false., reinitialize_uud=.false.
  logical :: ldustvelocity_shorttausd=.false., lvshear_dust_global_eps=.false.
  logical :: ldustcoagulation=.false., ldustcondensation=.false.
  logical :: lviscd_simplified=.false., lviscd_nud_const=.false.
  logical :: lviscd_shock=.false., lviscd_shock_simplified=.false.
  logical :: lviscd_hyper3_simplified=.false.
  logical :: lviscd_hyper3_rhod_nud_const=.false.
  logical :: lviscd_hyper3_nud_const=.false.
  logical :: lviscd_hyper3_polar=.false.
  logical :: lviscd_hyper3_mesh=.false.
  logical :: lstokes_highspeed_corr=.false.
  logical :: lpifactor1=.false., lpifactor2=.false.
  character (len=labellen), dimension(ninit) :: inituud='nothing'
  character (len=labellen) :: borderuud='nothing'
  character (len=labellen), dimension(nvisc_max) :: iviscd=''
  character (len=labellen) :: draglaw='epstein_cst', viscd_law='const'
  character (len=labellen) :: dust_geometry='sphere', dust_chemistry='nothing'
  character (len=labellen) :: iefficiency_type='nothing'
!
  namelist /dustvelocity_init_pars/ &
      uudx0, uudy0, uudz0, ampl_udx, ampl_udy, ampl_udz, &
      phase_udx, phase_udy, phase_udz, rhods, mu_ext, &
      md0, ad0, ad1, deltamd, draglaw, viscd_law, viscd_exponent, &
      ampluud, inituud, &
      kx_uud, ky_uud, kz_uud, Omega_pseudo, u0_gas_pseudo, &
      dust_chemistry, dust_geometry, tausd, gravx_dust, &
      beta_dPdr_dust, coeff,  ldustcoagulation, ldustcondensation, &
      llin_radiusbins, llog_massbins, &
      lvshear_dust_global_eps, cdtd, &
      ldustvelocity_shorttausd, scaleHtaus, z0taus, betad0,&
      lstokes_highspeed_corr, iefficiency_type
!
  namelist /dustvelocity_run_pars/ &
      nud, nud_all, iviscd, betad, betad_all, tausd, tausd_all, draglaw, &
      viscd_law, viscd_exponent, adref_nud, inituud, &
      ldragforce_dust, ldragforce_gas, ldustvelocity_shorttausd, mu_ext, &
      ladvection_dust, lcoriolisforce_dust, gravx_dust, &
      reinitialize_uud, beta_dPdr_dust, tausgmin, cdtd, nud_shock, &
      nud_hyper3, nud_hyper3_mesh, scaleHtaus, z0taus, widthtaus, shorttauslimit,&
      lstokes_highspeed_corr, lpifactor1, lpifactor2
!
  integer :: idiag_ekintot_dust=0
  integer, dimension(ndustspec) :: idiag_ud2m=0
  integer, dimension(ndustspec) :: idiag_udxm=0, idiag_udym=0, idiag_udzm=0
  integer, dimension(ndustspec) :: idiag_udx2m=0, idiag_udy2m=0, idiag_udz2m=0
  integer, dimension(ndustspec) :: idiag_udm2=0, idiag_oudm=0, idiag_od2m=0
  integer, dimension(ndustspec) :: idiag_udrms=0, idiag_udmax=0, idiag_odrms=0
  integer, dimension(ndustspec) :: idiag_odmax=0, idiag_rdudmax=0
  integer, dimension(ndustspec) :: idiag_udxmz=0, idiag_udymz=0, idiag_udzmz=0
  integer, dimension(ndustspec) :: idiag_udx2mz=0, idiag_udy2mz=0
  integer, dimension(ndustspec) :: idiag_udz2mz=0
  integer, dimension(ndustspec) :: idiag_udmx=0, idiag_udmy=0, idiag_udmz=0
  integer, dimension(ndustspec) :: idiag_udxmxy=0, idiag_udymxy=0
  integer, dimension(ndustspec) :: idiag_udzmxy=0
  integer, dimension(ndustspec) :: idiag_divud2m=0, idiag_epsKd=0
  integer, dimension(ndustspec) :: idiag_dtud=0, idiag_dtnud=0
  integer, dimension(ndustspec) :: idiag_rdudxm=0, idiag_rdudym=0
  integer, dimension(ndustspec) :: idiag_rdudzm=0, idiag_rdudx2m=0
!
! Auxiliary variables:
!
  real, dimension (nx) :: diffus_nud,diffus_nud3,advec_uud,advec_hypermesh_uud

  contains
!***********************************************************************
    subroutine register_dustvelocity()
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
!  18-mar-03/axel+anders: adapted from hydro
!
      use FArrayManager
      use General, only: itoa
      use SharedVariables, only: put_shared_variable
!
      integer :: k, uud_tmp
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
          "$Id$")
!
      call farray_index_append('nuud',ndustspec)
      call farray_register_pde('uud',uud_tmp,vector=3,array=ndustspec)
      do k=1, ndustspec
        iuud(k) = uud_tmp + (k-1)*3
        iudx(k) = iuud(k)
        iudy(k) = iuud(k)+1
        iudz(k) = iuud(k)+2
      enddo
!
!  Need deltamd for normalization purposes in dustdensity.
!
      call put_shared_variable('llin_radiusbins',llin_radiusbins,caller='register_dustvelocity')
      if (ldustdensity) call put_shared_variable('deltamd',deltamd)
!
    endsubroutine register_dustvelocity
!***********************************************************************
    subroutine initialize_dustvelocity(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  18-mar-03/axel+anders: adapted from hydro
!
      use EquationOfState, only: cs0
      use BorderProfiles, only: request_border_driving
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      integer :: i, j, k
      real :: gsurften, Eyoung, nu_Poisson, Eyoungred
!
!  Copy boundary condition on first dust species to all others.
!
      call copy_bcs_dust
!
!  Output grain mass discretization type. ldustcoagulation must also
!  be turned on for pure condensation, because we need the k bins.
!
      if (lroot .and. ldustcoagulation .and. .not.ldustdensity) then
        if (lmdvar) then
          print*, 'initialize_dustvelocity: variable grain mass'
        else
          print*, 'initialize_dustvelocity: constant grain mass'
        endif
      endif
!
!  Calculate inverse of minimum gas friction time.
!
      if (tausgmin/=0.0) then
        tausg1max=1.0/tausgmin
        if (lroot) print*,'initialize_dustvelocity: '//'minimum gas friction time tausgmin=',tausgmin
      endif
!
!  Define inverse of limiting friction time for short friction time
!  approximation.
!
      if (shorttauslimit/=0.0) shorttaus1limit=1/shorttauslimit
!
!  Grain chemistry
!
      if (lroot) print*, 'initialize_dustvelocity: dust_chemistry = ', dust_chemistry
!
      select case (dust_chemistry)

      case ('nothing')

        gsurften   = 0.0
        Eyoung     = 1.0
        nu_Poisson = 0.0
        Eyoungred  = 1.0
        unit_md = 1.0
        mumon   = 1.0
        mmon    = 1.0

      case ('ice')
!
!  Surface tension and Young's modulus for sticking velocity
!
        gsurften   = 370. ! erg cm^-2
        Eyoung     = 7e10 ! dyn cm^-2
        nu_Poisson = 0.25 !
        Eyoungred  = Eyoung/(2*(1-nu_Poisson**2))

        mumon = 18.0
        mmon  = mumon*1.6733e-24
        unit_md = mmon
!
!  for the following few items, no action is needed
!
      case ('pscalar')
      case ('hat(om*t)')
      case ('cos(om*t)')
      case ('simplified')

      case default
        call fatal_error('initialize_dustvelocity','no valid dust chemistry specified')

      endselect

      mumon1=1/mumon
!
!  Constant used in determination of sticking velocity
!    (extra factor 2 from Dominik & Tielens, 1997, end of Sec. 3.2)
!
      ustcst = sqrt(2* 2*9.6 * gsurften**(5/3.) * Eyoungred**(-2/3.))
!
!  Dust physics parameters.
!  Note that md(1) = md0*(1+deltamd)/2, so md0 = [2/(1+deltamd)] md(1).
!  With md(1)=4/3.*pi*ad0**3*rhods, we have md0=8*pi*ad0**3*rhods/[3(1+deltamd)]
!  So in practice we want to set ad1.
!
      if (ad0/=0.) md0 = 4/3.*pi*ad0**3*rhods/unit_md
      if (ad1/=0.) md0 = 8*pi/(3*(1.+deltamd))*ad1**3*rhods
      if (lroot) print*,'recalculated: md0=',md0
!
!  Choice between different spacings.
!  First, linearly spaced radius bins:
!
      if (llin_radiusbins) then
        do k=1,ndustspec
          ad(k)=ad0+ad1*(k-1)
        enddo
        md=4/3.*pi*ad**3*rhods
        llog_massbins=.false.
!
!  Logarithmically spaced mass bins:
!  (Do we really need unit_md? When would it not be 1?)
!
      elseif (llog_massbins) then
        do k=1,ndustspec
          mdminus(k) = md0*deltamd**(k-1)
          mdplus(k)  = md0*deltamd**k
          md(k) = 0.5*(mdminus(k)+mdplus(k))
        enddo
        ad=(0.75*md*unit_md/(pi*rhods))**onethird
        llin_radiusbins=.false.
      endif
      if (lroot) print*,'initialize_dustvelocity: ad=',ad
!
!  Reinitialize dustvelocity
!  'all_to_first' = reinitialize heavier particles to lightest one.
!
      if (reinitialize_uud) then
        do j=1,ninit
          select case (inituud(j))
          case ('all_to_first')
            do k=2,ndustspec
            do i=0,2
              f(:,:,:,iuud(k)+i)=f(:,:,:,iuud(1)+i)
            enddo
            enddo
          endselect
        enddo
      endif
!
!  Calculate betad.
!  By default (betad0=0), use Stokes formula, where Fd=6pi*mu_ext*ad*u.
!  Here, mu_ext is the dynamic viscosity of the gas, which is normally nearly
!  constant (better so than the kinematic one, which is given in the code
!  and which is often chosen larger based on numerical considerations).
!
      select case (draglaw)
      case ('stokes_varmass')
        if (lroot) print*,'initialize_dustvelocity: draglaw=',draglaw
        if (betad0/=0) then
          betad=betad0*md**(-2./3.)
        elseif (mu_ext/=0) then
          betad=4.5*mu_ext/(rhods*ad**2)
        else
          call fatal_error('initialize_dustvelocity','no betad calculation for draglaw='//trim(draglaw))
        endif
        if (lroot) print*,'initialize_dustvelocity: betad=',betad
!
!  Do nothing by default.
!
      case default
        if (lroot) print*, 'initialize_dustvelocity: No betad calculation for draglaw='//trim(draglaw)
      endselect
!
!  Grain geometry
!
      select case (dust_geometry)

      case ('sphere')
        dimd1 = onethird
        if (lroot) print*, 'initialize_dustvelocity: dust geometry = sphere'
        call get_dustsurface
        call get_dustcrosssection
        surfmon = surfd(1)*(mmon/(md(1)*unit_md))**(1.-dimd1)

      case default
        call fatal_error('initialize_dustvelocity','no such dust_geometry: '//trim(dust_geometry))
      endselect
!
!  Auxiliary variables necessary for different drag laws
!
      if (ldragforce_dust) then
        select case (draglaw)

        case ('epstein_var')
          rhodsad1 = 1./(rhods*ad)

        case ('epstein_cst')
          do k=1,ndustspec
            tausd1(:,k) = 1.0/tausd(k)
          enddo
!
!  Do nothing by default.
!
        case default
          if (lroot) print*, 'initialize_dustvelocity: doing nothing for draglaw='//trim(draglaw)

        endselect
      endif
!
!  If *_all set, make all primordial *(:) = *_all
!
      if (nud_all /= 0.) then
        select case (viscd_law)
        case ('const')
!         if (lroot .and. ip<6) &
!             print*, 'initialize_dustvelocity: nud_all=',nud_all
          do k=1,ndustspec
            if (nud(k) == 0.) nud(k)=nud_all
          enddo
        case ('md_exponential')
          nud=nud_all*md**viscd_exponent
        case ('ad_exponential')
          nud=nud_all*(ad/adref_nud)**viscd_exponent
        case default
          call fatal_error('initialize_dustvelocity','no such viscd_law: '//trim(viscd_law))
        endselect
        if (lroot) print*, 'initialize_dustvelocity: nud=',nud
      endif
!
      if (betad_all /= 0.) then
        if (lroot .and. ip<6) print*, 'initialize_dustvelocity: betad_all=',betad_all
        do k=1,ndustspec
          if (betad(k) == 0.) betad(k) = betad_all
        enddo
      endif
!
      if (tausd_all /= 0.) then
        if (lroot .and. ip<6) print*, 'initialize_dustvelocity: tausd_all=',tausd_all
        do k=1,ndustspec
          if (tausd(k) == 0.) tausd(k) = tausd_all
        enddo
      endif
!
      if (beta_dPdr_dust/=0.0) then
        beta_dPdr_dust_scaled=beta_dPdr_dust*Omega/cs0
        if (lroot) print*, 'initialize_dustvelocity: Global pressure '// &
                           'gradient with beta_dPdr_dust=', beta_dPdr_dust
      endif
!
      do i=1,nvisc_max
        select case (iviscd(i))
        case ('nothing','')
        case ('simplified', '0')
          if (lroot) print*, 'Viscous force (dust): nud*del2ud'
          lviscd_simplified=.true.
        case ('nud-const')
          if (lroot) print*, 'Viscous force (dust): nud*(del2ud+graddivud/3+2Sd.glnnd)'
          lviscd_nud_const=.true.
        case ('shock','nud-shock')
          lviscd_shock=.true.
        case ('shock_simplified','nud-shock_simplified')
          lviscd_shock_simplified=.true.
        case ('hyper3_simplified','hyper3-simplified')
          if (lroot) print*, 'Viscous force (dust): nud*del6ud'
          lviscd_hyper3_simplified=.true.
        case ('hyper3_rhod_nud-const','hyper3-rhod-nud-const')
          if (lroot) print*, 'Viscous force (dust): mud/rhod*del6ud'
          lviscd_hyper3_rhod_nud_const=.true.
        case ('hyper3_nud-const','hyper3-nud-const')
          if (lroot) print*, 'Viscous force (dust): nud*(del6ud+S.glnnd)'
          lviscd_hyper3_nud_const=.true.
        case ('hyper3-cyl','hyper3_cyl','hyper3-sph','hyper3_sph')
          if (lroot) print*,'viscous force: nud_hyper3/pi^4 *(Deltav)^6/Deltaq^2'
          lviscd_hyper3_polar=.true.
       case ('hyper3-mesh')
          if (lroot) print*,'viscous force: nud_hyper3_mesh/pi^5 *(Deltav)^6/Deltaq'
          lviscd_hyper3_mesh=.true.
        case default
          call fatal_error('initialize_dustvelocity','no such iviscd: '//trim(iviscd(i)))
        endselect
      enddo

      if (ldust_pressure.and.dust_pressure_factor==0.0) &
          call fatal_error('initialize_dustvelocity','dust_pressure_factor should not be 0')

      select case (borderuud)
      case ('zero','0','initial-condition')
!
!  Tell the BorderProfiles module if we intend to use border driving, so
!  that the module can request the right pencils.
!
        call request_border_driving(borderuud)
      case ('nothing')
        if (lroot.and.ip<=5) print*,"initialize_dustvelocity: borderuud='nothing'"
!
      case default
        call fatal_error('initialize_dustvelocity','no such borderuud: '//trim(borderuud))
      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_dustvelocity
!***********************************************************************
    subroutine copy_bcs_dust
!
!  Copy boundary conditions on first dust species to all others
!
!  27-feb-04/anders: Copied from initialize_dustvelocity
!
!  Copy boundary conditions on first dust species to all species
!
    integer :: k
!
    if (ndustspec>1) then
!
      bcx(iudx) = bcx(iudx(1))
      bcx(iudy) = bcx(iudy(1))
      bcx(iudz) = bcx(iudz(1))
      bcx(ind)  = bcx(ind(1))
      bcy(iudx) = bcy(iudx(1))
      bcy(iudy) = bcy(iudy(1))
      bcy(iudz) = bcy(iudz(1))
      bcy(ind)  = bcy(ind(1))
      bcz(iudx) = bcz(iudx(1))
      bcz(iudy) = bcz(iudy(1))
      bcz(iudz) = bcz(iudz(1))
      bcz(ind)  = bcz(ind(1))
!
      if (lmdvar) then
        bcx(imd) = bcx(imd(1))
        bcy(imd) = bcy(imd(1))
        bcz(imd) = bcz(imd(1))
      endif

      do k=1,2
        bcx12(iudx,k) = bcx12(iudx(1),k)
        bcx12(iudy,k) = bcx12(iudy(1),k)
        bcx12(iudz,k) = bcx12(iudz(1),k)
        bcx12(ind,k)  = bcx12(ind(1),k)
!
        bcy12(iudx,k) = bcy12(iudx(1),k)
        bcy12(iudy,k) = bcy12(iudy(1),k)
        bcy12(iudz,k) = bcy12(iudz(1),k)
        bcy12(ind,k)  = bcy12(ind(1),k)
!
        bcz12(iudx,k) = bcz12(iudx(1),k)
        bcz12(iudy,k) = bcz12(iudy(1),k)
        bcz12(iudz,k) = bcz12(iudz(1),k)
        bcz12(ind,k)  = bcz12(ind(1),k)
!
        if (lmdvar) then
          bcy12(imd,k) = bcy12(imd(1),k)
          bcx12(imd,k) = bcx12(imd(1),k)
          bcz12(imd,k) = bcz12(imd(1),k)
        endif
      enddo
!
      if (lroot) print*, 'copy_bcs_dust: Copied bcs on first dust species to all others'
    endif
!
    endsubroutine copy_bcs_dust
!***********************************************************************
    subroutine init_uud(f)
!
!  initialise uud; called from start.f90
!
!  18-mar-03/axel+anders: adapted from hydro
!  21-jan-15/MR: changes for use for reference state.
!
      use Sub
      use Gravity
      use Initcond
      use InitialCondition, only: initial_condition_uud
      use EquationOfState, only: pressure_gradient,cs20
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: lnrho,rho,cs2,rhod,cp1tilde
      real :: eps,cs,eta_glnrho,v_Kepler
      integer :: j,k,l
      logical :: lnothing
      real, dimension(:,:), pointer :: reference_state
      real, dimension(:), pointer :: beta_glnrho_global,beta_glnrho_scaled
!
      call get_shared_variable('beta_glnrho_global',beta_glnrho_global,caller='init_uud')
      call get_shared_variable('beta_glnrho_scaled',beta_glnrho_scaled)

      if (lreference_state) call get_shared_variable('reference_state',reference_state)
!
!  inituud corresponds to different initializations of uud (called from start).
!
      lnothing=.false.
      do j=1,ninit
        select case (inituud(j))

        case ('nothing')
          if (lroot .and. .not. lnothing) print*, 'init_uud: nothing'
          lnothing=.true.
        case ('zero', '0')
          do k=1,ndustspec; f(:,:,:,iudx(k):iudz(k))=0.0; enddo
          if (lroot) print*,'init_uud: zero dust velocity'
!
        case ('constant')
          do l=1,mx
            f(l,:,:,iudx(1)) = uudx0
            f(l,:,:,iudy(1)) = uudy0
            f(l,:,:,iudz(1)) = uudz0
          enddo
!
        case ('firsttwo')
          do l=1,mx
            f(l,:,:,iudx(1)) = uudx0
            f(l,:,:,iudy(1)) = uudy0
            f(l,:,:,iudz(1)) = uudz0
            if (ndustspec>1) then
              f(l,:,:,iudx(min(2,ndustspec))) = uudx0*2.
              f(l,:,:,iudy(min(2,ndustspec))) = uudy0*2.
              f(l,:,:,iudz(min(2,ndustspec))) = uudz0*2.
            endif
          enddo
!
        case ('gaussian-noise')
          do k=1,ndustspec; call gaunoise(ampluud,f,iudx(k),iudz(k)); enddo
        case ('sinwave-phase')
          do k=1,ndustspec
            call sinwave_phase(f,iudx(k),ampl_udx,kx_uud,ky_uud,kz_uud,phase_udx)
            call sinwave_phase(f,iudy(k),ampl_udy,kx_uud,ky_uud,kz_uud,phase_udy)
            call sinwave_phase(f,iudz(k),ampl_udz,kx_uud,ky_uud,kz_uud,phase_udz)
          enddo
        case ('udx_sinx')
          do l=1,mx; f(l,:,:,iudx(1)) = ampluud*sin(kx_uud*x(l)); enddo
        case ('udy_siny')
          do m=1,my; f(:,m,:,iudy(1)) = ampluud*sin(ky_uud*y(m)); enddo
        case ('sinwave-z-x')
          if (lroot) print*, 'init_uud: sinwave-z-x, ampluud=', ampluud
          call sinwave(ampluud,f,iudz(1),kx=kx_uud)
        case ('udz_sinz')
          do n=1,mz; f(:,:,n,iudz(1)) = ampluud*sin(kz_uud*z(n)); enddo
        case ('linear-z')
          do n=1,mz; f(:,:,n,iudz(1)) = ampluud*z(n); enddo
        case ('udz_siny')
          do m=m1,m2
            f(:,m,:,iudz(1)) = f(:,m,:,iudz(1)) + ampluud*sin(ky_uud*y(m))
          enddo
        case ('udx_sinxsinysinz')
          do l=1,mx; do m=1,my; do n=1,mz
            f(l,m,n,iudx(1)) = ampluud*sin(kx_uud*x(l))*sin(ky_uud*y(m))*sin(kz_uud*z(n))
          enddo; enddo; enddo
        case ('udy_sinxsinysinz')
          do l=1,mx; do m=1,my; do n=1,mz
            f(l,m,n,iudy(1)) = ampluud*sin(kx_uud*x(l))*sin(ky_uud*y(m))*sin(kz_uud*z(n))
          enddo; enddo; enddo
        case ('udz_sinxsinysinz')
          do l=1,mx; do m=1,my; do n=1,mz
            f(l,m,n,iudz(1)) = ampluud*sin(kx_uud*x(l))*sin(ky_uud*y(m))*sin(kz_uud*z(n))
          enddo; enddo; enddo
        case ('follow_gas','follow-gas')
          do k=1,ndustspec
            f(:,:,:,iudx(k):iudz(k))=f(:,:,:,iux:iuz)
          enddo
        case ('terminal_vz')
          if (lroot) print*, 'init_uud: terminal velocity'
          do k=1,ndustspec
            do m=m1,m2
              do n=n1,n2
                if (ldensity_nolog) then
                  if (lreference_state) then
                    rho = f(l1:l2,m,n,irho)+reference_state(:,iref_rho)
                  else
                    rho = f(l1:l2,m,n,irho)
                  endif
                  lnrho = log(rho)
                else
                  lnrho = f(l1:l2,m,n,ilnrho)
                  rho = exp(lnrho)
                endif
                if (ldustdensity_log) then
                  rhod = exp(f(l1:l2,m,n,ilnnd(k)))*md(k)
                else
                  rhod = f(l1:l2,m,n,ind(k))*md(k)
                endif
                call pressure_gradient(f,cs2,cp1tilde)
                call get_stoppingtime(f(l1:l2,m,n,iudx(k):iudz(k)),f(l1:l2,m,n,iux:iuz),rho,cs2,rhod,k)
                f(l1:l2,m,n,iudz(k)) = f(l1:l2,m,n,iudz(k)) - tausd1(:,k)**(-1)*nu_epicycle**2*z(n)
              enddo
            enddo
          enddo

        case ('vshear_dust')
!
!  Vertical shear due to global pressure gradient and back-reaction drag force
!  from dust on gas.
!
          if (lroot) then
            print*, 'init_uud: vertical shear due to dust'
            if (any(beta_glnrho_scaled/=0.0)) then
              print*, 'init_uud: beta_glnrho_scaled=', beta_glnrho_scaled
            elseif (beta_dPdr_dust_scaled/=0.0) then
              print*, 'init_uud: beta_dPdr_dust_scaled=', beta_dPdr_dust_scaled
            endif
          endif

          if (ldensity_nolog) then
            if (ldustdensity_log) then
              eps=sum(exp(f(l1:l2,m1:m2,n1:n2,ilnnd(1))))/sum(f(l1:l2,m1:m2,n1:n2,irho))
            else
              eps=sum(f(l1:l2,m1:m2,n1:n2,ind(1)))/sum(f(l1:l2,m1:m2,n1:n2,ilnrho))
            endif
          else
            if (ldustdensity_log) then
              eps=sum(exp(f(l1:l2,m1:m2,n1:n2,ilnnd(1))))/sum(exp(f(l1:l2,m1:m2,n1:n2,ilnrho)))
            else
              eps=sum(f(l1:l2,m1:m2,n1:n2,ind(1)))/sum(exp(f(l1:l2,m1:m2,n1:n2,ilnrho)))
            endif
          endif

          if (lroot) print*, 'init_uud: average dust-to-gas ratio=', eps

          do l=l1,l2; do m=m1,m2; do n=n1,n2
            cs=sqrt(cs20)

            if (.not. lvshear_dust_global_eps) then
              if (ldensity_nolog) then
                if (ldustdensity_log) then
                  eps=exp(f(l,m,n,ilnnd(1)))/f(l,m,n,ilnrho)
                else
                  eps=f(l,m,n,ind(1))/f(l,m,n,ilnrho)
                endif
              else
                if (ldustdensity_log) then
                  eps=exp(f(l,m,n,ilnnd(1)))/exp(f(l,m,n,ilnrho))
                else
                  eps=f(l,m,n,ind(1))/exp(f(l,m,n,ilnrho))
                endif
              endif
            endif

            if (beta_glnrho_scaled(1)/=0.0) then
              f(l,m,n,iux) = f(l,m,n,iux) - cs20*beta_glnrho_scaled(1)*eps*tausd(1)/ &
                                            (1.0+2*eps+eps**2+(Omega*tausd(1))**2)

              f(l,m,n,iuy) = f(l,m,n,iuy) + cs20*beta_glnrho_scaled(1)*(1+eps+(Omega*tausd(1))**2)/ &
                                            (2*Omega*(1.0+2*eps+eps**2+(Omega*tausd(1))**2))

              f(l,m,n,iudx(1)) = f(l,m,n,iudx(1)) + cs20*beta_glnrho_scaled(1)*tausd(1)/ &
                                                    (1.0+2*eps+eps**2+(Omega*tausd(1))**2)

              f(l,m,n,iudy(1)) = f(l,m,n,iudy(1)) + cs20*beta_glnrho_scaled(1)*(1+eps)/ &
                                                    (2*Omega*(1.0+2*eps+eps**2+(Omega*tausd(1))**2))
            elseif (beta_dPdr_dust_scaled/=0.0) then
              f(l,m,n,iux) = f(l,m,n,iux) - cs20*beta_dPdr_dust_scaled*eps*tausd(1)/ &
                                            (1.0+2*eps+eps**2+(Omega*tausd(1))**2)

              f(l,m,n,iuy) = f(l,m,n,iuy) - cs20*beta_dPdr_dust_scaled*(eps+eps**2)/ &
                                            (2*Omega*(1.0+2*eps+eps**2+(Omega*tausd(1))**2))

              f(l,m,n,iudx(1)) = f(l,m,n,iudx(1)) + cs20*beta_dPdr_dust_scaled*tausd(1)/ &
                                                    (1.0+2*eps+eps**2+(Omega*tausd(1))**2)

              f(l,m,n,iudy(1)) = f(l,m,n,iudy(1)) - cs20*beta_dPdr_dust_scaled* &
                                                    (eps+eps**2+(Omega*tausd(1))**2)/ &
                                                    (2*Omega*(1.0+2*eps+eps**2+(Omega*tausd(1))**2))
            endif
          enddo; enddo; enddo
!
        case ('vshear_dust_pseudo')
!
!  Vertical shear due to pseudo Coriolis force
!
          if (lroot) then
            print*, 'init_uud: vertical shear due to dust (pseudo)'
            print*, 'init_uud: u0_gas_pseudo=', u0_gas_pseudo
          endif
          do l=l1,l2; do m=m1,m2; do n=n1,n2
            if (ldensity_nolog) then
              if (ldustdensity_log) then
                eps=exp(f(l,m,n,ilnnd(1)))/f(l,m,n,ilnrho)
              else
                eps=f(l,m,n,ind(1))/f(l,m,n,ilnrho)
              endif
            else
              if (ldustdensity_log) then
                eps=exp(f(l,m,n,ilnnd(1)))/exp(f(l,m,n,ilnrho))
              else
                eps=f(l,m,n,ind(1))/exp(f(l,m,n,ilnrho))
              endif
            endif
            f(l,m,n,iux) = f(l,m,n,iux) + &
                u0_gas_pseudo*(1.0 + Omega_pseudo*tausd(1))/(1.0 + eps + Omega_pseudo*tausd(1))
            f(l,m,n,iudx) = f(l,m,n,iudx) + u0_gas_pseudo/(1.0 + eps + Omega_pseudo*tausd(1))
          enddo; enddo; enddo
!
        case ('streaming')
!
!  Mode unstable to streaming instability (Youdin & Goodman 2005)
!
          eta_glnrho = -0.5*abs(beta_glnrho_global(1))*beta_glnrho_global(1)
          v_Kepler   =  1.0/abs(beta_glnrho_global(1))

          if (lroot) print*, 'init_uud: eta, vK=', eta_glnrho, v_Kepler
!
          if (ldensity_nolog) then
            if (ldustdensity_log) then
              eps=sum(exp(f(l1:l2,m1:m2,n1:n2,ilnnd(1))))/sum(f(l1:l2,m1:m2,n1:n2,irho))
            else
              eps=sum(f(l1:l2,m1:m2,n1:n2,ind(1)))/sum(f(l1:l2,m1:m2,n1:n2,ilnrho))
            endif
          else
            if (ldustdensity_log) then
              eps=sum(exp(f(l1:l2,m1:m2,n1:n2,ilnnd(1))))/sum(exp(f(l1:l2,m1:m2,n1:n2,ilnrho)))
            else
              eps=sum(f(l1:l2,m1:m2,n1:n2,ind(1)))/sum(exp(f(l1:l2,m1:m2,n1:n2,ilnrho)))
            endif
          endif
!
          do m=m1,m2; do n=n1,n2
!
            f(l1:l2,m,n,ind(1)) = 0.0*f(l1:l2,m,n,ind(1)) + &
                eps*ampluud*cos(kz_uud*z(n))*cos(kx_uud*x(l1:l2))
!
            f(l1:l2,m,n,ilnrho) = f(l1:l2,m,n,ilnrho) + ampluud* &
                ( real(coeff(7))*cos(kx_uud*x(l1:l2)) - &
                 aimag(coeff(7))*sin(kx_uud*x(l1:l2)))*cos(kz_uud*z(n))
!
            f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) + eta_glnrho*v_Kepler*ampluud* &
                ( real(coeff(4))*cos(kx_uud*x(l1:l2)) - &
                 aimag(coeff(4))*sin(kx_uud*x(l1:l2)))*cos(kz_uud*z(n))
!
            f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + eta_glnrho*v_Kepler*ampluud* &
                ( real(coeff(5))*cos(kx_uud*x(l1:l2)) - &
                 aimag(coeff(5))*sin(kx_uud*x(l1:l2)))*cos(kz_uud*z(n))
!
            f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) + eta_glnrho*v_Kepler*(-ampluud)* &
                (aimag(coeff(6))*cos(kx_uud*x(l1:l2)) + &
                  real(coeff(6))*sin(kx_uud*x(l1:l2)))*sin(kz_uud*z(n))
!
            f(l1:l2,m,n,iudx(1)) = f(l1:l2,m,n,iudx(1)) + eta_glnrho*v_Kepler*ampluud* &
                ( real(coeff(1))*cos(kx_uud*x(l1:l2)) - &
                 aimag(coeff(1))*sin(kx_uud*x(l1:l2)))*cos(kz_uud*z(n))
!
            f(l1:l2,m,n,iudy(1)) = f(l1:l2,m,n,iudy(1)) + eta_glnrho*v_Kepler*ampluud* &
                ( real(coeff(2))*cos(kx_uud*x(l1:l2)) - &
                 aimag(coeff(2))*sin(kx_uud*x(l1:l2)))*cos(kz_uud*z(n))
!
            f(l1:l2,m,n,iudz(1)) = f(l1:l2,m,n,iudz(1)) + eta_glnrho*v_Kepler*(-ampluud)* &
                (aimag(coeff(3))*cos(kx_uud*x(l1:l2)) + &
                  real(coeff(3))*sin(kx_uud*x(l1:l2)))*sin(kz_uud*z(n))
!
          enddo; enddo
!
        case default
          call fatal_error('init_uud','no such inituud: '//trim(inituud(j)))

        endselect
!
!  End loop over initial conditions
!
      enddo
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_uud(f)
!
    endsubroutine init_uud
!***********************************************************************
    subroutine pencil_criteria_dustvelocity()
!
!  All pencils that the Dustvelocity module depends on are specified here.
!
!  20-11-04/anders: coded
!
      if (.not. lchemistry) then
        lpenc_requested(i_uud)=.true.
        if (ladvection_dust.and..not.ldustvelocity_shorttausd) &
            lpenc_requested(i_udgud)=.true.
        if (ldustvelocity_shorttausd) then
          if (lgrav) lpenc_requested(i_gg)=.true.
          lpenc_requested(i_cs2)=.true.
          lpenc_requested(i_jxbr)=.true.
          lpenc_requested(i_glnrho)=.true.
        endif
        if (ldragforce_dust.and..not.ldustvelocity_shorttausd) &
            lpenc_requested(i_rhod)=.true.
        if (ldragforce_gas) lpenc_requested(i_rho1)=.true.
        if (ldragforce_dust) then
          lpenc_requested(i_uu)=.true.
          if (draglaw=='epstein_var') then
            lpenc_requested(i_cs2)=.true.
            lpenc_requested(i_rho)=.true.
          endif
        endif
        if (ldust_pressure) then
          lpenc_requested(i_cs2)=.true.
          lpenc_requested(i_glnnd)=.true.
        endif
        if (lviscd_nud_const .or. lviscd_hyper3_nud_const .and. &  !MR: logics?
            ldustdensity) then
          lpenc_requested(i_sdij)=.true.
          lpenc_requested(i_glnnd)=.true.
        endif
        if (lviscd_simplified .or. lviscd_nud_const) &
            lpenc_requested(i_del2ud)=.true.
        if (lviscd_shock) then
          lpenc_requested(i_divud)=.true.
          lpenc_requested(i_glnrhod)=.true.
          lpenc_requested(i_graddivud)=.true.
          lpenc_requested(i_shock)=.true.
          lpenc_requested(i_gshock)=.true.
        endif
        if (lviscd_shock_simplified) then
          lpenc_requested(i_divud)=.true.
          lpenc_requested(i_graddivud)=.true.
          lpenc_requested(i_shock)=.true.
          lpenc_requested(i_gshock)=.true.
        endif
        if (lviscd_hyper3_simplified .or. lviscd_hyper3_nud_const .or. &
            lviscd_hyper3_rhod_nud_const) &
            lpenc_requested(i_del6ud)=.true.
        if (lviscd_nud_const .or. lviscd_hyper3_nud_const) &
            lpenc_requested(i_sdglnnd)=.true.
        if (lviscd_nud_const) lpenc_requested(i_graddivud)=.true.
        if (lviscd_hyper3_rhod_nud_const) lpenc_requested(i_rhod)=.true.
        if (beta_dPdr_dust/=0.) lpenc_requested(i_cs2)=.true.
        if (lstokes_highspeed_corr) lpenc_requested(i_rho)=.true.
!
        lpenc_diagnos(i_uud)=.true.
        if (maxval(idiag_divud2m)/=0) lpenc_diagnos(i_divud)=.true.
        if (maxval(idiag_rdudmax)/=0 .or. maxval(idiag_rdudxm)/=0 .or. &
            maxval(idiag_rdudym)/=0 .or. maxval(idiag_rdudzm)/=0 .or. &
            maxval(idiag_rdudx2m)/=0) &
            lpenc_diagnos(i_rhod)=.true.
        if (maxval(idiag_udrms)/=0 .or. maxval(idiag_udmax)/=0 .or. &
            maxval(idiag_rdudmax)/=0 .or. maxval(idiag_ud2m)/=0 .or. &
            maxval(idiag_udm2)/=0 .or. idiag_ekintot_dust/=0) &
            lpenc_diagnos(i_ud2)=.true.
        if (maxval(idiag_odrms)/=0 .or. maxval(idiag_odmax)/=0 .or. &
            maxval(idiag_od2m)/=0) lpenc_diagnos(i_od2)=.true.
        if (maxval(idiag_oudm)/=0) lpenc_diagnos(i_oud)=.true.
!
      endif
!
    endsubroutine pencil_criteria_dustvelocity
!***********************************************************************
    subroutine pencil_interdep_dustvelocity(lpencil_in)
!
!  Interdependency among pencils provided by the Dustvelocity module
!  is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_ud2)) lpencil_in(i_uud)=.true.
      if (lpencil_in(i_divud)) lpencil_in(i_udij)=.true.
      if (lpencil_in(i_udgud)) then
        lpencil_in(i_uud)=.true.
        lpencil_in(i_udij)=.true.
      endif
      if (lpencil_in(i_ood)) lpencil_in(i_udij)=.true.
      if (lpencil_in(i_od2)) lpencil_in(i_ood)=.true.
      if (lpencil_in(i_oud)) then
        lpencil_in(i_uud)=.true.
        lpencil_in(i_ood)=.true.
      endif
      if (lpencil_in(i_sdij)) then
        if (lviscd_nud_const) then
          lpencil_in(i_udij)=.true.
          lpencil_in(i_divud)=.true.
        endif
      endif
!
    endsubroutine pencil_interdep_dustvelocity
!***********************************************************************
    subroutine calc_pencils_dustvelocity(f,p)
!
!  Calculate Dustvelocity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  13-nov-04/anders: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      real, dimension (nx,3,3) :: tmp_pencil_3x3
      integer :: i,j,k
      real, dimension (nx) :: dot2_tmp
!
      intent(in) :: f
      intent(inout) :: p
!
      do k=1,ndustspec
! uud
        if (lpencil(i_uud)) p%uud(:,:,k)=f(l1:l2,m,n,iudx(k):iudz(k))
! ud2
        if (lpencil(i_ud2)) then
          call dot2_mn(p%uud(:,:,k),dot2_tmp)
          p%ud2(:,k)=dot2_tmp
        endif
! udij
        if (lpencil(i_udij)) call gij(f,iuud(k),p%udij(:,:,:,k),1)
! divud
        if (lpencil(i_divud)) p%divud(:,k) = p%udij(:,1,1,k) + p%udij(:,2,2,k) + p%udij(:,3,3,k)
! udgud
        if (lpencil(i_udgud)) then
          if (lspherical_coords.or.lcylindrical_coords) then
            call u_dot_grad(f,iuud(k),p%udij(:,:,:,k),p%uud(:,:,k),p%udgud(:,:,k))
          else
            call multmv_mn(p%udij(:,:,:,k),p%uud(:,:,k),p%udgud(:,:,k))
          endif
        endif
! ood
        if (lpencil(i_ood)) then
          p%ood(:,1,k)=p%udij(:,3,2,k)-p%udij(:,2,3,k)
          p%ood(:,2,k)=p%udij(:,1,3,k)-p%udij(:,3,1,k)
          p%ood(:,3,k)=p%udij(:,2,1,k)-p%udij(:,1,2,k)
        endif
! od2
        if (lpencil(i_od2)) call dot2_mn(p%ood(:,:,k),p%od2(:,k))
! oud
        if (lpencil(i_oud)) call dot_mn(p%ood(:,:,k),p%uud(:,:,k),p%oud(:,k))
! sdij
        if (lpencil(i_sdij)) then
          if (lviscd_nud_const) then
            do j=1,3
              p%sdij(:,j,j,k)=p%udij(:,j,j,k)
              do i=j+1,3
                p%sdij(:,i,j,k)=0.5*(p%udij(:,i,j,k)+p%udij(:,j,i,k))
                p%sdij(:,j,i,k)=p%sdij(:,i,j,k)
              enddo
              p%sdij(:,j,j,k)=p%sdij(:,j,j,k)-(1/3.0)*p%divud(:,k)
            enddo
          elseif (lviscd_hyper3_nud_const) then
            call gij(f,iuud(k),tmp_pencil_3x3,5)
            do i=1,3
              do j=1,3
                p%sdij(:,i,j,k)=tmp_pencil_3x3(:,i,j)
              enddo
            enddo
!          else
!              call warning('calc_pencils_dustvelocity','No rate-of-strain tensor matches iviscd=', iviscd
          endif
        endif
! del2ud
        if (lpencil(i_del2ud)) call del2v(f,iuud(k),p%del2ud(:,:,k))
! del6ud
        if (lpencil(i_del6ud)) call del6v(f,iuud(k),p%del6ud(:,:,k))
! graddivud
        if (lpencil(i_graddivud)) call del2v_etc(f,iuud(k),GRADDIV=p%graddivud(:,:,k))
      enddo
!
    endsubroutine calc_pencils_dustvelocity
!***********************************************************************
    subroutine duud_dt(f,df,p)

!  Dust velocity evolution
!  Calculate duud/dt = - uud.graduud - 2Omega x uud - 1/tausd*(uud-uu)
!
!  18-mar-03/axel+anders: adapted from hydro
!
      use Debug_IO
      use General
      use Sub
      use Deriv, only: der6
      use Diagnostics, only: max_mn_name
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: fviscd, AA_sfta, BB_sfta, tmp, tmp2
      real, dimension (nx) :: tausg1, mudrhod1, tmp3
      real :: c2, s2
      integer :: i, j, k, ju
!
      intent(in) :: f, p
      intent(out) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'duud_dt: SOLVE duud_dt'
      if (headtt) then
        call identify_bcs('udx',iudx(1))
        call identify_bcs('udy',iudy(1))
        call identify_bcs('udz',iudz(1))
      endif
!
!  Loop over dust species.
!
      do k=1,ndustspec
!
!  Inverse friction time pencil tausp1 is set in separate subroutine.
!
        call get_stoppingtime(p%uud(:,:,k),p%uu,p%rho,p%cs2,p%rhod(:,k),k)
!
!  Short stopping time approximation.
!  Calculated from master equation d(wx-ux)/dt = A + B*(wx-ux) = 0.
!
        if (ldustvelocity_shorttausd .and. any(tausd1(:,k)>=shorttaus1limit)) then
          do j=1,3
            if (lgrav) then
              AA_sfta(:,j)=p%gg(:,j)
            else
              AA_sfta(:,j)=0.0
            endif
          enddo
          if (ldensity) then
            do j=1,3; AA_sfta(:,j)=AA_sfta(:,j)+p%cs2(:)*p%glnrho(:,j); enddo
          endif
          if (lgrav) then
            if (lgravx_gas .neqv. lgravx_dust) then
              if (lgravx_gas) AA_sfta(:,1)=AA_sfta(:,1)-p%gg(:,1)
              if (lgravx_dust) AA_sfta(:,1)=AA_sfta(:,1)+p%gg(:,1)
            endif
            if (lgravz_gas .neqv. lgravz_dust) then
              if (lgravz_gas) AA_sfta(:,3)=AA_sfta(:,3)-p%gg(:,3)
              if (lgravz_dust) AA_sfta(:,3)=AA_sfta(:,3)+p%gg(:,3)
            endif
          endif
          if (lmagnetic) AA_sfta=AA_sfta-p%JxBr
          do j=1,3; BB_sfta(:,j)=-tausd1(:,k); enddo
          df(l1:l2,m,n,iudx(k):iudz(k)) = 1/dt_beta_ts(itsub)*( &
              f(l1:l2,m,n,iux:iuz)-f(l1:l2,m,n,iudx(k):iudz(k))-AA_sfta/BB_sfta)
        else
!
!  Direct integration of equation of motion.
!
          if (ladvection_dust) df(l1:l2,m,n,iudx(k):iudz(k)) = &
                               df(l1:l2,m,n,iudx(k):iudz(k)) - p%udgud(:,:,k)
!
!  Coriolis force, -2*Omega x ud
!  Omega=(-sin_theta, 0, cos_theta)
!  theta corresponds to latitude
!
          if (Omega/=0. .and. lcoriolisforce_dust) then
            if (theta==0) then
              if (headtt .and. k == 1) print*,'duud_dt: add Coriolis force; Omega=',Omega
              c2=2*Omega
              df(l1:l2,m,n,iudx(k)) = df(l1:l2,m,n,iudx(k)) + c2*p%uud(:,2,k)
              df(l1:l2,m,n,iudy(k)) = df(l1:l2,m,n,iudy(k)) - c2*p%uud(:,1,k)
            else
              if (headtt .and. k == 1) print*, 'duud_dt: Coriolis force; Omega,theta=',Omega,theta
              c2=2*Omega*cos(theta*pi/180.)
              s2=2*Omega*sin(theta*pi/180.)
              df(l1:l2,m,n,iudx(k)) = df(l1:l2,m,n,iudx(k)) + c2*p%uud(:,2,k)
              df(l1:l2,m,n,iudy(k)) = df(l1:l2,m,n,iudy(k)) - c2*p%uud(:,1,k) + s2*p%uud(:,3,k)
              df(l1:l2,m,n,iudz(k)) = df(l1:l2,m,n,iudz(k))                   + s2*p%uud(:,2,k)
            endif
          endif
!
!  Add drag force on dust
!
          if (ldragforce_dust) then
            do i=1,3
              df(l1:l2,m,n,iudx(k)-1+i)=df(l1:l2,m,n,iudx(k)-1+i) - tausd1(:,k)*(p%uud(:,i,k)-p%uu(:,i))
            enddo
!
!  Add drag force on gas (back-reaction from dust)
!
            if (ldragforce_gas) then
              tausg1 = p%rhod(:,k)*tausd1(:,k)*p%rho1
              if (tausgmin/=0.0) where (tausg1>=tausg1max) tausg1=tausg1max
              do i=1,3
                df(l1:l2,m,n,iux-1+i) = df(l1:l2,m,n,iux-1+i) - tausg1*(p%uu(:,i)-p%uud(:,i,k))
              enddo
              if (lfirst.and.ldt) dt1_max=max(dt1_max,(tausg1+tausd1(:,k))/cdtd)
            else
              if (lfirst.and.ldt) dt1_max=max(dt1_max,tausd1(:,k)/cdtd)
            endif
          endif
!
! Gravity force on dust in x direction
!
          if (gravx_dust/=0.0) df(l1:l2,m,n,iudx(k)) = df(l1:l2,m,n,iudx(k)) + gravx_dust
!
!  Add constant background pressure gradient beta=alpha*H0/r0, where alpha
!  comes from a global pressure gradient P = P0*(r/r0)^alpha.
!  (the term must be added to the dust equation of motion when measuring
!  velocities relative to the shear flow modified by the global pressure grad.)
!
          if (beta_dPdr_dust/=0.0) df(l1:l2,m,n,iudx(k)) = &
              df(l1:l2,m,n,iudx(k)) + p%cs2*beta_dPdr_dust_scaled
!
!  Artificial pressure force
!
          if (ldust_pressure) then
            do i=1,3
              df(l1:l2,m,n,iudx(k)-1+i) = df(l1:l2,m,n,iudx(k)-1+i) - &
                                          dust_pressure_factor*p%cs2*p%glnrho(:,i)
                  !dust_pressure_factor*p%cs2*p%glnnd(:,i,k)
            enddo
          endif
!
!  Add pseudo Coriolis force (to drive velocity difference between dust and gas)
!
          if (Omega_pseudo/=0.0) then
            df(l1:l2,m,n,iux)  = df(l1:l2,m,n,iux)  - Omega_pseudo*(p%uu(:,1)-u0_gas_pseudo)
            df(l1:l2,m,n,iudx) = df(l1:l2,m,n,iudx) - Omega_pseudo*p%uud(:,1,:)
          endif
!
!  Add viscosity on dust
!
          fviscd=0.0
          diffus_nud=0.0
          diffus_nud3=0.0
!
!  Viscous force: nud*del2ud
!     -- not physically correct (no momentum conservation)
!
          if (lviscd_simplified) then
            fviscd = fviscd + nud(k)*p%del2ud(:,:,k)
            if (lfirst.and.ldt) diffus_nud=diffus_nud+nud(k)*dxyz_2
          endif
!
!  Viscous force: nud*(del2ud+graddivud/3+2Sd.glnnd)
!    -- the correct expression for nud=const
!
          if (lviscd_nud_const) then
            if (ldustdensity) then
              fviscd = fviscd + 2*nud(k)*p%sdglnnd(:,:,k) + &
                       nud(k)*(p%del2ud(:,:,k)+1/3.0*p%graddivud(:,:,k))
            else
              fviscd = fviscd + nud(k)*(p%del2ud(:,:,k)+1/3.*p%graddivud(:,:,k))
            endif
            if (lfirst.and.ldt) diffus_nud=diffus_nud+nud(k)*dxyz_2
          endif
!
!  Viscous force: nud_shock
!
          if (lviscd_shock) then
            if (ldustdensity) then
              call multsv(p%divud(:,k),p%glnrhod(:,:,k),tmp2)
              tmp = tmp2 + p%graddivud(:,:,k)
            else
              tmp = p%graddivud(:,:,k)
            endif
            call multsv(nud_shock(k)*p%shock,tmp,tmp2)
            call multsv_add(tmp2,nud_shock(k)*p%divud(:,k),p%gshock,tmp)
            fviscd = fviscd + tmp
            if (lfirst.and.ldt) diffus_nud=diffus_nud+nud_shock(k)*p%shock*dxyz_2
          endif
!
!  Viscous force: nud_shock simplified (not momentum conserving)
!
          if (lviscd_shock_simplified) then
            tmp = p%graddivud(:,:,k)
            call multsv(nud_shock(k)*p%shock,tmp,tmp2)
            call multsv_add(tmp2,nud_shock(k)*p%divud(:,k),p%gshock,tmp)
            fviscd = fviscd + tmp
            if (lfirst.and.ldt) diffus_nud=diffus_nud+nud_shock(k)*p%shock*dxyz_2
          endif
!
!  Viscous force: nud*del6ud (not momentum-conserving)
!
          if (lviscd_hyper3_simplified) then
            fviscd = fviscd + nud_hyper3(k)*p%del6ud(:,:,k)
            if (lfirst.and.ldt) diffus_nud3=diffus_nud3+nud_hyper3(k)*dxyz_6
          endif
!
!  Viscous force: polar coordinates
!
          if (lviscd_hyper3_polar) then
            do j=1,3
              ju=j+iuud(k)-1
              do i=1,3
                call der6(f,ju,tmp3,i,IGNOREDX=.true.)
                fviscd(:,j) = fviscd(:,j) + nud_hyper3(k)*pi4_1*tmp3*dline_1(:,i)**2
              enddo
              if (lfirst.and.ldt) diffus_nud3=diffus_nud3+nud_hyper3(k)*pi4_1*dxmin_pencil**4
            enddo
          endif
!
!  Viscous force: Axel's mesh formulation
!
          if (lviscd_hyper3_mesh) then
            do j=1,3
              ju=j+iuud(k)-1
              do i=1,3
                call der6(f,ju,tmp3,i,IGNOREDX=.true.)
                fviscd(:,j) = fviscd(:,j) + nud_hyper3_mesh(k)*pi5_1/60.*tmp3*dline_1(:,i)
              enddo
            enddo
            if (lfirst .and. ldt) then
               advec_hypermesh_uud=nud_hyper3_mesh(k)*pi5_1*sqrt(dxyz_2)
               advec2_hypermesh=advec2_hypermesh+advec_hypermesh_uud**2
             endif
          endif
!
!  Viscous force: mud/rhod*del6ud
!
          if (lviscd_hyper3_rhod_nud_const) then
            mudrhod1=(nud_hyper3(k)*nd0*md0)/p%rhod(:,k)   ! = mud/rhod
            do i=1,3
              fviscd(:,i) = fviscd(:,i) + mudrhod1*p%del6ud(:,i,k)
            enddo
            if (lfirst.and.ldt) diffus_nud3=diffus_nud3+nud_hyper3(k)*dxyz_6
          endif
          if (lfirst.and.ldt) then
            maxdiffus3=max(maxdiffus3,diffus_nud3)
            maxdiffus=max(maxdiffus,diffus_nud)
          endif
!
!  Viscous force: nud*(del6ud+S.glnnd), where S_ij=d^5 ud_i/dx_j^5
!
          if (lviscd_hyper3_nud_const) then
            fviscd = fviscd + nud_hyper3(k)*(p%del6ud(:,:,k)+p%sdglnnd(:,:,k))
            if (lfirst.and.ldt) diffus_nud3=diffus_nud3+nud_hyper3(k)*dxyz_6
          endif
!
!  Add to dust equation of motion.
!
          df(l1:l2,m,n,iudx(k):iudz(k)) = df(l1:l2,m,n,iudx(k):iudz(k)) + fviscd
!
!  ``uud/dx'' for timestep
!
          if (lfirst .and. ldt) then
            advec_uud=sum(abs(p%uud(:,:,k))*dline_1,2)
            maxadvec=maxadvec+advec_uud
            if ((headtt.or.ldebug) .and. (ip<6)) then
              print*,'duud_dt: max(advec_uud) =',maxval(advec_uud)
              print*,'duud_dt: max(diffus_nud) =',maxval(diffus_nud)
            endif
          endif
!
!  Short friction time switch.
!
        endif
!
!  Apply border profile
!
        if (lborder_profiles) call set_border_dustvelocity(f,df,p,k)
!
!  End loop over dust species
!
      enddo

      call calc_diagnostics_dustvelocity(p)
!
    endsubroutine duud_dt
!***********************************************************************
    subroutine calc_diagnostics_dustvelocity(p)

      use Diagnostics
!
      type (pencil_case) :: p
!
      integer :: k
!
!  Calculate diagnostic variables
!
      if (ldiagnos) then
        do k=1,ndustspec
          if ((headtt.or.ldebug) .and. (ip<6)) print*, 'duud_dt: Calculate diagnostic values...'
          call sum_mn_name(p%ud2(:,k),idiag_udrms(k),lsqrt=.true.)
          call max_mn_name(p%ud2(:,k),idiag_udmax(k),lsqrt=.true.)
          if (idiag_rdudmax(k)/=0) &
              call max_mn_name(p%rhod(:,k)**2*p%ud2(:,k),idiag_rdudmax(k),lsqrt=.true.)
          call sum_mn_name(p%ud2(:,k),idiag_ud2m(k))
          call integrate_mn_name(p%ud2(:,k),idiag_ekintot_dust)
          call sum_mn_name(p%uud(:,1,k),idiag_udxm(k))
          call sum_mn_name(p%uud(:,2,k),idiag_udym(k))
          call sum_mn_name(p%uud(:,3,k),idiag_udzm(k))
          if (idiag_udx2m(k)/=0) call sum_mn_name(p%uud(:,1,k)**2,idiag_udx2m(k))
          if (idiag_udy2m(k)/=0) call sum_mn_name(p%uud(:,2,k)**2,idiag_udy2m(k))
          if (idiag_udz2m(k)/=0) call sum_mn_name(p%uud(:,3,k)**2,idiag_udz2m(k))
          call max_mn_name(p%ud2(:,k),idiag_udm2(k))
          if (idiag_divud2m(k)/=0) call sum_mn_name(p%divud(:,k)**2,idiag_divud2m(k))
          if (idiag_rdudxm(k)/=0) call sum_mn_name(p%rhod(:,k)*p%uud(:,1,k),idiag_rdudxm(k))
          if (idiag_rdudym(k)/=0) call sum_mn_name(p%rhod(:,k)*p%uud(:,2,k),idiag_rdudym(k))
          if (idiag_rdudzm(k)/=0) call sum_mn_name(p%rhod(:,k)*p%uud(:,3,k),idiag_rdudzm(k))
          if (idiag_rdudx2m(k)/=0) call sum_mn_name((p%rhod(:,k)*p%uud(:,1,k))**2,idiag_rdudx2m(k))
          call sum_mn_name(p%od2(:,k),idiag_odrms(k),lsqrt=.true.)
          call max_mn_name(p%od2(:,k),idiag_odmax(k),lsqrt=.true.)
          call sum_mn_name(p%od2(:,k),idiag_od2m(k))
          call sum_mn_name(p%oud(:,k),idiag_oudm(k))
          if (lfirst .and. ldt) then
            if (idiag_dtud(k)/=0) call max_mn_name(advec_uud/cdt,idiag_dtud(k),l_dt=.true.)
            if (idiag_dtnud(k)/=0) call max_mn_name(diffus_nud/cdtv,idiag_dtnud(k),l_dt=.true.)
          endif
        enddo
      endif
!
!  xy-averages
!
      if (l1davgfirst) then
        do k=1,ndustspec
          call xysum_mn_name_z(p%uud(:,1,k),idiag_udxmz(k))
          call xysum_mn_name_z(p%uud(:,2,k),idiag_udymz(k))
          call xysum_mn_name_z(p%uud(:,3,k),idiag_udzmz(k))
          if (idiag_udx2mz(k)/=0) call xysum_mn_name_z(p%uud(:,1,k)**2,idiag_udx2mz(k))
          if (idiag_udy2mz(k)/=0) call xysum_mn_name_z(p%uud(:,2,k)**2,idiag_udy2mz(k))
          if (idiag_udz2mz(k)/=0) call xysum_mn_name_z(p%uud(:,3,k)**2,idiag_udz2mz(k))
        enddo
      endif
!
!  2-D averages.
!
      if (l2davgfirst) then
        do k=1,ndustspec
          call zsum_mn_name_xy(p%uud(:,1,k),idiag_udxmxy(k))
          call zsum_mn_name_xy(p%uud(:,:,k),idiag_udymxy(k),(/0,1,0/))
          call zsum_mn_name_xy(p%uud(:,:,k),idiag_udzmxy(k),(/0,0,1/))
        enddo
      endif
!
    endsubroutine calc_diagnostics_dustvelocity
!***********************************************************************
    subroutine set_border_dustvelocity(f,df,p,k)
!
!  Calculates the driving term for the border profile
!  of the uud variable.
!
!  09-may-12/wlad: coded
!
      use BorderProfiles,  only: border_driving,set_border_initcond
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: f_target
      integer :: j,ju,k
!
      select case (borderuud)
!
      case ('zero','0')
        f_target=0.
!
      case ('initial-condition')
        do j=1,3
          ju=j+iuud(k)-1
          call set_border_initcond(f,ju,f_target(:,j))
        enddo
!
      case ('nothing')
        return
      endselect
!
      do j=1,3
        ju=j+iuud(k)-1
        call border_driving(f,df,p,f_target(:,j),ju)
      enddo
!
    endsubroutine set_border_dustvelocity
!***********************************************************************
    subroutine get_dustsurface
!
!  Calculate surface of dust particles
!
    integer :: i
!
      ad(1)    = (0.75*md(1)*unit_md/(pi*rhods))**dimd1
      surfd(1) = 4*pi*ad(1)**2
      do i=2,ndustspec
        ad(i)  = ad(1)*(md(i)/md(1))**dimd1
        surfd(i) = surfd(1)*(md(i)/md(1))**(1.-dimd1)
      enddo
!
    endsubroutine get_dustsurface
!***********************************************************************
    subroutine get_dustcrosssection
!
!  Calculate surface of dust particles. Add collision efficiency.
!
!  12-dec-14/Xiangyu: coded
!   5-mar-15/nils+Xiangyu+axel: used actual interval size for ratio. radius
!
      integer, parameter :: max_rows = 200, max_cols = 110
      
      integer :: i,j, row,col,ex=1,ey=1
      real, dimension(max_rows,max_cols) :: efficiency
      real, dimension(max_cols) :: radius
      real, dimension(max_rows) :: ratio
      real ::  e,radius_ratio, adx, ady, radiusx, radiusy,  ratiox, ratioy

      logical :: luse_table
!
!  select efficiency type
!
      select case (iefficiency_type)
!
      case ('nothing')
        luse_table=.false.
        efficiency=1.
      case ('read_table')
        luse_table=.true.
!
!  read file (can make more general)
!
        open(unit=11, file="Interpolation.txt")
        open(unit=12, file="radius.txt")
        open(unit=13, file="ratio.txt")
!
!  read table
!
        do row = 1,max_rows
          read(11,*) (efficiency(row,col),col=1,max_cols)
        enddo
        close(unit=11)
!
!  read corresponding radius
!
        do row = 1,max_cols
          read(12,*) radius(row)
        enddo
        close(unit=12)
!
!  read corresponding radius ratio
!
        do row = 1,max_rows
          read(13,*) ratio(row)
        enddo     
        close(unit=13)
      endselect
!
!  compute cross section
!
      do i=1,ndustspec
        do j=1,ndustspec
          adx = ad(i)
          ady = ad(j)
!
          if (luse_table) then
            if (adx<=ady) then
              adx=ad(j)
              ady=ad(i)
            endif
            radius_ratio = ady/adx
!
!  radius within interval
!
            do col = 2,max_cols-1
              radiusx = .5*(radius(col)+radius(col+1))
              radiusy = .5*(radius(col)+radius(col-1))
              if (adx<=radiusy .and. adx>radiusx) ey = col
            end do
!
!  lower end
!
            col = 1
            radiusx = .5*(radius(col)+radius(col+1))
            if (adx>radiusx) ey = col
!
!  upper end
!
            col = max_cols
            radiusy = .5*(radius(col)+radius(col-1))
            if (adx<=radiusy) ey = col
!
!  ratio within interval
!
            do row = 2,max_rows-1
              ratiox = .5*(ratio(row)+ratio(row+1))
              ratioy = .5*(ratio(row)+ratio(row-1))
              if (radius_ratio>=ratioy .and. radius_ratio<ratiox) ex = row
            end do
!
!  lower end
!
            row = 1
            ratiox = .5*(ratio(row)+ratio(row+1))
            if (radius_ratio<ratiox) ex = row
!
!  upper end
!
            row = max_rows
            ratioy = .5*(ratio(row)+ratio(row-1))
            if (radius_ratio>=ratioy) ex = row
!
!  set efficiency
!
            e = efficiency(ex,ey)
          else
            e=1.
          endif
          scolld(i,j) = e*pi*(adx+ady)**2
        enddo
      enddo
!
    endsubroutine get_dustcrosssection
!***********************************************************************
    subroutine get_stoppingtime(uud,uu,rho,cs2,rhod,k)
!
!  Calculate stopping time depending on choice of drag law.
!
      use Sub, only: dot2

      !real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: uu
      real, dimension (nx,3) :: uud
      real, dimension (nx) :: rho,rhod,csrho,cs2,deltaud2, Rep
      real :: pifactor1, pifactor2
      integer :: k
!
      select case (draglaw)

      case ('epstein_cst')
        ! Do nothing, initialized in initialize_dustvelocity
      case ('epstein_cst_b')
        tausd1(:,k) = betad(k)/rhod
      case ('stokes_cst_tausd')
        tausd1(:,k) = betad(k)
      case ('stokes_varmass')
        tausd1(:,k) = betad(k)
!
!  Correction term that account also for larger particle Reynolds
!  numbers (see e.g. Haugen and Kragset 2010)
!
        if (lstokes_highspeed_corr) then
          call dot2(uud-uu,deltaud2)
          Rep=2*ad(k)*rho*sqrt(deltaud2)/mu_ext
          tausd1(:,k) = tausd1(:,k)*(1+0.15*Rep**0.687)
        endif
!
!  Allow for a sqrt(8/pi) factor (=pifactor1) and a
!  9*pi/128 factor (=pifactor2) from Mattsson+Fynbo+Villarroel19
!
      case ('epstein_var')
        call dot2(uud-uu,deltaud2)
        if (lpifactor1) then
          pifactor1=sqrt(8./pi)
        else
          pifactor1=1.
        endif
        if (lpifactor2) then
          pifactor2=9.*pi/128.
        else
          pifactor2=1.
        endif
        csrho=pifactor1*sqrt(cs2+pifactor2*deltaud2)*rho
        tausd1(:,k) = csrho*rhodsad1(k)
      case ('epstein_gaussian_z')
        tausd1(:,k) = (1/tausd(k))*exp(-z(n)**2/(2*scaleHtaus**2))
        if (z0taus/=0.0) tausd1(:,k)=tausd1(:,k)/( &
            0.5*(tanh((z(n)+z0taus)/widthtaus)+tanh((-z(n)+z0taus)/widthtaus)))
      case default
        call fatal_error("get_stoppingtime","no such drag law: "//trim(draglaw))

      endselect
!
    endsubroutine get_stoppingtime
!***********************************************************************
    subroutine read_dustvelocity_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=dustvelocity_init_pars, IOSTAT=iostat)
!
    endsubroutine read_dustvelocity_init_pars
!***********************************************************************
    subroutine write_dustvelocity_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=dustvelocity_init_pars)
!
    endsubroutine write_dustvelocity_init_pars
!***********************************************************************
    subroutine read_dustvelocity_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=dustvelocity_run_pars, IOSTAT=iostat)
!
    endsubroutine read_dustvelocity_run_pars
!***********************************************************************
    subroutine write_dustvelocity_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=dustvelocity_run_pars)
!
    endsubroutine write_dustvelocity_run_pars
!***********************************************************************
    subroutine rprint_dustvelocity(lreset,lwrite)
!
!  Reads and registers print parameters relevant for dust velocity.
!
!   3-may-02/axel: coded
!
      use Diagnostics
      use General, only: itoa
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname, inamez, inamexy, k
      character (len=intlen) :: sdust
!
!  Reset everything in case of reset.
!
      if (lreset) then
        idiag_dtud=0; idiag_dtnud=0; idiag_ud2m=0; idiag_udx2m=0
        idiag_udxm=0; idiag_udym=0; idiag_udzm=0
        idiag_udy2m=0; idiag_udz2m=0; idiag_udm2=0; idiag_oudm=0; idiag_od2m=0
        idiag_udrms=0
        idiag_ekintot_dust=0
        idiag_udmax=0; idiag_odrms=0; idiag_odmax=0; idiag_rdudmax=0
        idiag_udmx=0; idiag_udmy=0; idiag_udmz=0; idiag_divud2m=0
        idiag_epsKd=0; idiag_rdudxm=0;idiag_rdudym=0; idiag_rdudzm=0;
        idiag_rdudx2m=0; idiag_udx2mz=0; idiag_udy2mz=0; idiag_udz2mz=0
      endif
!
!  Loop over dust layers
!
      do k=1,ndustspec
!
!  iname runs through all possible names that may be listed in print.in
!
        if (lroot.and.ip<14) print*,'rprint_dustvelocity: run through parse list'
        do iname=1,nname
          if (ndustspec == 1) then
            sdust=''
          else
            sdust=itoa(k)
          endif
          call parse_name(iname,cname(iname),cform(iname),'dtud'//trim(sdust),idiag_dtud(k))
          call parse_name(iname,cname(iname),cform(iname),'dtnud'//trim(sdust),idiag_dtnud(k))
          call parse_name(iname,cname(iname),cform(iname),'udxm'//trim(sdust),idiag_udxm(k))
          call parse_name(iname,cname(iname),cform(iname),'udym'//trim(sdust),idiag_udym(k))
          call parse_name(iname,cname(iname),cform(iname),'udzm'//trim(sdust),idiag_udzm(k))
          call parse_name(iname,cname(iname),cform(iname),'ud2m'//trim(sdust),idiag_ud2m(k))
          call parse_name(iname,cname(iname),cform(iname),'ekintot_dust'//trim(sdust),idiag_ekintot_dust)
          call parse_name(iname,cname(iname),cform(iname),'udx2m'//trim(sdust),idiag_udx2m(k))
          call parse_name(iname,cname(iname),cform(iname),'udy2m'//trim(sdust),idiag_udy2m(k))
          call parse_name(iname,cname(iname),cform(iname),'udz2m'//trim(sdust),idiag_udz2m(k))
          call parse_name(iname,cname(iname),cform(iname),'udm2'//trim(sdust),idiag_udm2(k))
          call parse_name(iname,cname(iname),cform(iname),'od2m'//trim(sdust),idiag_od2m(k))
          call parse_name(iname,cname(iname),cform(iname),'oudm'//trim(sdust),idiag_oudm(k))
          call parse_name(iname,cname(iname),cform(iname),'udrms'//trim(sdust),idiag_udrms(k))
          call parse_name(iname,cname(iname),cform(iname),'udmax'//trim(sdust),idiag_udmax(k))
          call parse_name(iname,cname(iname),cform(iname),'rdudmax'//trim(sdust),idiag_rdudmax(k))
          call parse_name(iname,cname(iname),cform(iname),'rdudxm'//trim(sdust),idiag_rdudxm(k))
          call parse_name(iname,cname(iname),cform(iname),'rdudym'//trim(sdust),idiag_rdudym(k))
          call parse_name(iname,cname(iname),cform(iname),'rdudzm'//trim(sdust),idiag_rdudzm(k))
          call parse_name(iname,cname(iname),cform(iname),'rdudx2m'//trim(sdust),idiag_rdudx2m(k))
          call parse_name(iname,cname(iname),cform(iname),'odrms'//trim(sdust),idiag_odrms(k))
          call parse_name(iname,cname(iname),cform(iname),'odmax'//trim(sdust),idiag_odmax(k))
          call parse_name(iname,cname(iname),cform(iname),'udmx'//trim(sdust),idiag_udmx(k))
          call parse_name(iname,cname(iname),cform(iname),'udmy'//trim(sdust),idiag_udmy(k))
          call parse_name(iname,cname(iname),cform(iname),'udmz'//trim(sdust),idiag_udmz(k))
          call parse_name(iname,cname(iname),cform(iname),'divud2m'//trim(sdust),idiag_divud2m(k))
          call parse_name(iname,cname(iname),cform(iname),'epsKd'//trim(sdust),idiag_epsKd(k))
        enddo
!
!  Check for those quantities for which we want xy-averages.
!
        do inamez=1,nnamez
          call parse_name(inamez,cnamez(inamez),cformz(inamez),'udxmz'//trim(sdust),idiag_udxmz(k))
          call parse_name(inamez,cnamez(inamez),cformz(inamez),'udymz'//trim(sdust),idiag_udymz(k))
          call parse_name(inamez,cnamez(inamez),cformz(inamez),'udzmz'//trim(sdust),idiag_udzmz(k))
          call parse_name(inamez,cnamez(inamez),cformz(inamez),'udx2mz'//trim(sdust),idiag_udx2mz(k))
          call parse_name(inamez,cnamez(inamez),cformz(inamez),'udy2mz'//trim(sdust),idiag_udy2mz(k))
          call parse_name(inamez,cnamez(inamez),cformz(inamez),'udz2mz'//trim(sdust),idiag_udz2mz(k))
        enddo
!
!  Check for those quantities for which we want z-averages.
!
        do inamexy=1,nnamexy
          call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'udxmxy'//trim(sdust),idiag_udxmxy(k))
          call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'udymxy'//trim(sdust),idiag_udymxy(k))
          call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'udzmxy'//trim(sdust),idiag_udzmxy(k))
        enddo
!
!  End loop over dust layers
!
      enddo
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then 
        where(cnamev=='uud') cformv='DEFINED'
      endif
!
    endsubroutine rprint_dustvelocity
!***********************************************************************
    subroutine get_slices_dustvelocity(f,slices)
!
!  Write slices for animation of Dustvelocity variables.
!
!  26-jul-06/tony: coded
!
      use Slices_methods, only: assign_slices_vec

      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Dustvelocity.
!
        case ('uud'); call assign_slices_vec(slices,f,iudx(1))
!
      endselect
!
    endsubroutine get_slices_dustvelocity
!***********************************************************************
endmodule Dustvelocity
