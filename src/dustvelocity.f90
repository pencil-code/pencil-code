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
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'dustvelocity.h'
!ajwm SHOULDN'T REALLY BE SHARED
!ajwm but are used consistently with the Dustdensity module
!ajwm - not good but for reasons of dust density / velocity interaction
  public :: dust_geometry, dimd1, rhods, surfd, mdplus, mdminus
  public :: ad, scolld, ustcst, tausd1, tausd
  public :: unit_md, dust_chemistry, mumon, mmon, mi, md
!
  complex, dimension (7) :: coeff
  real, dimension(ndustspec,ndustspec) :: scolld
  real, dimension(nx,ndustspec) :: tausd1
  real, dimension(ndustspec) :: md=1.0,mdplus,mdminus,ad,surfd,mi,rhodsad1
  real, dimension(ndustspec) :: tausd=1.,betad=0.,nud=0.,nud_hyper3=0.
  real :: ampluud=0.,ampl_udx=0.0,ampl_udy=0.0,ampl_udz=0.0
  real :: phase_udx=0.0, phase_udy=0.0, phase_udz=0.0
  real :: kx_uud=1.,ky_uud=1.,kz_uud=1.
  real :: rhods=1.,nd0=1.,md0=1.,rhod0=1.
  real :: ad0=0.,ad1=0.,dimd1=0.333333,deltamd=1.0
  real :: nud_all=0.,betad_all=0.,tausd_all=0.
  real :: mmon,mumon,mumon1,surfmon,ustcst,unit_md=1.0
  real :: beta_dPdr_dust=0.0, beta_dPdr_dust_scaled=0.0,cdtd=0.2
  real :: Omega_pseudo=0.0, u0_gas_pseudo=0.0, tausgmin=0.0, tausg1max=0.0
  real :: shorttauslimit=0.0, shorttaus1limit=0.0
  real :: scaleHtaus=1.0, z0taus=0.0, widthtaus=1.0
  logical :: ladvection_dust=.true.,lcoriolisforce_dust=.true.
  logical :: ldragforce_dust=.true.,ldragforce_gas=.false.
  logical :: lviscosity_dust=.true.
  logical :: ldustvelocity_shorttausd=.false., lvshear_dust_global_eps=.false.
  logical :: ldustcoagulation=.false., ldustcondensation=.false.
  character (len=labellen), dimension(ninit) :: inituud='nothing'
  character (len=labellen) :: draglaw='epstein_cst',iviscd='simplified'
  character (len=labellen) :: dust_geometry='sphere', dust_chemistry='nothing'
!
  namelist /dustvelocity_init_pars/ &
      ampl_udx, ampl_udy, ampl_udz, phase_udx, phase_udy, phase_udz, &
      rhods, md0, ad0, ad1, deltamd, draglaw, ampluud, inituud, &
      kx_uud, ky_uud, kz_uud, Omega_pseudo, u0_gas_pseudo, &
      dust_chemistry, dust_geometry, tausd, beta_dPdr_dust, coeff, &
      ldustcoagulation, ldustcondensation, lvshear_dust_global_eps, cdtd, &
      ldustvelocity_shorttausd, scaleHtaus, z0taus
!
  namelist /dustvelocity_run_pars/ &
      nud, nud_all, iviscd, betad, betad_all, tausd, tausd_all, draglaw, &
      ldragforce_dust, ldragforce_gas, ldustvelocity_shorttausd, &
      ladvection_dust, lcoriolisforce_dust, beta_dPdr_dust, tausgmin, cdtd, &
      nud_hyper3, scaleHtaus, z0taus, widthtaus, shorttauslimit
!
  integer :: idiag_ekintot_dust=0
  integer, dimension(ndustspec) :: idiag_ud2m=0
  integer, dimension(ndustspec) :: idiag_udxm=0,idiag_udym=0,idiag_udzm=0
  integer, dimension(ndustspec) :: idiag_udx2m=0,idiag_udy2m=0,idiag_udz2m=0
  integer, dimension(ndustspec) :: idiag_udm2=0,idiag_oudm=0,idiag_od2m=0
  integer, dimension(ndustspec) :: idiag_udrms=0,idiag_udmax=0,idiag_odrms=0
  integer, dimension(ndustspec) :: idiag_odmax=0,idiag_rdudmax=0
  integer, dimension(ndustspec) :: idiag_udxmz=0,idiag_udymz=0,idiag_udzmz=0
  integer, dimension(ndustspec) :: idiag_udx2mz=0,idiag_udy2mz=0,idiag_udz2mz=0
  integer, dimension(ndustspec) :: idiag_udmx=0,idiag_udmy=0,idiag_udmz=0
  integer, dimension(ndustspec) :: idiag_udxmxy=0,idiag_udymxy=0,idiag_udzmxy=0
  integer, dimension(ndustspec) :: idiag_divud2m=0,idiag_epsKd=0
  integer, dimension(ndustspec) :: idiag_dtud=0,idiag_dtnud=0
  integer, dimension(ndustspec) :: idiag_rdudxm=0,idiag_rdudym=0,idiag_rdudzm=0
  integer, dimension(ndustspec) :: idiag_rdudx2m=0
!
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
      use General, only: chn
!
      integer :: k, uud_tmp
      character(len=5) :: sdust
!
!  Write dust index in short notation
!
      if (lroot .and. ndustspec/=1) then
        open(3,file=trim(datadir)//'/index.pro', position='append')
        call chn(ndustspec,sdust)
        write(3,*) 'iuud=intarr('//trim(sdust)//')'
        close(3)
      endif
!
      do k=1,ndustspec
        call chn(k-1,sdust)
        sdust='['//trim(sdust)//']'
        if (ndustspec==1) sdust=''
        call farray_register_pde('uud'//sdust,uud_tmp,vector=3)
        iuud(k) = uud_tmp
        iudx(k) = iuud(k)
        iudy(k) = iuud(k)+1
        iudz(k) = iuud(k)+2
      enddo
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
          "$Id$")
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
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      integer :: k
      real :: gsurften,Eyoung,nu_Poisson,Eyoungred
!
!  Copy boundary condition on first dust species to all others.
!
      call copy_bcs_dust
!
!  Output grain mass discretization type
!
      if (lroot .and. ldustcoagulation) then
        if (lmdvar) then
          print*, 'register_dustvelocity: variable grain mass'
        else
          print*, 'register_dustvelocity: constant grain mass'
        endif
      endif
!
!  Turn off dust viscosity if zero viscosity
!
      if (maxval(nud) == 0.) lviscosity_dust=.false.
      if (lroot) print*, &
          'initialize_dustvelocity: lviscosity_dust=',lviscosity_dust
!
!  Calculate inverse of minimum gas friction time.
!
      if (tausgmin/=0.0) then
        tausg1max=1.0/tausgmin
        if (lroot) print*, 'initialize_dustvelocity: '// &
            'minimum gas friction time tausgmin=', tausgmin
      endif
!
!  Define inverse of limiting friction time for short friction time
!  approximation.
!
      if (shorttauslimit/=0.0) shorttaus1limit=1/shorttauslimit
!
      if (ldustcoagulation .or. ldustcondensation) then
!
!  Grain chemistry
!
        if (lroot) &
            print*, 'initialize_dustvelocity: dust_chemistry = ', dust_chemistry
!
        select case (dust_chemistry)

        case ('nothing')
          gsurften   = 0.
          Eyoung     = 1.
          nu_Poisson = 0.
          Eyoungred  = 1.
          unit_md = 1.
          mumon   = 1.
          mmon    = 1.

        case ('ice')
!
!  Surface tension and Young's modulus for sticking velocity
!
          gsurften   = 370. ! erg cm^-2
          Eyoung     = 7e10 ! dyn cm^-2
          nu_Poisson = 0.25 !
          Eyoungred  = Eyoung/(2*(1-nu_Poisson**2))

          mumon = 18.
          mmon  = mumon*1.6733e-24
          unit_md = mmon

        case default
          call fatal_error &
              ('initialize_dustvelocity','No valid dust chemistry specified.')

        endselect

        mumon1=1/mumon
!
!  Constant used in determination of sticking velocity
!    (extra factor 2 from Dominik & Tielens, 1997, end of Sec. 3.2)
!
        ustcst = sqrt(2* 2*9.6 * gsurften**(5/3.) * Eyoungred**(-2/3.))
!
!  Dust physics parameters
!
        if (ad0/=0.) md0 = 4/3.*pi*ad0**3*rhods/unit_md
        if (ad1/=0.) md0 = 8*pi/(3*(1.+deltamd))*ad1**3*rhods
!
!  Mass bins
!
        do k=1,ndustspec
          mdminus(k) = md0*deltamd**(k-1)
          mdplus(k)  = md0*deltamd**k
          md(k) = 0.5*(mdminus(k)+mdplus(k))
        enddo
!
!  Grain geometry
!
        select case (dust_geometry)

        case ('sphere')
          dimd1 = 0.333333
          if (lroot) print*, 'initialize_dustvelocity: dust geometry = sphere'
          call get_dustsurface
          call get_dustcrosssection
          surfmon = surfd(1)*(mmon/(md(1)*unit_md))**(1.-dimd1)

        case default
          call fatal_error( &
              'initialize_dustvelocity','No valid dust geometry specified.')

        endselect
      endif
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

        endselect
      endif
!
!  If *_all set, make all primordial *(:) = *_all
!
      if (nud_all /= 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustvelocity: nud_all=',nud_all
        do k=1,ndustspec
          if (nud(k) == 0.) nud(k)=nud_all
        enddo
      endif
!
      if (betad_all /= 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustvelocity: betad_all=',betad_all
        do k=1,ndustspec
          if (betad(k) == 0.) betad(k) = betad_all
        enddo
      endif
!
      if (tausd_all /= 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustvelocity: tausd_all=',tausd_all
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
!
!  Copy boundary conditions on first dust species to all species
!
      if (ndustspec>1) then
         bcx(iudx) =  bcx(iudx(1))
        bcx1(iudx) = bcx1(iudx(1))
        bcx2(iudx) = bcx2(iudx(1))
         bcx(iudy) =  bcx(iudy(1))
        bcx1(iudy) = bcx1(iudy(1))
        bcx2(iudy) = bcx2(iudy(1))
         bcx(iudz) =  bcx(iudz(1))
        bcx1(iudz) = bcx1(iudz(1))
        bcx2(iudz) = bcx2(iudz(1))
         bcx(ind)  =  bcx(ind(1))
        bcx1(ind)  = bcx1(ind(1))
        bcx2(ind)  = bcx2(ind(1))
        if (lmdvar) then
           bcx(imd) =  bcx(imd(1))
          bcx1(imd) = bcx1(imd(1))
          bcx2(imd) = bcx2(imd(1))
        endif
!
         bcy(iudx) =  bcy(iudx(1))
        bcy1(iudx) = bcy1(iudx(1))
        bcy2(iudx) = bcy2(iudx(1))
         bcy(iudy) =  bcy(iudy(1))
        bcy1(iudy) = bcy1(iudy(1))
        bcy2(iudy) = bcy2(iudy(1))
         bcy(iudz) =  bcy(iudz(1))
        bcy1(iudz) = bcy1(iudz(1))
        bcy2(iudz) = bcy2(iudz(1))
         bcy(ind)  =  bcy(ind(1))
        bcy1(ind)  = bcy1(ind(1))
        bcy2(ind)  = bcy2(ind(1))
        if (lmdvar) then
           bcy(imd) =  bcy(imd(1))
          bcy1(imd) = bcy1(imd(1))
          bcy2(imd) = bcy2(imd(1))
        endif
!
         bcz(iudx) =  bcz(iudx(1))
        bcz1(iudx) = bcz1(iudx(1))
        bcz2(iudx) = bcz2(iudx(1))
         bcz(iudy) =  bcz(iudy(1))
        bcz1(iudy) = bcz1(iudy(1))
        bcz2(iudy) = bcz2(iudy(1))
         bcz(iudz) =  bcz(iudz(1))
        bcz1(iudz) = bcz1(iudz(1))
        bcz2(iudz) = bcz2(iudz(1))
         bcz(ind)  =  bcz(ind(1))
        bcz1(ind)  = bcz1(ind(1))
        bcz2(ind)  = bcz2(ind(1))
        if (lmdvar) then
           bcz(imd) =  bcz(imd(1))
          bcz1(imd) = bcz1(imd(1))
          bcz2(imd) = bcz2(imd(1))
        endif
        if (lroot) print*, 'copy_bcs_dust: '// &
            'Copied bcs on first dust species to all others'
      endif
!
    endsubroutine copy_bcs_dust
!***********************************************************************
    subroutine init_uud(f)
!
!  initialise uud; called from start.f90
!
!  18-mar-03/axel+anders: adapted from hydro
!
      use EquationOfState, only: gamma, beta_glnrho_global, beta_glnrho_scaled
      use Sub
      use Gravity
      use Initcond
      use InitialCondition, only: initial_condition_uud
      use EquationOfState, only: pressure_gradient,cs20
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: lnrho,rho,cs2,rhod,cp1tilde
      real :: eps,cs,eta_glnrho,v_Kepler
      integer :: j,k,l
      logical :: lnothing
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
        case ('udz_siny')
          do m=m1,m2
            f(:,m,:,iudz(1)) = f(:,m,:,iudz(1)) + ampluud*sin(ky_uud*y(m))
          enddo
        case ('udx_sinxsinysinz')
          do l=1,mx; do m=1,my; do n=1,mz
            f(l,m,n,iudx(1)) = &
                ampluud*sin(kx_uud*x(l))*sin(ky_uud*y(m))*sin(kz_uud*z(n))
          enddo; enddo; enddo
        case ('udy_sinxsinysinz')
          do l=1,mx; do m=1,my; do n=1,mz
            f(l,m,n,iudy(1)) = &
                ampluud*sin(kx_uud*x(l))*sin(ky_uud*y(m))*sin(kz_uud*z(n))
          enddo; enddo; enddo
        case ('udz_sinxsinysinz')
          do l=1,mx; do m=1,my; do n=1,mz
            f(l,m,n,iudz(1)) = &
                ampluud*sin(kx_uud*x(l))*sin(ky_uud*y(m))*sin(kz_uud*z(n))
          enddo; enddo; enddo
        case ('follow_gas')
          do k=1,ndustspec
            f(:,:,:,iudx(k):iudz(k))=f(:,:,:,iux:iuz)
          enddo
        case ('terminal_vz')
          if (lroot) print*, 'init_uud: terminal velocity'
          do k=1,ndustspec
            do m=m1,m2
              do n=n1,n2
                if (ldensity_nolog) then
                  rho = f(l1:l2,m,n,irho)
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
                call get_stoppingtime(f,rho,cs2,rhod,k)
                f(l1:l2,m,n,iudz(k)) = &
                    f(l1:l2,m,n,iudz(k)) - tausd1(:,k)**(-1)*nu_epicycle**2*z(n)
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
            if (maxval(abs(beta_glnrho_scaled))/=0.0) then
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
              f(l,m,n,iux) = f(l,m,n,iux) - &
                  cs20*beta_glnrho_scaled(1)*eps*tausd(1)/ &
                  (1.0+2*eps+eps**2+(Omega*tausd(1))**2)
              f(l,m,n,iuy) = f(l,m,n,iuy) + &
                  cs20*beta_glnrho_scaled(1)* &
                  (1+eps+(Omega*tausd(1))**2)/ &
                  (2*Omega*(1.0+2*eps+eps**2+(Omega*tausd(1))**2))
              f(l,m,n,iudx(1)) = f(l,m,n,iudx(1)) + &
                  cs20*beta_glnrho_scaled(1)*tausd(1)/ &
                  (1.0+2*eps+eps**2+(Omega*tausd(1))**2)
              f(l,m,n,iudy(1)) = f(l,m,n,iudy(1)) + &
                  cs20*beta_glnrho_scaled(1)*(1+eps)/ &
                  (2*Omega*(1.0+2*eps+eps**2+(Omega*tausd(1))**2))
            elseif (beta_dPdr_dust_scaled/=0.0) then
              f(l,m,n,iux) = f(l,m,n,iux) - &
                  cs20*beta_dPdr_dust_scaled*eps*tausd(1)/ &
                  (1.0+2*eps+eps**2+(Omega*tausd(1))**2)
              f(l,m,n,iuy) = f(l,m,n,iuy) - &
                  cs20*beta_dPdr_dust_scaled*(eps+eps**2)/ &
                  (2*Omega*(1.0+2*eps+eps**2+(Omega*tausd(1))**2))
              f(l,m,n,iudx(1)) = f(l,m,n,iudx(1)) + &
                  cs20*beta_dPdr_dust_scaled*tausd(1)/ &
                  (1.0+2*eps+eps**2+(Omega*tausd(1))**2)
              f(l,m,n,iudy(1)) = f(l,m,n,iudy(1)) - &
                  cs20*beta_dPdr_dust_scaled* &
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
                u0_gas_pseudo*(1.0 + Omega_pseudo*tausd(1))/ &
                (1.0 + eps + Omega_pseudo*tausd(1))
            f(l,m,n,iudx) = f(l,m,n,iudx) + &
                u0_gas_pseudo/(1.0 + eps + Omega_pseudo*tausd(1))
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
              eps=sum(exp(f(l1:l2,m1:m2,n1:n2,ilnnd(1))))/ &
                  sum(f(l1:l2,m1:m2,n1:n2,irho))
            else
              eps=sum(f(l1:l2,m1:m2,n1:n2,ind(1)))/ &
                  sum(f(l1:l2,m1:m2,n1:n2,ilnrho))
            endif
          else
            if (ldustdensity_log) then
              eps=sum(exp(f(l1:l2,m1:m2,n1:n2,ilnnd(1))))/ &
                  sum(exp(f(l1:l2,m1:m2,n1:n2,ilnrho)))
            else
              eps=sum(f(l1:l2,m1:m2,n1:n2,ind(1)))/ &
                  sum(exp(f(l1:l2,m1:m2,n1:n2,ilnrho)))
            endif
          endif
!
          do m=m1,m2; do n=n1,n2
!
            f(l1:l2,m,n,ind(1)) = 0.0*f(l1:l2,m,n,ind(1)) + &
                eps*ampluud*cos(kz_uud*z(n))*cos(kx_uud*x(l1:l2))
!
            f(l1:l2,m,n,ilnrho) = f(l1:l2,m,n,ilnrho) + &
                ampluud* &
                ( real(coeff(7))*cos(kx_uud*x(l1:l2)) - &
                 aimag(coeff(7))*sin(kx_uud*x(l1:l2)))*cos(kz_uud*z(n))
!
            f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) + &
                eta_glnrho*v_Kepler*ampluud* &
                ( real(coeff(4))*cos(kx_uud*x(l1:l2)) - &
                 aimag(coeff(4))*sin(kx_uud*x(l1:l2)))*cos(kz_uud*z(n))
!
            f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + &
                eta_glnrho*v_Kepler*ampluud* &
                ( real(coeff(5))*cos(kx_uud*x(l1:l2)) - &
                 aimag(coeff(5))*sin(kx_uud*x(l1:l2)))*cos(kz_uud*z(n))
!
            f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) + &
                eta_glnrho*v_Kepler*(-ampluud)* &
                (aimag(coeff(6))*cos(kx_uud*x(l1:l2)) + &
                  real(coeff(6))*sin(kx_uud*x(l1:l2)))*sin(kz_uud*z(n))
!
            f(l1:l2,m,n,iudx(1)) = f(l1:l2,m,n,iudx(1)) + &
                eta_glnrho*v_Kepler*ampluud* &
                ( real(coeff(1))*cos(kx_uud*x(l1:l2)) - &
                 aimag(coeff(1))*sin(kx_uud*x(l1:l2)))*cos(kz_uud*z(n))
!
            f(l1:l2,m,n,iudy(1)) = f(l1:l2,m,n,iudy(1)) + &
                eta_glnrho*v_Kepler*ampluud* &
                ( real(coeff(2))*cos(kx_uud*x(l1:l2)) - &
                 aimag(coeff(2))*sin(kx_uud*x(l1:l2)))*cos(kz_uud*z(n))
!
            f(l1:l2,m,n,iudz(1)) = f(l1:l2,m,n,iudz(1)) + &
                eta_glnrho*v_Kepler*(-ampluud)* &
                (aimag(coeff(3))*cos(kx_uud*x(l1:l2)) + &
                  real(coeff(3))*sin(kx_uud*x(l1:l2)))*sin(kz_uud*z(n))
!
          enddo; enddo
!
!  Catch unknown values
!
        case default
          write (unit=errormsg,fmt=*) 'No such such value for inituu: ', trim(inituud(j))
          call fatal_error('init_uud',errormsg)

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
      if (lviscosity_dust) then
        if ((iviscd=='nud-const' .or. iviscd=='hyper3_nud-const') &
            .and. ldustdensity) then
          lpenc_requested(i_sdij)=.true.
          lpenc_requested(i_glnnd)=.true.
        endif
        if (iviscd=='simplified' .or. iviscd=='nud-const') &
            lpenc_requested(i_del2ud)=.true.
        if (iviscd=='hyper3_simplified' .or. iviscd=='hyper3_nud-const' .or. &
            iviscd=='hyper3_rhod_nud-const') &
            lpenc_requested(i_del6ud)=.true.
        if (iviscd=='nud-const' .or. iviscd=='hyper3_nud-const') &
            lpenc_requested(i_sdglnnd)=.true.
        if (iviscd=='nud-const') lpenc_requested(i_graddivud)=.true.
        if (iviscd=='hyper3_rhod_nud-const') lpenc_requested(i_rhod)=.true.
      endif
      if (beta_dPdr_dust/=0.) lpenc_requested(i_cs2)=.true.
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
        if (iviscd=='nud-const') then
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
!
      intent(in) :: f
      intent(inout) :: p
!
      do k=1,ndustspec
! uud
        if (lpencil(i_uud)) p%uud(:,:,k)=f(l1:l2,m,n,iudx(k):iudz(k))
! ud2
        if (lpencil(i_ud2)) call dot2_mn(p%uud(:,:,k),p%ud2(:,k))
! udij
        if (lpencil(i_udij)) call gij(f,iuud(k),p%udij(:,:,:,k),1)
! divud
        if (lpencil(i_divud)) &
            p%divud(:,k) = p%udij(:,1,1,k) + p%udij(:,2,2,k) + p%udij(:,3,3,k)
! udgud
        if (lpencil(i_udgud)) call multmv_mn(p%udij(:,:,:,k),p%uud(:,:,k),p%udgud(:,:,k))
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
          select case (iviscd)
          case ('nud-const')
            do j=1,3
              do i=1,3
                p%sdij(:,i,j,k)=.5*(p%udij(:,i,j,k)+p%udij(:,j,i,k))
              enddo
              p%sdij(:,j,j,k)=p%sdij(:,j,j,k)-.333333*p%divud(:,k)
            enddo
          case ('hyper3_nud-const')
            call gij(f,iuud(k),tmp_pencil_3x3,5)
            do i=1,3
              do j=1,3
                p%sdij(:,i,j,k)=tmp_pencil_3x3(:,i,j)
              enddo
            enddo
          case default
            if (headtt) then
              write (unit=errormsg,fmt=*) 'No rate-of-strain tensor matches iviscd=', iviscd
              call warning('calc_pencils_dustvelocity',errormsg)
            endif
          endselect
        endif
! del2ud
        if (lpencil(i_del2ud)) call del2v(f,iuud(k),p%del2ud(:,:,k))
! del6ud
        if (lpencil(i_del6ud)) call del6v(f,iuud(k),p%del6ud(:,:,k))
! graddivud
        if (lpencil(i_graddivud)) &
            call del2v_etc(f,iuud(k),GRADDIV=p%graddivud(:,:,k))
      enddo
!
    endsubroutine calc_pencils_dustvelocity
!***********************************************************************
    subroutine duud_dt(f,df,p)
!
!  Dust velocity evolution
!  Calculate duud/dt = - uud.graduud - 2Omega x uud - 1/tausd*(uud-uu)
!
!  18-mar-03/axel+anders: adapted from hydro
!
      use Diagnostics
      use EquationOfState, only: gamma
      use General
      use IO
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: fviscd,AA_sfta,BB_sfta
      real, dimension (nx) :: tausg1,mudrhod1
      real :: c2,s2 !(coefs for Coriolis force with inclined Omega)
      integer :: i,j,k
!
      intent(in) :: f,p
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
        call get_stoppingtime(f,p%rho,p%cs2,p%rhod(:,k),k)
!
!  Short stopping time approximation.
!  Calculated from master equation d(wx-ux)/dt = A + B*(wx-ux) = 0.
!
        if (ldustvelocity_shorttausd .and. &
            any(tausd1(:,k)>=shorttaus1limit)) then
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
              if (headtt .and. k == 1) &
                  print*,'duud_dt: add Coriolis force; Omega=',Omega
              c2=2*Omega
              df(l1:l2,m,n,iudx(k)) = df(l1:l2,m,n,iudx(k)) + c2*p%uud(:,2,k)
              df(l1:l2,m,n,iudy(k)) = df(l1:l2,m,n,iudy(k)) - c2*p%uud(:,1,k)
            else
              if (headtt .and. k == 1) print*, &
                  'duud_dt: Coriolis force; Omega,theta=',Omega,theta
              c2=2*Omega*cos(theta*pi/180.)
              s2=2*Omega*sin(theta*pi/180.)
              df(l1:l2,m,n,iudx(k)) = &
                  df(l1:l2,m,n,iudx(k)) + c2*p%uud(:,2,k)
              df(l1:l2,m,n,iudy(k)) = &
                  df(l1:l2,m,n,iudy(k)) - c2*p%uud(:,1,k) + s2*p%uud(:,3,k)
              df(l1:l2,m,n,iudz(k)) = &
                  df(l1:l2,m,n,iudz(k))                   + s2*p%uud(:,2,k)
            endif
          endif
!
!  Add drag force on dust
!
          if (ldragforce_dust) then
            do i=1,3
              df(l1:l2,m,n,iudx(k)-1+i) = df(l1:l2,m,n,iudx(k)-1+i) - &
                  tausd1(:,k)*(p%uud(:,i,k)-p%uu(:,i))
            enddo
!
!  Add drag force on gas (back-reaction from dust)
!
            if (ldragforce_gas) then
              tausg1 = p%rhod(:,k)*tausd1(:,k)*p%rho1
              if (tausgmin/=0.0) where (tausg1>=tausg1max) tausg1=tausg1max
              do i=1,3
                df(l1:l2,m,n,iux-1+i) = df(l1:l2,m,n,iux-1+i) - &
                   tausg1*(p%uu(:,i)-p%uud(:,i,k))
              enddo
              if (lfirst.and.ldt) dt1_max=max(dt1_max,(tausg1+tausd1(:,k))/cdtd)
            else
              if (lfirst.and.ldt) dt1_max=max(dt1_max,tausd1(:,k)/cdtd)
            endif
          endif
!
!  Add constant background pressure gradient beta=alpha*H0/r0, where alpha
!  comes from a global pressure gradient P = P0*(r/r0)^alpha.
!  (the term must be added to the dust equation of motion when measuring
!  velocities relative to the shear flow modified by the global pressure grad.)
!
          if (beta_dPdr_dust/=0.0) df(l1:l2,m,n,iudx(k)) = &
              df(l1:l2,m,n,iudx(k)) + p%cs2*beta_dPdr_dust_scaled
!
!  Add pseudo Coriolis force (to drive velocity difference between dust and gas)
!
          if (Omega_pseudo/=0.0) then
            df(l1:l2,m,n,iux) = &
                df(l1:l2,m,n,iux) - Omega_pseudo*(p%uu(:,1)-u0_gas_pseudo)
            df(l1:l2,m,n,iudx(:)) = &
                df(l1:l2,m,n,iudx(:)) - Omega_pseudo*p%uud(:,1,:)
          endif
!
!  Add viscosity on dust
!
          if (lviscosity_dust) then
!
          fviscd=0.0
! AJ: this only works if viscosity coefficient is same for all species:
          diffus_nud=0.0  ! Do not sum viscosity from all dust species
!
            select case (iviscd)
!
!  Viscous force: nud*del2ud
!     -- not physically correct (no momentum conservation)
!
            case ('simplified')
              if (headtt) print*, 'Viscous force (dust): nud*del2ud'
              fviscd = fviscd + nud(k)*p%del2ud(:,:,k)
              if (lfirst.and.ldt) diffus_nud=diffus_nud+nud(k)*dxyz_2
!
!  Viscous force: nud*(del2ud+graddivud/3+2Sd.glnnd)
!    -- the correct expression for nud=const
!
            case ('nud-const')
              if (headtt) print*, &
                  'Viscous force (dust): nud*(del2ud+graddivud/3+2Sd.glnnd)'
              if (ldustdensity) then
                fviscd = fviscd + 2*nud(k)*p%sdglnnd(:,:,k) + &
                    nud(k)*(p%del2ud(:,:,k)+1/3.*p%graddivud(:,:,k))
              else
                fviscd = fviscd + nud(k)*(p%del2ud(:,:,k)+1/3.*p%graddivud(:,:,k))
              endif
              if (lfirst.and.ldt) diffus_nud=diffus_nud+nud(k)*dxyz_2
!
!  Viscous force: nud*del6ud (not momentum-conserving)
!
            case ('hyper3_simplified')
              if (headtt) print*, 'Viscous force (dust): nud*del6ud'
              fviscd = fviscd + nud_hyper3(k)*p%del6ud(:,:,k)
              if (lfirst.and.ldt) diffus_nud3=diffus_nud3+nud_hyper3(k)*dxyz_6

            case ('hyper3_rhod_nud-const')
!
!  Viscous force: mud/rhod*del6ud
!
              if (headtt) print*, 'Viscous force (dust): mud/rhod*del6ud'
              mudrhod1=(nud_hyper3(k)*nd0*md0)/p%rhod(:,k)   ! = mud/rhod
              do i=1,3
                fviscd(:,i) = fviscd(:,i) + mudrhod1*p%del6ud(:,i,k)
              enddo
              if (lfirst.and.ldt) diffus_nud3=diffus_nud3+nud_hyper3(k)*dxyz_6

            case ('hyper3_nud-const')
!
!  Viscous force: nud*(del6ud+S.glnnd), where S_ij=d^5 ud_i/dx_j^5
!
              if (headtt) print*, 'Viscous force (dust): nud*(del6ud+S.glnnd)'
              fviscd = fviscd + nud_hyper3(k)*(p%del6ud(:,:,k)+p%sdglnnd(:,:,k))
              if (lfirst.and.ldt) diffus_nud3=diffus_nud3+nud_hyper3(k)*dxyz_6

            case default

              write (unit=errormsg,fmt=*) 'No such value for iviscd: ', trim(iviscd)
              call fatal_error('duud_dt',errormsg)

            endselect

          df(l1:l2,m,n,iudx(k):iudz(k)) = df(l1:l2,m,n,iudx(k):iudz(k)) + fviscd

          endif
!
!  ``uud/dx'' for timestep
!
          if (lfirst .and. ldt) then
            advec_uud=max(advec_uud,abs(p%uud(:,1,k))*dx_1(l1:l2)+ &
                                    abs(p%uud(:,2,k))*dy_1(  m  )+ &
                                    abs(p%uud(:,3,k))*dz_1(  n  ))
            if (idiag_dtud(k)/=0) &
                call max_mn_name(advec_uud/cdt,idiag_dtud(k),l_dt=.true.)
            if (idiag_dtnud(k)/=0) &
                 call max_mn_name(diffus_nud/cdtv,idiag_dtnud(k),l_dt=.true.)
          endif
          if ((headtt.or.ldebug) .and. (ip<6)) then
            print*,'duud_dt: max(advec_uud) =',maxval(advec_uud)
            print*,'duud_dt: max(diffus_nud) =',maxval(diffus_nud)
          endif
!
!  Short friction time switch.
!
        endif
!
!  End loop over dust species
!
      enddo
!
!  Calculate diagnostic variables
!
      if (ldiagnos) then
        do k=1,ndustspec
          if ((headtt.or.ldebug) .and. (ip<6)) &
              print*, 'duud_dt: Calculate diagnostic values...'
          if (idiag_udrms(k)/=0) &
              call sum_mn_name(p%ud2(:,k),idiag_udrms(k),lsqrt=.true.)
          if (idiag_udmax(k)/=0) &
              call max_mn_name(p%ud2(:,k),idiag_udmax(k),lsqrt=.true.)
          if (idiag_rdudmax(k)/=0) &
              call max_mn_name(p%rhod(:,k)**2*p%ud2(:,k),idiag_rdudmax(k), &
              lsqrt=.true.)
          if (idiag_ud2m(k)/=0) call sum_mn_name(p%ud2(:,k),idiag_ud2m(k))
          if (idiag_ekintot_dust/=0) &
             call integrate_mn_name(p%ud2,idiag_ekintot_dust)
          if (idiag_udxm(k)/=0) call sum_mn_name(p%uud(:,1,k),idiag_udxm(k))
          if (idiag_udym(k)/=0) call sum_mn_name(p%uud(:,2,k),idiag_udym(k))
          if (idiag_udzm(k)/=0) call sum_mn_name(p%uud(:,3,k),idiag_udzm(k))
          if (idiag_udx2m(k)/=0) &
              call sum_mn_name(p%uud(:,1,k)**2,idiag_udx2m(k))
          if (idiag_udy2m(k)/=0) &
              call sum_mn_name(p%uud(:,2,k)**2,idiag_udy2m(k))
          if (idiag_udz2m(k)/=0) &
              call sum_mn_name(p%uud(:,3,k)**2,idiag_udz2m(k))
          if (idiag_udm2(k)/=0) call max_mn_name(p%ud2(:,k),idiag_udm2(k))
          if (idiag_divud2m(k)/=0) &
              call sum_mn_name(p%divud(:,k)**2,idiag_divud2m(k))
          if (idiag_rdudxm(k)/=0) &
              call sum_mn_name(p%rhod(:,k)*p%uud(:,1,k),idiag_rdudxm(k))
          if (idiag_rdudym(k)/=0) &
              call sum_mn_name(p%rhod(:,k)*p%uud(:,2,k),idiag_rdudym(k))
          if (idiag_rdudzm(k)/=0) &
              call sum_mn_name(p%rhod(:,k)*p%uud(:,3,k),idiag_rdudzm(k))
          if (idiag_rdudx2m(k)/=0) &
              call sum_mn_name((p%rhod(:,k)*p%uud(:,1,k))**2,idiag_rdudx2m(k))
          if (idiag_odrms(k)/=0) &
              call sum_mn_name(p%od2,idiag_odrms(k),lsqrt=.true.)
          if (idiag_odmax(k)/=0) &
              call max_mn_name(p%od2,idiag_odmax(k),lsqrt=.true.)
          if (idiag_od2m(k)/=0) call sum_mn_name(p%od2,idiag_od2m(k))
          if (idiag_oudm(k)/=0) call sum_mn_name(p%oud,idiag_oudm(k))
        enddo
      endif
!
!  xy-averages
!
      if (l1davgfirst) then
        do k=1,ndustspec
          if (idiag_udxmz(k)/=0) &
              call xysum_mn_name_z(p%uud(:,1,k),idiag_udxmz(k))
          if (idiag_udymz(k)/=0) &
              call xysum_mn_name_z(p%uud(:,2,k),idiag_udymz(k))
          if (idiag_udzmz(k)/=0) &
              call xysum_mn_name_z(p%uud(:,3,k),idiag_udzmz(k))
          if (idiag_udx2mz(k)/=0) &
              call xysum_mn_name_z(p%uud(:,1,k)**2,idiag_udx2mz(k))
          if (idiag_udy2mz(k)/=0) &
              call xysum_mn_name_z(p%uud(:,2,k)**2,idiag_udy2mz(k))
          if (idiag_udz2mz(k)/=0) &
              call xysum_mn_name_z(p%uud(:,3,k)**2,idiag_udz2mz(k))
        enddo
      endif
!
!  2-D averages.
!
      if (l2davgfirst) then
        do k=1,ndustspec
          if (idiag_udxmxy(k)/=0) &
              call zsum_mn_name_xy(p%uud(:,1,k),idiag_udxmxy(k))
          if (idiag_udymxy(k)/=0) &
              call zsum_mn_name_xy(p%uud(:,2,k),idiag_udymxy(k))
          if (idiag_udzmxy(k)/=0) &
              call zsum_mn_name_xy(p%uud(:,3,k),idiag_udzmxy(k))
        enddo
      endif
!
    endsubroutine duud_dt
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
!  Calculate surface of dust particles
!
      integer :: i,j
!
      do i=1,ndustspec
        do j=1,ndustspec
          scolld(i,j) = pi*(ad(i)+ad(j))**2
        enddo
      enddo
!
    endsubroutine get_dustcrosssection
!***********************************************************************
    subroutine get_stoppingtime(f,rho,cs2,rhod,k)
!
!  Calculate stopping time depending on choice of drag law.
!
      use Sub, only: dot2

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: rho,rhod,csrho,cs2,deltaud2
      integer :: k
!
      select case (draglaw)

      case ('epstein_cst')
        ! Do nothing, initialized in initialize_dustvelocity
      case ('epstein_cst_b')
        tausd1(:,k) = betad(k)/rhod
      case ('stokes_cst_tausd')
        tausd1(:,k) = betad(k)
      case ('epstein_var')
        call dot2(f(l1:l2,m,n,iudx(k):iudz(k))-f(l1:l2,m,n,iux:iuz),deltaud2)
        csrho       = sqrt(cs2+deltaud2)*rho
        tausd1(:,k) = csrho*rhodsad1(k)
      case ('epstein_gaussian_z')
        tausd1(:,k) = (1/tausd(k))*exp(-z(n)**2/(2*scaleHtaus**2))
        if (z0taus/=0.0) then
          tausd1(:,k)=tausd1(:,k)/( &
            0.5*(tanh((z(n)+z0taus)/widthtaus)+tanh((-z(n)+z0taus)/widthtaus)))
        endif
      case default
        call fatal_error("get_stoppingtime","No valid drag law specified.")

      endselect
!
    endsubroutine get_stoppingtime
!***********************************************************************
    subroutine read_dustvelocity_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=dustvelocity_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=dustvelocity_init_pars,ERR=99)
      endif


99    return
    endsubroutine read_dustvelocity_init_pars
!***********************************************************************
    subroutine write_dustvelocity_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=dustvelocity_init_pars)

    endsubroutine write_dustvelocity_init_pars
!***********************************************************************
    subroutine read_dustvelocity_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=dustvelocity_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=dustvelocity_run_pars,ERR=99)
      endif


99    return
    endsubroutine read_dustvelocity_run_pars
!***********************************************************************
    subroutine write_dustvelocity_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=dustvelocity_run_pars)

    endsubroutine write_dustvelocity_run_pars
!***********************************************************************
    subroutine rprint_dustvelocity(lreset,lwrite)
!
!  reads and registers print parameters relevant for hydro part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Diagnostics
      use General, only: chn
!
      integer :: iname,inamez,k
      logical :: lreset,lwr
      logical, optional :: lwrite
      character (len=5) :: sdust
!
!  Write information to index.pro that should not be repeated for i
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite

      if (lwr) then
        write(3,*) 'ndustspec=',ndustspec
        write(3,*) 'nname=',nname
      endif
!
!  reset everything in case of reset
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
          call chn(k,sdust)
          if (ndustspec == 1) sdust=''
          call parse_name(iname,cname(iname),cform(iname), &
              'dtud'//trim(sdust),idiag_dtud(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'dtnud'//trim(sdust),idiag_dtnud(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udxm'//trim(sdust),idiag_udxm(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udym'//trim(sdust),idiag_udym(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udzm'//trim(sdust),idiag_udzm(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'ud2m'//trim(sdust),idiag_ud2m(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'ekintot_dust'//trim(sdust),idiag_ekintot_dust)
          call parse_name(iname,cname(iname),cform(iname), &
              'udx2m'//trim(sdust),idiag_udx2m(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udy2m'//trim(sdust),idiag_udy2m(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udz2m'//trim(sdust),idiag_udz2m(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udm2'//trim(sdust),idiag_udm2(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'od2m'//trim(sdust),idiag_od2m(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'oudm'//trim(sdust),idiag_oudm(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udrms'//trim(sdust),idiag_udrms(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udmax'//trim(sdust),idiag_udmax(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'rdudmax'//trim(sdust),idiag_rdudmax(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'rdudxm'//trim(sdust),idiag_rdudxm(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'rdudym'//trim(sdust),idiag_rdudym(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'rdudzm'//trim(sdust),idiag_rdudzm(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'rdudx2m'//trim(sdust),idiag_rdudx2m(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'odrms'//trim(sdust),idiag_odrms(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'odmax'//trim(sdust),idiag_odmax(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udmx'//trim(sdust),idiag_udmx(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udmy'//trim(sdust),idiag_udmy(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udmz'//trim(sdust),idiag_udmz(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'divud2m'//trim(sdust),idiag_divud2m(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'epsKd'//trim(sdust),idiag_epsKd(k))
        enddo
!
!  check for those quantities for which we want xy-averages
!
        do inamez=1,nnamez
          call parse_name(inamez,cnamez(inamez),cformz(inamez), &
              'udxmz'//trim(sdust),idiag_udxmz(k))
          call parse_name(inamez,cnamez(inamez),cformz(inamez), &
              'udymz'//trim(sdust),idiag_udymz(k))
          call parse_name(inamez,cnamez(inamez),cformz(inamez), &
              'udzmz'//trim(sdust),idiag_udzmz(k))
          call parse_name(inamez,cnamez(inamez),cformz(inamez), &
              'udx2mz'//trim(sdust),idiag_udx2mz(k))
          call parse_name(inamez,cnamez(inamez),cformz(inamez), &
              'udy2mz'//trim(sdust),idiag_udy2mz(k))
          call parse_name(inamez,cnamez(inamez),cformz(inamez), &
              'udz2mz'//trim(sdust),idiag_udz2mz(k))
        enddo
!
!  End loop over dust layers
!
      enddo
!
    endsubroutine rprint_dustvelocity
!***********************************************************************
    subroutine get_slices_dustvelocity(f,slices)
!
!  Write slices for animation of Dustvelocity variables.
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
!  Dustvelocity.
!
        case ('uud')
          if (slices%index>=3) then
            slices%ready=.false.
          else
            slices%index=slices%index+1
            slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,iudx(1)-1+slices%index)
            slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,iudx(1)-1+slices%index)
            slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,iudx(1)-1+slices%index)
            slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,iudx(1)-1+slices%index)
            if (lwrite_slice_xy3) &
                slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,iudx(1)-1+slices%index)
            if (lwrite_slice_xy4) &
                slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,iudx(1)-1+slices%index)
            if (slices%index<=3) slices%ready=.true.
          endif
      endselect
!
    endsubroutine get_slices_dustvelocity
!***********************************************************************
endmodule Dustvelocity
