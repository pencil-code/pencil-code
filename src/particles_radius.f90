! $Id$
!
!  This module takes care of everything related to particle radius.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 1
! CPARAM logical, parameter :: lparticles_radius=.true.
!
!***************************************************************
module Particles_radius
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_sub
  use Particles_chemistry
!
  implicit none
!
  include 'particles_radius.h'
!
  real :: vthresh_sweepup=-1.0, deltavp12_floor=0.0
  real, dimension(ndustrad) :: ap0=0.0
  real, dimension(ndustrad) :: radii_distribution=0.0
  real :: tstart_sweepup_par=0.0, cdtps=0.2, cdtpc=0.2
  real :: tstart_condensation_par=0.0
  real :: apmin=0.0, latent_heat_SI=2.257e6, alpha_cond=1.0, alpha_cond1=1.0
  real :: diffusion_coefficient=1.0, diffusion_coefficient1=1.0
  real :: tau_damp_evap=0.0, tau_damp_evap1=0.0
  real :: tau_ocean_driving=0.0, tau_ocean_driving1=0.0
  real :: ztop_ocean=0.0, TTocean=300.0
  real :: aplow=1.0, apmid=1.5, aphigh=2.0, mbar=1.0
  real :: ap1=1.0, qplaw=0.0,vapor_mixing_ratio_qvs=0., GS_condensation=0.0
  real :: modified_vapor_diffusivity=6.74e-6, n0_mean = 1.e8
  real, pointer :: G_condensation, ssat0 
  real :: sigma_initdist=0.2, a0_initdist=5e-6, rpbeta0=0.0
  real :: xi_accretion = 0.
  integer :: nbin_initdist=20, ip1=npar/2
  logical :: lsweepup_par=.false., lcondensation_par=.false.
  logical :: llatent_heat=.true., lborder_driving_ocean=.false.
  logical :: lcondensation_simplified=.false.
  logical :: lcondensation_rate=.false., ldt_evaporation=.false.
  logical :: ldt_condensation=.false., ldt_condensation_off=.false.
  logical :: lconstant_radius_w_chem=.false.
  logical :: lfixed_particles_radius=.false.
  logical :: reinitialize_ap=.false.
  logical :: ltauascalar = .false., ldust_condensation=.false.
  logical :: ldust_accretion = .false.
  character(len=labellen), dimension(ninit) :: initap='nothing'
  character(len=labellen) :: condensation_coefficient_type='constant'
!
  namelist /particles_radius_init_pars/ &
      initap, ap0, rhopmat, vthresh_sweepup, deltavp12_floor, &
      lsweepup_par, lcondensation_par, tstart_sweepup_par, cdtps, apmin, &
      condensation_coefficient_type, alpha_cond, diffusion_coefficient, &
      tau_damp_evap, llatent_heat, cdtpc, tau_ocean_driving, &
      lborder_driving_ocean, ztop_ocean, radii_distribution, TTocean, &
      aplow, apmid, aphigh, mbar, ap1, ip1, qplaw, eps_dtog, nbin_initdist, &
      sigma_initdist, a0_initdist, rpbeta0, lparticles_radius_rpbeta, &
      lfixed_particles_radius
!
  namelist /particles_radius_run_pars/ &
      rhopmat, vthresh_sweepup, deltavp12_floor, &
      lsweepup_par, lcondensation_par, tstart_sweepup_par, cdtps, apmin, &
      condensation_coefficient_type, alpha_cond, diffusion_coefficient, &
      tau_damp_evap, llatent_heat, cdtpc, tau_ocean_driving, &
      lborder_driving_ocean, ztop_ocean, TTocean, &
      lcondensation_simplified, GS_condensation, rpbeta0, &
      lfixed_particles_radius, &
      lconstant_radius_w_chem, &
      reinitialize_ap, initap, &
      lcondensation_rate, vapor_mixing_ratio_qvs, &
      ltauascalar, modified_vapor_diffusivity, ldt_evaporation, &
      ldt_condensation, ldt_condensation_off, &
      ldust_condensation, xi_accretion, ldust_accretion, &
      tstart_condensation_par
!
  integer :: idiag_apm=0, idiag_ap2m=0, idiag_apmin=0, idiag_apmax=0
  integer :: idiag_dvp12m=0, idiag_dtsweepp=0, idiag_npswarmm=0
  integer :: idiag_ieffp=0
!
  contains
!***********************************************************************
    subroutine register_particles_radius()
!
!  Set up indices for access to the fp and dfp arrays.
!
!  22-aug-05/anders: coded
!
      if (lroot) call svn_id( "$Id$")
!
!  Index for particle radius.
!
      call append_npvar('iap',iap)
!
! Index for the effectiveness factor of surface reactions
!
      if (lparticles_chemistry) call append_npaux('ieffp',ieffp)
!
      if (lparticles_radius_rpbeta) call append_npvar('irpbeta',irpbeta)
!
    endsubroutine register_particles_radius
!***********************************************************************
    subroutine initialize_particles_radius(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  22-aug-05/anders: coded
!
      use SharedVariables, only: put_shared_variable, get_shared_variable
!
      real, dimension(mx,my,mz,mfarray) :: f
!
!  Calculate the number density of bodies within a superparticle.
!
      if (npart_radii > 1 .and. &
          (.not. lcartesian_coords .or. lparticles_number .or. lparticles_spin)) then
        call fatal_error('initialize_particles_radius: npart_radii > 1','')
      else
        mpmat = 4/3.0*pi*rhopmat*ap0(1)**3
        if (lroot) print*, 'initialize_particles_radius: '// &
            'mass per dust grain mpmat=', mpmat
      endif
!
      if ((lsweepup_par .or. lcondensation_par).and. .not. lpscalar &
          .and. .not. lascalar &
          .and. .not. lcondensation_simplified) then
        call fatal_error('initialize_particles_radius', &
            'must have passive scalar module for sweep-up and condensation')
      endif
!
!  Short hand for spherical particle prefactor.
!
      four_pi_rhopmat_over_three = four_pi_over_three*rhopmat
!
!  Inverse coefficients.
!
      alpha_cond1 = 1/alpha_cond
      diffusion_coefficient1 = 1/diffusion_coefficient
      if (tau_damp_evap /= 0.0) tau_damp_evap1 = 1/tau_damp_evap
      if (tau_ocean_driving /= 0.0) tau_ocean_driving1 = 1/tau_ocean_driving
!
      call put_shared_variable('ap0',ap0,caller='initialize_particles_radius')
!
! If we have decided to hold the radius of the particles to be fixed then
! we should not have any process that changes the radius.
!
      if (lfixed_particles_radius) then
        if (lsweepup_par) &
            call fatal_error('initialize_particles_radius', &
            'incosistency: lfixed_particles_radius and lsweepup_par cannot both be true')
        if (lcondensation_par) &
            call fatal_error('initialize_particles_radius', &
            'incosistency: lfixed_particles_radius and lcondensation_par cannot both be true')
        if (lparticles_chemistry) &
            call fatal_error('initialize_particles_radius', &
            'incosistency: lfixed_particles_radius and lparticles_chemistry cannot both be true')
      endif
!
      call keep_compiler_quiet(f)
!
      if (lascalar) then
        call get_shared_variable('G_condensation', G_condensation)
        if (lcondensation_rate) call get_shared_variable('ssat0', ssat0)
      endif
!
    endsubroutine initialize_particles_radius
!***********************************************************************
    subroutine set_particle_radius(f,fp,npar_low,npar_high,init)
!
!  Set radius of new particles.
!
!  18-sep-09/nils: adapted from init_particles_radius
!
      use General, only: random_number_wrapper
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mparray) :: fp
      integer :: npar_low, npar_high
      logical, optional :: init
!
      real, dimension(mpar_loc) :: r_mpar_loc, p_mpar_loc, tmp_mpar_loc
      real, dimension(ndustrad) :: radii_cumulative
      real, dimension(nbin_initdist) :: n_initdist, a_initdist
      integer, dimension(nbin_initdist) :: nn_initdist
      real :: radius_fraction, mcen, mmin, mmax, fcen, p
      real :: lna0, lna1, lna, lna0_initdist
      integer :: i, j, k, kend, ind, ibin
      logical :: initial
!
      initial = .false.
      if (present(init)) then
        if (init) initial = .true.
      endif
!
      do j = 1,ninit
!
        select case (initap(j))
!
        case ('nothing')
          if (initial .and. lroot .and. j == 1)  print*, 'set_particles_radius: nothing'
!
        case ('constant')
          if (initial .and. lroot) print*, 'set_particles_radius: constant radius'
          ind = 1
          do k = npar_low,npar_high
            if (npart_radii > 1) then
              call random_number_wrapper(radius_fraction)
              ind = ceiling(npart_radii*radius_fraction)
            endif
            fp(k,iap) = ap0(ind)
          enddo
!
        case ('constant-1')
          if (initial .and. lroot) print*, 'set_particles_radius: set particle 1 radius'
          do k = npar_low,npar_high
            if (ipar(k) == 1) fp(k,iap) = ap1
          enddo
!
        case ('2-size')
          if (initial .and. lroot) print*, 'set_particles_radius: give particles two radii'
          do k = npar_low,npar_high
            if (ipar(k)<=ip1) then
              fp(k,iap)=aplow
            else
              fp(k,iap)=aphigh
            endif
          enddo
!
        case ('2-size-alternate')
          if (initial .and. lroot) print*, 'set_particles_radius: give particles alternating two radii'
          do k = npar_low,npar_high
            if (mod(k,2)==0) then
              fp(k,iap)=aplow
            else
              fp(k,iap)=aphigh
            endif
          enddo
!
        case ('3-size-alternate')
          if (initial .and. lroot) print*, 'set_particles_radius: give particles alternating three radii'
          do k = npar_low,npar_high
            if (mod(k,3)==0) then
              fp(k,iap)=aplow
            else if (mod(k,3)==1) then
              fp(k,iap)=apmid
            else
              fp(k,iap)=aphigh
            endif
          enddo
!
        case ('random')
          if (initial .and. lroot) print*, 'set_particles_radius: random radius'
          do k = npar_low,npar_high
            call random_number_wrapper(fp(k,iap))
            fp(k,iap) = fp(k,iap)*(aphigh-aplow)+aplow
          enddo
!
        case ('logarithmically-spaced')
          if (initial .and. lroot) print*, 'set_particles_radius: '// &
              'logarithmically spaced with ap0, ap1=', ap0(j), ap1
          do k = npar_low,npar_high
            fp(k,iap) = 10**(alog10(ap0(j))+(ipar(k)-0.5)* &
                (alog10(ap1)-alog10(ap0(j)))/npar)
          enddo
!
          lspd: if (lparticles_density) then
            lsqp: if (qplaw /= 4) then
              lsk: do k = npar_low, npar_high
                aplow  = 10**(alog10(fp(k,iap)) - 0.5 * (alog10(ap1) - alog10(ap0(j))) / npar)
                aphigh = 10**(alog10(fp(k,iap)) + 0.5 * (alog10(ap1) - alog10(ap0(j))) / npar)
                fp(k,irhopswarm) = (aphigh**(4-qplaw) - aplow**(4-qplaw)) / &
                                   (ap1**(4-qplaw) - ap0(j)**(4-qplaw)) * rhop_swarm * real(npar)
              enddo lsk
            else lsqp
              fp(npar_low:npar_high,irhopswarm) = rhop_swarm
            endif lsqp
          endif lspd
!
!  Lognormal distribution. Here, ap1 is the largest value in the distribution.
!  Initialize particle radii by a direct probabilistic calculation using
!  gaussian noise for ln(a/a0)/sigma.
!
        case ('lognormal')
!
          if (initial .and. lroot) print*, 'set_particles_radius: '// &
              'lognormal=', a0_initdist
          call random_number_wrapper(r_mpar_loc)
          call random_number_wrapper(p_mpar_loc)
          tmp_mpar_loc = sqrt(-2*log(r_mpar_loc))*sin(2*pi*p_mpar_loc)
          fp(:,iap) = a0_initdist*exp(sigma_initdist*tmp_mpar_loc)
!
!  Lognormal distribution. Here, ap1 is the largest value in the distribution
!  and ap0 is the smallest radius initially.
!
        case ('old_lognormal')
!
          if (initial .and. lroot) print*, 'set_particles_radius: '// &
              'lognormal=', ap0(1), ap1
          lna0 = log(ap0(1))
          lna1 = log(ap1)
          lna0_initdist = log(a0_initdist)
          do ibin = 1,nbin_initdist
            lna = lna0+(lna1-lna0)*(ibin-1)/(nbin_initdist-1)
            a_initdist(ibin) = exp(lna)
            n_initdist(ibin) = exp(-0.5*(lna-lna0_initdist)**2/sigma_initdist**2) &
                /(sqrt(twopi)*sigma_initdist*a_initdist(ibin))
          enddo
          nn_initdist = nint(n_initdist)
          nn_initdist = (npar_high-npar_low+1)*nn_initdist/(sum(nn_initdist)-1)
!
!  is now normalized to the number of particles,
!  so set the corresponding distribution.
!  Normally, the number of bins is less than the number of available ones,
!  so then we patch the rest with fp(kend+1:,iap)=a0_initdist.
!
          k = npar_low
          do ibin = 1,nbin_initdist
            kend = k+nn_initdist(ibin)
            if (kend > k) then
              fp(k:kend,iap) = a_initdist(ibin)
              k = kend
            endif
          enddo
!
!  put all the remaining particles (from kend+1 to the end of the array)
!  in the bin corresponding to the middle of the distribution.
!
          fp(kend+1:,iap) = a0_initdist
!
        case ('specify')
!
!  User specified particle size distribution with constant radii.
!
          if (initial .and. lroot) &
              print*, 'set_particles_radius: constant radius, user specified distribution'
          radii_cumulative = 0.0
          radii_cumulative(1) = radii_distribution(1)
          do i = 2,npart_radii
            radii_cumulative(i) = radii_cumulative(i-1) + radii_distribution(i)
          enddo
          if (radii_cumulative(npart_radii) /= 1.0) then
!
!  Renormalize.
!
            do i = 1,npart_radii
              radii_cumulative(i) = radii_cumulative(i)/radii_cumulative(npart_radii)
            enddo
          endif
!
          do k = npar_low,npar_high
            call random_number_wrapper(radius_fraction)
            do i = 1,npart_radii
              if (radius_fraction <= radii_cumulative(i)) then
                fp(k,iap) = ap0(i)
                exit
              endif
            enddo
          enddo
!
!  Coagulation test with linear kernel. We initially put particles according
!  to the distribution function
!
!    fk = dn/dm = n_0/mbar_0*exp(-m/mbar_0)
!        => rhok = m_k*n_0/mbar_0*exp(-m/mbar_0)
!
!  Integrating rhok over all mass and normalising by n_0*mbar_0 gives
!
!    I = int_0^m[rhok]/int_0^oo[rhok] = 1 - exp(-x) - x*exp(-x)
!
!  where x=mk/mbar_0. We place particles equidistantly along the integral
!  to obtain the right distribution function.
!
        case ('kernel-lin')
          if (initial .and. lroot) print*, 'set_particles_radius: '// &
              'initial condition for linear kernel test'
          do k = npar_low,npar_high
            p = (ipar(k)-0.5)/float(npar)
            mmin = 0.0
            mmax = 1.0
            do while (.true.)
              if ((1.0-exp(-mmax)*(1+mmax)) < p) then
                mmax = mmax+1.0
              else
                exit
              endif
            enddo
!
            mcen = 0.5*(mmin+mmax)
            fcen = 1.0-exp(-mcen)*(1+mcen)
!
            do while (abs(p-fcen) > 1.0e-6)
              if (fcen < p) then
                mmin = mcen
              else
                mmax = mcen
              endif
              mcen = 0.5*(mmin+mmax)
              fcen = 1.0-exp(-mcen)*(1+mcen)
            enddo
!
            fp(k,iap) = (mcen*mbar/four_pi_rhopmat_over_three)**(1.0/3.0)
!
          enddo
!
        case ('power-law')
          call random_number_wrapper(fp(npar_low:npar_high,iap))
          fp(npar_low:npar_high,iap) = ((aphigh**(qplaw+1.0)-aplow**(qplaw+1.0)) &
              *fp(npar_low:npar_high,iap)+aplow**(qplaw+1.0))**(1.0/(qplaw+1.0))
!
        case default
          if (lroot) print*, 'init_particles_radius: '// &
              'No such such value for initap: ', trim(initap(j))
          call fatal_error('init_particles_radius','')
        endselect
      enddo
!
!  Set initial particle radius if lparticles_mass=T
!
      if (lparticles_mass) fp(:,iapinit) = fp(:,iap)
!
! Reinitialize particle radius if reinitialize_ap=T
! 09-Feb-17/Xiangyu: adapted from set_particles_radius
      if (reinitialize_ap) then
        do j = 1,ninit
          select case (initap(j))
!
          case ('constant')
            if (initial .and. lroot) print*, 'set_particles_radius: constant radius'
            ind = 1
            do k = npar_low,npar_high
              if (npart_radii > 1) then
                call random_number_wrapper(radius_fraction)
                ind = ceiling(npart_radii*radius_fraction)
              endif
              fp(k,iap) = ap0(ind)
            enddo
!
          case ('constant-1')
            if (initial .and. lroot) print*, 'set_particles_radius: set particle 1 radius'
            do k = npar_low,npar_high
              if (ipar(k) == 1) fp(k,iap) = ap1
            enddo
!
          case ('lognormal')
!
            if (initial .and. lroot) print*, 'set_particles_radius: '// &
                'lognormal=', a0_initdist
            call random_number_wrapper(r_mpar_loc)
            call random_number_wrapper(p_mpar_loc)
            tmp_mpar_loc = sqrt(-2*log(r_mpar_loc))*sin(2*pi*p_mpar_loc)
            fp(:,iap) = a0_initdist*exp(sigma_initdist*tmp_mpar_loc)
          endselect
        enddo
      endif
!
      if (lparticles_radius_rpbeta) &
        fp(npar_low:npar_high,irpbeta) = rpbeta0/(fp(npar_low:npar_high,iap)*rhopmat)
!
      call keep_compiler_quiet(f)
!
    endsubroutine set_particle_radius
!***********************************************************************
    subroutine pencil_criteria_par_radius()
!
!  All pencils that the Particles_radius module depends on are specified here.
!
!  21-nov-06/anders: coded
!
      if (lsweepup_par) then
        lpenc_requested(i_uu) = .true.
        lpenc_requested(i_rho) = .true.
        lpenc_requested(i_cc) = .true.
      endif
      if (lcondensation_par) then
        lpenc_requested(i_csvap2) = .true.
        lpenc_requested(i_TT1) = .true.
        lpenc_requested(i_ppvap) = .true.
        lpenc_requested(i_rho) = .true.
        lpenc_requested(i_rho1) = .true.
        lpenc_requested(i_cc) = .true.
        lpenc_requested(i_cc1) = .true.
        lpenc_requested(i_np) = .true.
        lpenc_requested(i_rhop) = .true.
        if (ltemperature) then
          lpenc_requested(i_cv1) = .true.
          lpenc_requested(i_TT1) = .true.
          lpenc_requested(i_TT) = .true.
        endif
        if (lascalar) then
          lpenc_requested(i_acc) = .true.
          lpenc_requested(i_ssat) = .true.
        endif
      endif
!
    endsubroutine pencil_criteria_par_radius
!***********************************************************************
    subroutine dap_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle radius.
!
!  22-aug-05/anders: coded
!
      use Particles_number
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mpar_loc,mpvar) :: dfp
      type (pencil_case) :: p
      integer, dimension(mpar_loc,3) :: ineargrid
      logical :: lfirstcall=.true., lheader
      integer :: k, k1, k2
      real :: mass_per_radius, rho
      real, dimension(:), allocatable :: effectiveness_factor, mass_loss
!
      intent(in) :: f
      intent(out) :: dfp
      intent(inout) :: fp
!
!  Print out header information in first time step.
!
      lheader = lfirstcall .and. lroot .and.(.not. lpencil_check_at_work)
!
!  Identify module.
!
      if (lheader) print*,'dap_dt_pencil: Calculate dap/dt'
!
      if (lsweepup_par) call dap_dt_sweepup_pencil(f,df,fp,dfp,p,ineargrid)
      if (lcondensation_par) &
          call dap_dt_condensation_pencil(f,df,fp,dfp,p,ineargrid)
!
!
      lfirstcall = .false.
!
!  Decrease particle radius with dmass if the density in the outer shell
!  is wholly consumed (equation 51 in Equations to be solved for DNS
!  reactive particles in a turbulent flow)
!
      if (lparticles_chemistry .and. npar_imn(imn) /= 0) then
        k1 = k1_imn(imn)
        k2 = k2_imn(imn)
!
!  Particles can lose radius under consumption
!
        if (.not. lconstant_radius_w_chem) then
!
!  mass loss has a positive value -> particle is losing mass
!  (the mass vector is pointing out of the particle)
!
          allocate(mass_loss(k1:k2))
          allocate(effectiveness_factor(k1:k2))
!
          call get_radius_chemistry(mass_loss,effectiveness_factor)
!
          if (.not. lsurface_nopores) then
            do k = k1,k2
              if (fp(k,irhosurf) < 0) then
                rho = fp(k,imp) / (fp(k,iap)**3 * 4./3. * pi )
                mass_per_radius = 4. * pi * rho * fp(k,iap)**2
                dfp(k,iap) = dfp(k,iap) - mass_loss(k) *(1-effectiveness_factor(k))/mass_per_radius
              endif
              fp(k,ieffp) = effectiveness_factor(k)
            enddo
          else
            do k = k1,k2
              rho = fp(k,imp) / (fp(k,iap)**3 * 4./3. * pi )
              mass_per_radius = 4. * pi * rho * fp(k,iap)**2
              dfp(k,iap) = dfp(k,iap) - mass_loss(k)/mass_per_radius
            enddo
          endif
          deallocate(mass_loss)
          deallocate(effectiveness_factor)
        else
!
!  Constant particle radius with activated chemistry
!
          dfp(k1:k2,iap) = 0.0
!
        endif
      endif
!
    endsubroutine dap_dt_pencil
!***********************************************************************
    subroutine dap_dt_sweepup_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Increase in particle radius due to sweep-up of small grains in the gas.
!
!  22-aug-05/anders: coded
!
      use Diagnostics, only: max_mn_name
      use Particles_number
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mpar_loc,mpvar) :: dfp
      type (pencil_case) :: p
      integer, dimension(mpar_loc,3) :: ineargrid
!
      real, dimension(nx) :: dt1_sweepup
      real :: deltavp
      integer :: k, ix0, ix
!
      intent(in) :: f, fp
      intent(inout) :: dfp
!
!  Increase in particle radius due to sweep-up of small grains in the gas.
!
      if (t >= tstart_sweepup_par) then
!
!
        if (npar_imn(imn) /= 0) then
          do k = k1_imn(imn),k2_imn(imn)
            ix0 = ineargrid(k,1)
            ix = ix0-nghost
!  No interpolation needed here.
!  Relative speed.
            deltavp = sqrt( &
                (fp(k,ivpx)-p%uu(ix,1))**2 + &
                (fp(k,ivpy)-p%uu(ix,2))**2 + &
                (fp(k,ivpz)-p%uu(ix,3))**2 )
            if (deltavp12_floor /= 0.0) &
                deltavp = sqrt(deltavp**2+deltavp12_floor**2)
!  Allow boulders to sweep up small grains if relative velocity not too high.
            if (deltavp <= vthresh_sweepup .or. vthresh_sweepup < 0.0) then
!  Radius increase due to sweep-up.
              dfp(k,iap) = dfp(k,iap) + 0.25*deltavp*p%cc(ix)*p%rho(ix)*rhopmat1
!
!  Deplete gas of small grains.
!
              if (lparticles_number) np_swarm = fp(k,inpswarm)
              if (lpscalar_nolog) then
                df(ix0,m,n,icc) = df(ix0,m,n,icc) - &
                    np_swarm*pi*fp(k,iap)**2*deltavp*p%cc(ix)
              else
                df(ix0,m,n,ilncc) = df(ix0,m,n,ilncc) - &
                    np_swarm*pi*fp(k,iap)**2*deltavp
              endif
!
!  Time-step contribution of sweep-up.
!
              if (lfirst .and. ldt) then
                dt1_sweepup(ix) = dt1_sweepup(ix) + &
                    np_swarm*pi*fp(k,iap)**2*deltavp
              endif
!
            endif
!
            if (ldiagnos) then
              if (idiag_dvp12m /= 0) call sum_par_name((/deltavp/),idiag_dvp12m)
            endif
          enddo
        endif
!
!  Time-step contribution of sweep-up.
!
        if (lfirst .and. ldt) then
          dt1_sweepup = dt1_sweepup/cdtps
          dt1_max = max(dt1_max,dt1_sweepup)
          if (ldiagnos .and. idiag_dtsweepp /= 0) &
              call max_mn_name(dt1_sweepup,idiag_dtsweepp,l_dt=.true.)
        endif
!
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine dap_dt_sweepup_pencil
!***********************************************************************
    subroutine dap_dt_condensation_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Change in particle radius due to condensation of monomers from the gas and
!  evaporation of solid/liquid droplets.
!
!  15-jan-10/anders: coded
!
      use EquationOfState, only: gamma, rho0
      use Particles_number
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mpar_loc,mpvar) :: dfp
      type (pencil_case) :: p
      integer, dimension(mpar_loc,3) :: ineargrid
!
      real, dimension(nx) :: ap_equi, vth, dt1_condensation, rhovap
      real, dimension(nx) :: total_surface_area, ppsat
      real, dimension(nx) :: rhocond_tot, rhosat, np_total
      real, dimension(nx) :: tau_phase1, tau_evaporation1, tau_evaporation
      real :: dapdt, drhocdt, alpha_cond_par
      integer :: k, ix, ix0
!
      intent(in) :: f, fp
      intent(inout) :: dfp
!
!  Change in particle radius due to condensation and evaporation.
!
      if (t >= tstart_condensation_par) then
!
        if (npar_imn(imn) /= 0) then
!          rhovap=p%cc(:,1)*p%rho
!DMDM
          rhovap = p%cc*p%rho
          ppsat = 6.035e11*exp(-5938*p%TT1)  ! Valid for water
          vth = sqrt(p%csvap2)
          rhosat = gamma*ppsat/p%csvap2
          rhocond_tot = p%rhop+rhovap
          if (lfirst .and. ldt) then
            np_total = 0.0
            total_surface_area = 0.0
            dt1_condensation = 0.0
            tau_phase1 = 0.0
            tau_evaporation1 = 0.0
            tau_evaporation = 0.0
          endif
          do k = k1_imn(imn),k2_imn(imn)
            ix0 = ineargrid(k,1)
            ix = ix0-nghost
!
!  The condensation/evaporation mass flux is
!
!    F = vth*rhovap*alpha
!
!  where vth is the thermal speed, rhopvap is the mass density of vapor,
!  and alpha \in [0,1] is the condensation coefficient. Various authors argue
!  for different choices of alpha. Following Barkstrom 1978 we take
!
!    alpha = 1/(vth*ap/D + 1/alpha0)
!
!  where D is the diffusion coefficient of vapor and ap the particle radius.
!  Small particles, with ap<D/(alpha0*vth), have alpha=alpha0, while the
!  condenation coefficient decreases linearly with particle size for larger
!  particles. Barkstrom 1978 use alpha0=0.04.
!
            select case (condensation_coefficient_type)
            case ('constant')
              alpha_cond_par = alpha_cond
            case ('size-dependent')
              alpha_cond_par = 1/(vth(ix)*fp(k,iap)*diffusion_coefficient1+ &
                  alpha_cond1)
            case default
              if (lroot) print*, 'dap_dt_condensation_pencil: '// &
                  'invalid condensation coefficient type'
              call fatal_error('dap_dt_condensation_pencil','')
              alpha_cond_par = 0.
            endselect
!
!  Radius increase by condensation or decrease by evaporation.
!
            if (fp(k,iap) < apmin) then
!
!  Do not allow particles to become smaller than a minimum radius.
!
              dapdt = -tau_damp_evap1*(fp(k,iap)-apmin)
              if (lfirst .and. ldt) then
                dt1_condensation(ix) = max(dt1_condensation(ix),tau_damp_evap1)
              endif
            else
!                    
              if (lcondensation_simplified) then
                if (ldust_accretion) then
                  dapdt = xi_accretion*p%rho(ix)/rho0
                else
                  dapdt = GS_condensation/fp(k,iap)
                endif
              elseif (lascalar) then
                if (ltauascalar) dapdt = G_condensation*f(ix,m,n,iacc)/fp(k,iap)
                if (lcondensation_rate) dapdt = G_condensation*f(ix,m,n,issat)/fp(k,iap)
                if (lcondensation_rate .and. ldust_condensation) dapdt = G_condensation*f(ix,m,n,issat)
              else
                dapdt = 0.25*vth(ix)*rhopmat1* &
                (rhovap(ix)-rhosat(ix))*alpha_cond_par
              endif
!            
!  Damp approach to minimum size. The radius decreases linearly with time in
!  the limit of small particles; therefore we need to damp the evaporation to
!  avoid small time-steps.
!
              if (dapdt < 0.0) then
                dapdt = dapdt*min(1.0,(fp(k,iap)/apmin-1.0)**2)
                if (lfirst .and. ldt) then
                  dt1_condensation(ix) = max(dt1_condensation(ix), &
                      abs(dapdt/(fp(k,iap)-apmin)))
                endif
              endif
            endif
!
            dfp(k,iap) = dfp(k,iap)+dapdt
!
!  Vapor monomers are added to the gas or removed from the gas.
!
            if (lparticles_number) np_swarm = fp(k,inpswarm)
            if (lcondensation_simplified .or. lcondensation_rate) then
              drhocdt = 0.
            else
              drhocdt = -dapdt*4*pi*fp(k,iap)**2*rhopmat*np_swarm
            endif
!
!  Drive the vapor pressure towards the saturated pressure due to contact
!  with "ocean" at the box bottom.
!
            if (lborder_driving_ocean) then
              if (fp(k,izp) < ztop_ocean) then
                drhocdt = drhocdt + tau_ocean_driving1*(rhosat(ix)-rhovap(ix))
              endif
            endif
!
!  feedback, but should not be used if we don't have density
!
            if (ldensity) then
              if (ldensity_nolog) then
                df(ix0,m,n,irho)   = df(ix0,m,n,irho)   + drhocdt
              else
                df(ix0,m,n,ilnrho) = df(ix0,m,n,ilnrho) + drhocdt*p%rho1(ix)
              endif
            endif
!
            if (lpscalar_nolog) then
              df(ix0,m,n,icc)   = df(ix0,m,n,icc)   + &
                  (1.0-p%cc(ix))*p%rho1(ix)*drhocdt
            elseif (lpscalar) then
              df(ix0,m,n,ilncc) = df(ix0,m,n,ilncc) + &
                  (p%cc1(ix)-1.0)*p%rho1(ix)*drhocdt
            endif
!
!  Release latent heat to gas / remove heat from gas.
!
            if (ltemperature .and. llatent_heat) then
              df(ix0,m,n,ilnTT) = df(ix0,m,n,ilnTT) - &
                  latent_heat_SI*p%rho1(ix)*p%TT1(ix)*p%cv1(ix)*drhocdt
            endif
!
!            if (lfirst .and. ldt) then
            if (lfirst .and. ldt .and. .not. lascalar .and. .not. lcondensation_simplified) then
              total_surface_area(ix) = total_surface_area(ix)+ &
                  4*pi*fp(k,iap)**2*np_swarm*alpha_cond_par
              np_total(ix) = np_total(ix)+np_swarm
            elseif (lfirst .and. ldt .and. lascalar .and. ldt_condensation) then
              tau_phase1(ix) = tau_phase1(ix)+4.0*pi*modified_vapor_diffusivity*fp(k,iap)*fp(k, inpswarm)
              if (ldt_evaporation) tau_evaporation1(ix) = -(2*G_condensation*f(ix,m,n,issat))/fp(k,iap)**2
            elseif (lfirst .and. ldt .and. lcondensation_simplified .and. ldt_condensation) then
              tau_phase1(ix) = tau_phase1(ix)+4.0*pi*modified_vapor_diffusivity*fp(k,iap)*fp(k, inpswarm)
            endif
          enddo
!
!  Time-step contribution of condensation.
!
!          if (lfirst .and. ldt) then
          if (lfirst .and. ldt .and. .not. lascalar .and. .not. lcondensation_simplified) then

            ap_equi = ((p%rhop+(rhovap-rhosat))/ &
                (4.0/3.0*pi*rhopmat*np_swarm*p%np))**(1.0/3.0)
            do ix = 1,nx
              if (rhocond_tot(ix) > rhosat(ix)) then
                dt1_condensation(ix) = max(total_surface_area(ix)*vth(ix), &
                    pi*vth(ix)*np_total(ix)*ap_equi(ix)**2*alpha_cond)
              endif
            enddo
          elseif (lfirst .and. ldt .and. lascalar .and. ldt_condensation) then
            do ix = 1,nx
              dt1_condensation(ix) = tau_phase1(ix)
              if (ldt_evaporation) dt1_condensation(ix) = max(tau_phase1(ix), tau_evaporation1(ix))
            enddo
          elseif (lfirst .and. ldt .and. lcondensation_simplified .and. ldt_condensation) then  
            do ix = 1,nx
              if (ldt_evaporation) then
                tau_evaporation1(ix) = -(2*GS_condensation)/ap0(1)**2
                dt1_condensation(ix) = max(tau_phase1(ix), tau_evaporation1(ix))
              else
                dt1_condensation(ix) = tau_phase1(ix)
              endif
            enddo
          elseif (lfirst .and. ldt .and. ldt_condensation_off) then
              dt1_condensation = 0.
          endif
        endif
!
        if (lborder_driving_ocean) then
          if (z(n) < ztop_ocean) then
            df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - &
                (1.0-TTocean*p%TT1)*tau_ocean_driving1
          endif
        endif
!
!  Time-step contribution of condensation.
!
        if (lfirst .and. ldt) then
          dt1_condensation = dt1_condensation/cdtpc
          dt1_max = max(dt1_max,dt1_condensation)
        endif
!
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine dap_dt_condensation_pencil
!***********************************************************************
    subroutine dap_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle radius.
!
!  21-nov-06/anders: coded
      use Sub
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mpar_loc,mpvar) :: dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
!  Diagnostic output.
!
      if (ldiagnos) then
        if (idiag_apm /= 0) call sum_par_name(fp(1:npar_loc,iap),idiag_apm)
        if (idiag_ap2m /= 0) call sum_par_name(fp(1:npar_loc,iap)**2,idiag_ap2m)
        if (idiag_apmin /= 0) &
            call max_par_name(-fp(1:npar_loc,iap),idiag_apmin,lneg=.true.)
        if (idiag_apmax /= 0) call max_par_name(fp(1:npar_loc,iap),idiag_apmax)
        if (idiag_npswarmm /= 0) &
            call sum_par_name(rhop_swarm/ &
            (four_pi_rhopmat_over_three*fp(1:npar_loc,iap)**3),idiag_npswarmm)
        if (idiag_ieffp /= 0) &
            call sum_par_name(fp(1:npar_loc,ieffp),idiag_ieffp)
      endif
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dap_dt
!***********************************************************************
    subroutine read_particles_rad_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
      integer :: pos
!
      read (parallel_unit, NML=particles_radius_init_pars, IOSTAT=iostat)
!
!  Find how many different particle radii we are using. This must be done
!  because not all parts of the code are adapted to work with more than one
!  particle radius.
!
      do pos = 1,ndustrad
        if (ap0(pos) /= 0) then
          npart_radii = npart_radii+1
        endif
      enddo
!
    endsubroutine read_particles_rad_init_pars
!***********************************************************************
    subroutine write_particles_rad_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write (unit, NML=particles_radius_init_pars)
!
    endsubroutine write_particles_rad_init_pars
!***********************************************************************
    subroutine read_particles_rad_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read (parallel_unit, NML=particles_radius_run_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_rad_run_pars
!***********************************************************************
    subroutine write_particles_rad_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write (unit, NML=particles_radius_run_pars)
!
    endsubroutine write_particles_rad_run_pars
!***********************************************************************
    subroutine rprint_particles_radius(lreset,lwrite)
!
!  Read and register print parameters relevant for particles radius.
!
!  22-aug-05/anders: coded
!
      use Diagnostics, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
!
!  Reset everything in case of reset.
!
      if (lreset) then
        idiag_apm = 0
        idiag_ap2m = 0
        idiag_apmin = 0
        idiag_apmax = 0
        idiag_dvp12m = 0
        idiag_dtsweepp = 0
        idiag_npswarmm = 0
        idiag_ieffp = 0
      endif
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot .and. ip < 14) &
          print*, 'rprint_particles_radius: run through parse list'
      do iname = 1,nname
        call parse_name(iname,cname(iname),cform(iname),'apm',idiag_apm)
        call parse_name(iname,cname(iname),cform(iname),'ap2m',idiag_ap2m)
        call parse_name(iname,cname(iname),cform(iname),'apmin',idiag_apmin)
        call parse_name(iname,cname(iname),cform(iname),'apmax',idiag_apmax)
        call parse_name(iname,cname(iname),cform(iname),'dvp12m',idiag_dvp12m)
        call parse_name(iname,cname(iname),cform(iname),'dtsweepp', &
            idiag_dtsweepp)
        call parse_name(iname,cname(iname),cform(iname),'ieffp',idiag_ieffp)
        if (.not. lparticles_number) call parse_name(iname,cname(iname), &
            cform(iname),'npswarmm',idiag_npswarmm)
      enddo
!
    endsubroutine rprint_particles_radius
!***********************************************************************
    subroutine get_stbin(api,iStbin)
      integer, intent(out) :: iStbin
      integer :: k=0
      real :: api
      k = 1
      if (lfixed_particles_radius) then
        do while((api  >=  ap0(k)).and.(k  <=  ndustrad))
          iStbin = k
          k = k+1
        enddo
      endif
    endsubroutine get_stbin
!***********************************************************************
    subroutine get_mass_from_radius(mpi,fpwn,ip)
      real, dimension (lpar_max,mparray) :: fpwn
      integer,intent(in) :: ip
      real,intent(out) :: mpi
      real :: api
      api = fpwn(ip,iap)
      mpi=(4./3.)*pi*rhopmat*(api**3)      
    endsubroutine get_mass_from_radius
!***********************************************************************
    subroutine get_maxrad(apmax)
      real :: apmax
      apmax=maxval(ap0)
    endsubroutine get_maxrad
!***********************************************************************
  endmodule Particles_radius
