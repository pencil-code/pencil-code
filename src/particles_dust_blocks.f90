! $Id$
!
!  This module takes care of everything related to dust particles.
!
!  This version is for block domain decomposition of particles. See
!  particles_map_blocks.f90 for documentation.
!
!  The file is relatively stripped down relative to particles_dust.f90, as
!  there is no guarantee that everything would have worked for block domain
!  decomposition. Feel free to add functionality from particles_dust.f90, but
!  be careful and test that it works properly for block domain decomposition.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 6
! MAUX CONTRIBUTION 2
! CPARAM logical, parameter :: lparticles=.true.
! CPARAM character (len=20), parameter :: particles_module="dust_blocks"
!
! PENCILS PROVIDED np; rhop; epsp; grhop(3);peh
! PENCILS PROVIDED uup(3)
!
!***************************************************************
module Particles
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
  use Particles_sub
  use Particles_radius
!
  implicit none
!
  include 'particles.h'
  include 'particles_common.h'
!
  real, dimension(mz) :: rho0z = 1.0
  real, dimension(npar_species), target :: tausp_species = 0.0, tausp1_species = 0.0
  real, dimension (3) :: pos_sphere=(/0.0,0.0,0.0/)
  real, dimension (3) :: pos_ellipsoid=(/0.0,0.0,0.0/)
  real, target :: tausp = 0.0
  real :: xp0=0.0, yp0=0.0, zp0=0.0, vpx0=0.0, vpy0=0.0, vpz0=0.0
  real :: xp1=0.0, yp1=0.0, zp1=0.0, vpx1=0.0, vpy1=0.0, vpz1=0.0
  real :: xp2=0.0, yp2=0.0, zp2=0.0, vpx2=0.0, vpy2=0.0, vpz2=0.0
  real :: xp3=0.0, yp3=0.0, zp3=0.0, vpx3=0.0, vpy3=0.0, vpz3=0.0
  real :: Lx0=0.0, Ly0=0.0, Lz0=0.0
  real :: delta_vp0=1.0, tausp1=0.0, tausp01=0.0
  real :: nu_epicycle=0.0, nu_epicycle2=0.0
  real :: beta_dPdr_dust=0.0, beta_dPdr_dust_scaled=0.0
  real :: epsp_friction_increase=0.0
  real :: cdtp=0.2, cdtpgrav=0.1
  real :: gravx=0.0, gravz=0.0, gravr=1.0, kx_gg=1.0, kz_gg=1.0
  real :: gravsmooth=0.0, gravsmooth2=0.0, Ri0=0.25, eps1=0.5
  real :: kx_xxp=0.0, ky_xxp=0.0, kz_xxp=0.0, amplxxp=0.0
  real :: kx_vvp=0.0, ky_vvp=0.0, kz_vvp=0.0, amplvvp=0.0
  real :: kx_vpx=0.0, kx_vpy=0.0, kx_vpz=0.0
  real :: ky_vpx=0.0, ky_vpy=0.0, ky_vpz=0.0
  real :: kz_vpx=0.0, kz_vpy=0.0, kz_vpz=0.0
  real :: phase_vpx=0.0, phase_vpy=0.0, phase_vpz=0.0
  real :: tstart_dragforce_par=0.0, tstart_grav_par=0.0
  real :: rad_sphere=0.0
  real :: a_ellipsoid=0.0, b_ellipsoid=0.0, c_ellipsoid=0.0
  real :: a_ell2=0.0, b_ell2=0.0, c_ell2=0.0
  real :: xsinkpoint=0.0, ysinkpoint=0.0, zsinkpoint=0.0, rsinkpoint=0.0
  real :: compensate_sedimentation=1. ! Is this still being used?
  real :: mean_free_path_gas=0.0
  real :: cs2_powerlaw
  logical, pointer :: lshearadvection_as_shift
  logical, target :: ldragforce_gas_par=.false.
  logical :: lcalc_uup = .false.
  logical :: ldragforce_dust_par=.false.
  logical :: lpar_spec=.false., learly_particle_map=.true.
  logical :: ldragforce_equi_global_eps=.false.
  logical :: ldraglaw_epstein=.true., ldraglaw_variable=.false.
  logical :: ldraglaw_variable_density=.false.
  logical :: ldraglaw_eps_stk_transonic=.false.,luse_tau_ap=.true.
  logical :: ldt_grav_par=.false., ldt_adv_par=.true.
  logical :: lsinkpoint=.false., lglobalrandom=.false.
  logical :: lcoriolis_force_par=.true., lcentrifugal_force_par=.false.
  logical :: lcylindrical_gravity_par=.false.
  logical :: lreassign_strat_rhom=.true.
  logical :: lcompensate_sedimentation=.false.
!
  character (len=labellen) :: interp_pol_uu ='ngp'
  character (len=labellen) :: interp_pol_oo ='ngp'
  character (len=labellen) :: interp_pol_TT ='ngp'
  character (len=labellen) :: interp_pol_rho='ngp'
!
  character (len=labellen), dimension (ninit) :: initxxp='nothing'
  character (len=labellen), dimension (ninit) :: initvvp='nothing'
  character (len=labellen) :: gravx_profile='', gravz_profile=''
  character (len=labellen) :: gravr_profile=''
!
  namelist /particles_init_pars/ &
      initxxp, initvvp, xp0, yp0, zp0, vpx0, vpy0, vpz0, delta_vp0, &
      ldragforce_gas_par, ldragforce_dust_par, bcpx, bcpy, bcpz, tausp, &
      beta_dPdr_dust, mpmat, np_swarm, rhop_swarm, eps_dtog, nu_epicycle, &
      rp_int, rp_ext, gravx_profile, gravz_profile, gravr_profile, gravx, &
      gravz, gravr, gravsmooth, kx_gg, kz_gg, Ri0, eps1, lmigration_redo, &
      ldragforce_equi_global_eps, kx_vvp, ky_vvp, kz_vvp, amplvvp, kx_xxp, &
      ky_xxp, kz_xxp, amplxxp, kx_vpx, kx_vpy, kx_vpz, ky_vpx, ky_vpy, ky_vpz, &
      kz_vpx, kz_vpy, kz_vpz, phase_vpx, phase_vpy, phase_vpz, &
      particle_mesh, lparticlemesh_cic, lparticlemesh_tsc, linterpolate_spline, &
      tstart_dragforce_par, tstart_grav_par, tausp_species, &
      learly_particle_map, epsp_friction_increase, lmigration_real_check, &
      ldraglaw_epstein, lcheck_exact_frontier, dustdensity_powerlaw, ldt_grav_par, &
      lsinkpoint,xsinkpoint, ysinkpoint, zsinkpoint, rsinkpoint, &
      lcoriolis_force_par, lcentrifugal_force_par, ldt_adv_par, Lx0, Ly0, &
      Lz0, lglobalrandom, linsert_particles_continuously, &
      remove_particle_at_time, remove_particle_criteria, remove_particle_criteria_size, &
      remove_particle_criteria_edtog, rad_sphere, pos_sphere, rhopmat, &
      a_ellipsoid, b_ellipsoid, c_ellipsoid, pos_ellipsoid, &
      lrandom_particle_pencils, it1_loadbalance, &
      lcalc_uup, lnocalc_np, lnocalc_rhop, lcommunicate_rhop, lcommunicate_np, &
      np_const, rhop_const, lrandom_particle_blocks, lreblock_particles_run, &
      lbrick_partition, ldraglaw_variable, ladopt_own_light_bricks, &
      xp1, yp1, zp1, vpx1, vpy1, vpz1, xp2, yp2, zp2, vpx2, vpy2, vpz2, &
      xp3, yp3, zp3, vpx3, vpy3, vpz3, lreassign_strat_rhom, &
      ldraglaw_eps_stk_transonic, luse_tau_ap, mean_free_path_gas
!
  namelist /particles_run_pars/ &
      bcpx, bcpy, bcpz, tausp, dsnap_par_minor, beta_dPdr_dust, &
      ldragforce_gas_par, ldragforce_dust_par, mp_swarm, np_swarm, &
      eps_dtog, cdtp, cdtpgrav, lpar_spec, linterp_reality_check, &
      nu_epicycle, gravx_profile, gravz_profile, gravr_profile, gravx, gravz, &
      gravr, gravsmooth, kx_gg, kz_gg, lmigration_redo, &
      tstart_dragforce_par, tstart_grav_par, &
      particle_mesh, lparticlemesh_cic, lparticlemesh_tsc, &
      epsp_friction_increase, learly_particle_map, lmigration_real_check, &
      ldraglaw_epstein, lcheck_exact_frontier, ldraglaw_variable_density, &
      remove_particle_at_time, remove_particle_criteria, remove_particle_criteria_size, &
      remove_particle_criteria_edtog, ldt_grav_par, lsinkpoint, rhopmat, &
      xsinkpoint, ysinkpoint, zsinkpoint, rsinkpoint, lcoriolis_force_par, &
      lcentrifugal_force_par, ldt_adv_par, linsert_particles_continuously, &
      lrandom_particle_pencils, lnocalc_np, lnocalc_rhop, it1_loadbalance, &
      np_const, rhop_const, lrandom_particle_blocks, lreblock_particles_run, &
      lbrick_partition, ldraglaw_variable, ladopt_own_light_bricks, &
      lcylindrical_gravity_par, lcommunicate_rhop, lcommunicate_np, &
      lcompensate_sedimentation,compensate_sedimentation, &
      ldraglaw_eps_stk_transonic, luse_tau_ap, mean_free_path_gas
!
  integer :: idiag_xpm=0, idiag_ypm=0, idiag_zpm=0
  integer :: idiag_xp2m=0, idiag_yp2m=0, idiag_zp2m=0
  integer :: idiag_rpm=0, idiag_rp2m=0
  integer :: idiag_vpxm=0, idiag_vpym=0, idiag_vpzm=0
  integer :: idiag_vpx2m=0, idiag_vpy2m=0, idiag_vpz2m=0, idiag_ekinp=0
  integer :: idiag_vpxmax=0, idiag_vpymax=0, idiag_vpzmax=0, idiag_vpmax=0
  integer :: idiag_lpxm=0, idiag_lpym=0, idiag_lpzm=0
  integer :: idiag_lpx2m=0, idiag_lpy2m=0, idiag_lpz2m=0
  integer :: idiag_npm=0, idiag_np2m=0, idiag_npmax=0, idiag_npmin=0
  integer :: idiag_dtdragp=0
  integer :: idiag_nparmin=0, idiag_nparmax=0, idiag_npargone=0
  integer :: idiag_rhopm=0, idiag_rhoprms=0, idiag_rhop2m=0, idiag_rhopmax=0
  integer :: idiag_rhopmin=0, idiag_decollp=0, idiag_rhopmphi=0
  integer :: idiag_epspmin=0, idiag_epspmax=0
  integer :: idiag_npmx=0, idiag_npmy=0, idiag_npmz=0
  integer :: idiag_rhopmx=0, idiag_rhopmy=0, idiag_rhopmz=0
  integer :: idiag_epspmx=0, idiag_epspmy=0, idiag_epspmz=0
  integer :: idiag_rpvpxmx = 0, idiag_rpvpymx = 0, idiag_rpvpzmx = 0
  integer :: idiag_rpvpx2mx = 0, idiag_rpvpy2mx = 0, idiag_rpvpz2mx = 0
  integer :: idiag_mpt=0, idiag_dedragp=0
  integer :: idiag_rhopmxy=0, idiag_rhopmxz=0, idiag_rhopmr=0
  integer :: idiag_sigmap=0
  integer :: idiag_dvpx2m=0, idiag_dvpy2m=0, idiag_dvpz2m=0
  integer :: idiag_dvpm=0, idiag_dvpmax=0, idiag_epotpm=0
  integer :: idiag_nparbmax=0, idiag_nblockmin=0, idiag_nblockmax=0
  integer :: idiag_eccpxm=0, idiag_eccpym=0, idiag_eccpzm=0
  integer :: idiag_eccpx2m=0, idiag_eccpy2m=0, idiag_eccpz2m=0
  integer :: idiag_deshearbcsm=0
!
  real :: gamma

  contains
!***********************************************************************
    subroutine register_particles()
!
!  Set up indices for access to the fp and dfp arrays
!
!  29-dec-04/anders: coded
!
      use FArrayManager
      use SharedVariables, only: put_shared_variable
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Indices for particle position.
!
      call append_npvar('ixp',ixp)
      call append_npvar('iyp',iyp)
      call append_npvar('izp',izp)
!
!  Indices for particle velocity.
!
      call append_npvar('ivpx',ivpx)
      call append_npvar('ivpy',ivpy)
      call append_npvar('ivpz',ivpz)
!
!  Set indices for auxiliary variables
!
      if (lcalc_uup) then
        call farray_register_auxiliary("uup", iuup, communicated=.true., vector=3)
        iupx = iuup
        iupy = iuup + 1
        iupz = iuup + 2
      endif
!
      if (.not. lnocalc_np) call farray_register_auxiliary('np',inp,communicated=lcommunicate_np)
!
      if (.not. lnocalc_rhop) call farray_register_auxiliary('rhop',irhop, &
          communicated=lparticles_sink.or.lcommunicate_rhop)
!
!  Share friction time (but only if Epstein drag regime!).
!
      call put_shared_variable("ldragforce_gas_par", ldragforce_gas_par, caller="register_particles")
!
      if (ldraglaw_epstein) then
        call put_shared_variable("tausp", tausp)
        call put_shared_variable("tausp_species", tausp_species)
        call put_shared_variable("tausp1_species", tausp1_species)
      endif 
!
!  Share Keplerian gravity.
!
      call put_shared_variable('gravr',gravr,caller='register_particles')
!
    endsubroutine register_particles
!***********************************************************************
    subroutine initialize_particles(f,fp)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  29-dec-04/anders: coded
!   5-mar-15/MR: reference state included in calculation of mean density
!
      use EquationOfState, only: cs0, rho0, get_stratz, get_gamma_etc
      use FArrayManager
      use SharedVariables, only: get_shared_variable
      use Density, only: mean_density
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray), intent (in) :: fp
!
      real :: rhom
      integer :: ierr, jspec
      real, pointer :: reference_state_mass
      real, pointer :: temperature_power_law
!
!  This module is incompatible with normal domain decomposition.
!
      if (.not.lparticles_blocks) &
        call fatal_error('initialize_particles', &
                         'must use PARTICLES=PARTICLES_DUST for normal domain decomposition')

      call get_gamma_etc(gamma)
!
!  The inverse stopping time is needed for drag force and collisional cooling.
!
      if (tausp/=0.0) tausp1=1/tausp
!
!  Get density stratification.
!
      if (lstratz) then
        call get_stratz(z, rho0z)
      else
        rho0z = rho0
      endif
!
!  For drag force calculation we need to fill blocks with information about
!  the gas density field and the gas velocity field.
!
      if (lhydro.and.ldragforce_dust_par) lfill_blocks_velocity=.true.
      if (ldragforce_gas_par)  lfill_blocks_density=.true.
      if (ldragforce_gas_par)  lfill_bricks_velocity=.true.
!
!  Multiple dust species. Friction time is given in the array tausp_species.
!
      if (npar_species>1) then
        if (lroot) then
          print*, 'initialize_particles: Number of particle species = ', npar_species
          print*, 'initialize_particles: tausp_species = ', tausp_species
        endif
!
!  Must have set tausp_species for drag force.
!
        if (ldragforce_dust_par .or. ldragforce_gas_par) then
          if (any(tausp_species==0.0)) &
                call fatal_error('initialize_particles','drag force must have tausp_species/=0')
!
!  Inverse friction time is needed for drag force.
!
          do jspec=1,npar_species
            if (tausp_species(jspec)/=0.0) tausp1_species(jspec)=1/tausp_species(jspec)
          enddo
        endif
      else
        tausp_species(1)=tausp
        if (tausp_species(1)/=0.0) tausp1_species(1)=1/tausp_species(1)
      endif
!
      if (ldraglaw_eps_stk_transonic.and.ldraglaw_epstein) &
        call fatal_error('initialize_particles', &
             'both epstein and epstein-stokes transonic drag laws are on - only one allowed')
!
      if (ldraglaw_eps_stk_transonic) then
        if (llocal_iso) then
          call get_shared_variable('itemperature_power_law',temperature_power_law)
          cs2_powerlaw=temperature_power_law
        endif
      endif
!
!  Normalization for Epstein drag with particles radius.
!
      if (ldraglaw_epstein .and. iap /= 0) tausp01 = rho0 * cs0 / (rhopmat * Omega)
!
!  Global gas pressure gradient seen from the perspective of the dust.
!
      if (beta_dPdr_dust/=0.0) then
        beta_dPdr_dust_scaled=beta_dPdr_dust*Omega/cs0
        if (lroot) print*, 'initialize_particles: Global pressure '// &
            'gradient with beta_dPdr_dust=', beta_dPdr_dust
      endif
!
!  Calculate mass density per particle (for back-reaction drag force on gas)
!  following the formula
!    rhop_swarm*N_cell = eps*rhom
!  where rhop_swarm is the mass density per superparticle, N_cell is the number
!  of particles per grid cell and rhom is the mean gas density in the box.
!
      if (rhop_swarm==0.0 .or. mp_swarm==0.0) then
! For stratification, take into account gas present outside the simulation box.
        if (lreassign_strat_rhom.and.((lgravz .and. lgravz_gas) .or. gravz_profile=='linear')) then
          ! rhom = (total mass) / (box volume) = Sigma / Lz
          ! Sigma = sqrt(2pi) * rho0 * H
          !   rho0 = mid-plane density, H = (sound speed) / (epicycle frequency)
          rhom = sqrt(2.0 * pi) / Lz
          if (nu_epicycle > 0.0) rhom = rhom * (rho0 * cs0 / nu_epicycle)
        else
          rhom = rho0
          if (lreference_state) then
            call get_shared_variable('reference_state_mass',reference_state_mass, &
                                     caller='initialize_particles')
            rhom=rhom+reference_state_mass/box_volume
          endif
        endif
        if (rhop_swarm==0.0) rhop_swarm = eps_dtog*rhom/(real(npar)/(nxgrid*nygrid*nzgrid))
        if (mp_swarm==0.0) mp_swarm   = eps_dtog*rhom*box_volume/(real(npar))
        if (lroot) print*, 'initialize_particles: dust-to-gas ratio eps_dtog=', eps_dtog
      endif
!
      if (lroot) then
        print*, 'initialize_particles: mass per constituent particle mpmat=', mpmat
        print*, 'initialize_particles: mass per superparticle mp_swarm =', mp_swarm
        print*, 'initialize_particles: number density per superparticle np_swarm=', np_swarm
        print*, 'initialize_particles: mass density per superparticle rhop_swarm=', rhop_swarm
      endif
!
!  Calculate nu_epicycle**2 for gravity.
!
      if (gravz_profile=='' .and. nu_epicycle/=0.0) gravz_profile='linear'
      nu_epicycle2=nu_epicycle**2
!
!  Calculate gravsmooth**2 for gravity.
!
      if (gravsmooth/=0.0) gravsmooth2=gravsmooth**2
!
!  Gas density is needed for back-reaction friction force.
!
      if (ldragforce_gas_par .and. .not. ldensity) &
        call fatal_error('initialize_particles', &
                         'friction force on gas only works together with gas density module')
!
!  Need to map particles on the grid for dragforce on gas.
!
      if (ldragforce_gas_par) then
!
!  When drag force is smoothed, df is also set in the first ghost zone. This
!  region needs to be folded back into the df array after pde is finished,
!
        if (lparticlemesh_cic .or. lparticlemesh_tsc) lfold_df=.true.
      endif
!
!  Initialize storage of energy gain released by shearing boundaries.
!
      if (idiag_deshearbcsm/=0) energy_gain_shear_bcs=0.0
!
!  Drag force on gas right now assumed rhop_swarm is the same for all particles.
!
      if (ldragforce_gas_par.and.(lparticles_radius.or.lparticles_number) &
          .and..not.lparticles_density) &
        call not_implemented('initialize_particles', &
                             'drag force on gas for variable particle radius or number')
!
!  Check if shear advection is on for time-step condition.
!
      if (lshear) call get_shared_variable('lshearadvection_as_shift', lshearadvection_as_shift)
!
!  Write constants to disk.
!
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro',position='append')
          write (1,*) 'np_swarm=', np_swarm
          write (1,*) 'mpmat=', mpmat
          write (1,*) 'mp_swarm=', mp_swarm
          write (1,*) 'rhop_swarm=', rhop_swarm
        close (1)
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles
!***********************************************************************
    subroutine init_particles(f,fp,ineargrid)
!
!  Initial positions and velocities of dust particles.
!
!  29-dec-04/anders: coded
!
      use EquationOfState, only: cs20
      use SharedVariables, only: get_shared_variable
      use General, only: random_number_wrapper
      use Mpicomm, only: mpireduce_sum, mpibcast_real
      use Sub
      use InitialCondition, only: initial_condition_xxp,initial_condition_vvp
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      intent (out) :: f, fp, ineargrid
!
      real, dimension (3) :: Lxyz_par, xyz0_par, xyz1_par
      real :: vpx_sum, vpy_sum, vpz_sum
      real :: r, p, q, px, py, pz, eps, cs, k2_xxp, rp2
      real :: dim1, npar_loc_x, npar_loc_y, npar_loc_z, dx_par, dy_par, dz_par
      real :: rad,rad_scl,tht,phi,tmp,OO,xx0,yy0,r2
      integer :: l, j, k, ix0, iy0, iz0, ib
      logical :: lequidistant=.false.
!
      real, dimension(:), pointer :: beta_glnrho_global
!
      call get_shared_variable('beta_glnrho_global',beta_glnrho_global,caller='init_particles')
!
!  Use either a local random position or a global random position for certain
!  initial conditions. The default is a local random position, but the equal
!  number of particles per processors means that this is not completely random.
!
      if (lglobalrandom) then
        Lxyz_par=Lxyz
        xyz0_par=xyz0
        xyz1_par=xyz1
      else
        Lxyz_par=Lxyz_loc
        xyz0_par=xyz0_loc
        xyz1_par=xyz1_loc
      endif
!
!  Initial particle position.
!
      do j=1,ninit
!
        select case (initxxp(j))
!
        case ('nothing')
          if (lroot .and. j==1) print*, 'init_particles: nothing'
!
        case ('origin')
          if (lroot) print*, 'init_particles: All particles at origin'
          fp(1:npar_loc,ixp:izp)=0.
!
        case ('zero-z')
          if (lroot) print*, 'init_particles: Zero z coordinate'
          fp(1:npar_loc,izp)=0.
!
        case ('constant')
          if (lroot) print*, 'init_particles: All particles at x,y,z=', xp0, yp0, zp0
          fp(1:npar_loc,ixp)=xp0
          fp(1:npar_loc,iyp)=yp0
          fp(1:npar_loc,izp)=zp0
!
        case ('constant-1')
          if (lroot) print*, 'init_particles: Particle 1 at x,y,z=', xp1, yp1, zp1
          do k=1,npar_loc
            if (ipar(k)==1) then
              fp(k,ixp)=xp1
              fp(k,iyp)=yp1
              fp(k,izp)=zp1
            endif
          enddo
!
        case ('constant-2')
          if (lroot) print*, 'init_particles: Particle 2 at x,y,z=', xp2, yp2, zp2
          do k=1,npar_loc
            if (ipar(k)==2) then
              fp(k,ixp)=xp2
              fp(k,iyp)=yp2
              fp(k,izp)=zp2
            endif
          enddo
!
        case ('constant-3')
          if (lroot) print*, 'init_particles: Particle 2 at x,y,z=', xp3, yp3, zp3
          do k=1,npar_loc
            if (ipar(k)==3) then
              fp(k,ixp)=xp3
              fp(k,iyp)=yp3
              fp(k,izp)=zp3
            endif
          enddo
!
        case ('random')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k=1,npar_loc
            if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
            if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
            if (nzgrid/=1) call random_number_wrapper(fp(k,izp))
          enddo
          if (nxgrid/=1) &
              fp(1:npar_loc,ixp)=xyz0_par(1)+fp(1:npar_loc,ixp)*Lxyz_par(1)
          if (nygrid/=1) &
              fp(1:npar_loc,iyp)=xyz0_par(2)+fp(1:npar_loc,iyp)*Lxyz_par(2)
          if (nzgrid/=1) &
              fp(1:npar_loc,izp)=xyz0_par(3)+fp(1:npar_loc,izp)*Lxyz_par(3)
!
        case ('random-circle')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k=1,npar_loc
            call random_number_wrapper(r)
            if (zp0>yp0) then
              fp(k,ixp)=xp0*cos((zp0-yp0)*r+yp0)
              fp(k,iyp)=xp0*sin((zp0-yp0)*r+yp0)
            else
              fp(k,ixp)=xp0*cos(2*pi*r)
              fp(k,iyp)=xp0*sin(2*pi*r)
            endif
          enddo
!
        case ('random-sphere')
          if (lroot) print*, 'init_particles: Random particle positions '// &
              'in a sphere around ', pos_sphere, ' with radius=',rad_sphere
          if (rad_sphere==0) then
            call fatal_error('init_particles','random-sphere '// &
                  'radius needs to be larger than zero')
          endif
          if (-rad_sphere+pos_sphere(1)<xyz0(1) .or. &
               rad_sphere+pos_sphere(1)>xyz1(1) .or. &
              -rad_sphere+pos_sphere(2)<xyz0(2) .or. &
               rad_sphere+pos_sphere(2)>xyz1(2) .or. &
              -rad_sphere+pos_sphere(3)<xyz0(3) .or. &
               rad_sphere+pos_sphere(3)>xyz1(3)) then
            call fatal_error('init_particles','random-sphere needs to fit in the box')
          endif
          if (lcartesian_coords) then
            do k=1,npar_loc
              rp2=2.*rad_sphere**2
              do while (rp2>rad_sphere**2)
                call random_number_wrapper(fp(k,ixp))
                call random_number_wrapper(fp(k,iyp))
                call random_number_wrapper(fp(k,izp))
                fp(k,ixp)=(fp(k,ixp)-0.5)*2.*rad_sphere
                fp(k,iyp)=(fp(k,iyp)-0.5)*2.*rad_sphere
                fp(k,izp)=(fp(k,izp)-0.5)*2.*rad_sphere
                rp2=fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2
              enddo
              fp(k,ixp)=fp(k,ixp)+pos_sphere(1)
              fp(k,iyp)=fp(k,iyp)+pos_sphere(2)
              fp(k,izp)=fp(k,izp)+pos_sphere(3)
            enddo
          else
            call fatal_error('init_particles','random-sphere '// &
                 'only implemented for cartesian coordinates')
          endif
!
        case ('random-ellipsoid')
          if (lroot) print*, 'init_particles: Random particle positions '// &
              'in an ellipsoid around ', pos_ellipsoid, ' with ' // &
              'semi-principal axes a,b,c =',a_ellipsoid,b_ellipsoid,c_ellipsoid
          if ((a_ellipsoid==0) .or. (b_ellipsoid==0) .or. (c_ellipsoid==0)) then
            call fatal_error('init_particles','random-ellipsoid '// &
                'all semi-principal axes need to be larger than zero')
          endif
          if (-a_ellipsoid+pos_ellipsoid(1)<xyz0(1) .or. &
               a_ellipsoid+pos_ellipsoid(1)>xyz1(1) .or. &
              -b_ellipsoid+pos_ellipsoid(2)<xyz0(2) .or. &
               b_ellipsoid+pos_ellipsoid(2)>xyz1(2) .or. &
              -c_ellipsoid+pos_ellipsoid(3)<xyz0(3) .or. &
               c_ellipsoid+pos_ellipsoid(3)>xyz1(3)) then
            call fatal_error('init_particles','random-ellipsoid needs to fit in the box')
          endif
          if (lcartesian_coords) then
            a_ell2=a_ellipsoid**2
            b_ell2=b_ellipsoid**2
            c_ell2=c_ellipsoid**2
            do k=1,npar_loc
              rp2=2.
              do while (rp2>1.)
                call random_number_wrapper(fp(k,ixp))
                call random_number_wrapper(fp(k,iyp))
                call random_number_wrapper(fp(k,izp))
                fp(k,ixp)=(fp(k,ixp)-0.5)*2.*a_ellipsoid
                fp(k,iyp)=(fp(k,iyp)-0.5)*2.*b_ellipsoid
                fp(k,izp)=(fp(k,izp)-0.5)*2.*c_ellipsoid
                rp2=fp(k,ixp)**2/a_ell2+fp(k,iyp)**2/b_ell2+fp(k,izp)**2/c_ell2
              enddo
              fp(k,ixp)=fp(k,ixp)+pos_ellipsoid(1)
              fp(k,iyp)=fp(k,iyp)+pos_ellipsoid(2)
              fp(k,izp)=fp(k,izp)+pos_ellipsoid(3)
            enddo
          else
            call fatal_error('init_particles','random-ellipsoid '// &
                 'only implemented for cartesian coordinates')
          endif
!
        case ('random-hole')
          if (lroot) print*, 'init_particles: Random particle positions with inner hole'
          do k=1,npar_loc
            rp2=-1.0
            do while (rp2<rp_int**2)
              if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
              if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
              if (nzgrid/=1) call random_number_wrapper(fp(k,izp))
              if (nxgrid/=1) fp(k,ixp)=xyz0(1)+fp(k,ixp)*Lxyz(1)
              if (nygrid/=1) fp(k,iyp)=xyz0(2)+fp(k,iyp)*Lxyz(2)
              if (nzgrid/=1) fp(k,izp)=xyz0(3)+fp(k,izp)*Lxyz(3)
              rp2=fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2
            enddo
          enddo
!
        case ('random-box')
          if (lroot) print*, 'init_particles: Random particle positions within a box'
          do k=1,npar_loc
            if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
            if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
            if (nzgrid/=1) call random_number_wrapper(fp(k,izp))
            if (lcylindrical_coords) then
              xx0=xp0+fp(k,ixp)*Lx0
              yy0=yp0+fp(k,iyp)*Ly0
              r2=xx0**2+yy0**2
              if (nxgrid/=1) fp(k,ixp)=sqrt(r2)
              if (nygrid/=1) fp(k,iyp)=atan(yy0/xx0)+pi*(xx0/abs(xx0)-1)*0.5
              if (nzgrid/=1) fp(k,izp)=zp0+fp(k,izp)*Lz0
            else
              if (nxgrid/=1) fp(k,ixp)=xp0+fp(k,ixp)*Lx0
              if (nygrid/=1) fp(k,iyp)=yp0+fp(k,iyp)*Ly0
              if (nzgrid/=1) fp(k,izp)=zp0+fp(k,izp)*Lz0
            endif
          enddo
!
       case ('random-cylindrical','random-cyl','random-spherical','random-sph')
          if (lroot) print*, 'init_particles: Random particle '// &
               'cylindrical positions with power-law =',dustdensity_powerlaw
!
          do k=1,npar_loc
!
! Start the particles obeying a power law
!
            if (lcylindrical_coords.or.lcartesian_coords) then
              tmp=2-dustdensity_powerlaw
            elseif (lspherical_coords) then
              tmp=3-dustdensity_powerlaw
            else
              call fatal_error("init_particles","The world is flat, and we never got here")
            endif
            call random_number_wrapper(rad_scl)
            rad_scl = rp_int**tmp + rad_scl*(rp_ext**tmp-rp_int**tmp)
            rad = rad_scl**(1./tmp)
!
! Random in azimuth
!
            call random_number_wrapper(phi)
!
            if (lcartesian_coords) then
              phi = 2*pi*phi
              if (nxgrid/=1) fp(k,ixp)=rad*cos(phi)
              if (nygrid/=1) fp(k,iyp)=rad*sin(phi)
            elseif (lcylindrical_coords) then
              phi = xyz0_par(2)+phi*Lxyz_par(2)
              if (nxgrid/=1) fp(k,ixp)=rad
              if (nygrid/=1) fp(k,iyp)=phi
            elseif (lspherical_coords) then
              if (nxgrid/=1) fp(k,ixp)=rad
              if (nygrid/=1) then
                call random_number_wrapper(tht)
                fp(k,iyp) = xyz0_par(2)+tht*Lxyz_par(2)
              endif
              if (nzgrid/=1) then
                call random_number_wrapper(phi)
                fp(k,izp) = xyz0_par(3)+phi*Lxyz_par(3)
              endif
            endif
!
            if (nzgrid/=1) call random_number_wrapper(fp(k,izp))
            if (nzgrid/=1) fp(k,izp)=xyz0_par(3)+fp(k,izp)*Lxyz_par(3)
!
          enddo
!
        case ('np-constant')
          if (lroot) print*, 'init_particles: Constant number density'
          k=1
k_loop:   do while (.not. (k>npar_loc))
            do l=l1,l2; do m=m1,m2; do n=n1,n2
              if (nxgrid/=1) call random_number_wrapper(px)
              if (nygrid/=1) call random_number_wrapper(py)
              if (nzgrid/=1) call random_number_wrapper(pz)
              fp(k,ixp)=x(l)+(px-0.5)*dx
              fp(k,iyp)=y(m)+(py-0.5)*dy
              fp(k,izp)=z(n)+(pz-0.5)*dz
              k=k+1
              if (k>npar_loc) exit k_loop
            enddo; enddo; enddo
          enddo k_loop
!
        case ('equidistant')
          if (lroot) print*, 'init_particles: Particles placed equidistantly'
          dim1=1.0/dimensionality
!
!  Number of particles per direction. Found by solving the equation system
!
!    npar_loc_x/npar_loc_y = Lx_loc/Ly_loc
!    npar_loc_x/npar_loc_z = Lx_loc/Lz_loc
!    npar_loc_y/npar_loc_z = Ly_loc/Lz_loc
!    npar_loc_x*npar_loc_y*npar_loc_z = npar_loc
!
!  Found it to be easier to separate in all possible dimensionalities.
!  For a missing direction i, set npar_loc_i=1 in the above equations and
!  ignore any equation that has Li_loc in it.
!
          if (dimensionality==3) then
!  3-D
            npar_loc_x=(npar_loc*Lxyz_loc(1)**2/(Lxyz_loc(2)*Lxyz_loc(3)))**dim1
            npar_loc_y=(npar_loc*Lxyz_loc(2)**2/(Lxyz_loc(1)*Lxyz_loc(3)))**dim1
            npar_loc_z=(npar_loc*Lxyz_loc(3)**2/(Lxyz_loc(1)*Lxyz_loc(2)))**dim1
          elseif (dimensionality==2) then
!  2-D
            if (nxgrid==1) then
              npar_loc_x=1
              npar_loc_y=(npar_loc*Lxyz_loc(2)/Lxyz_loc(3))**dim1
              npar_loc_z=(npar_loc*Lxyz_loc(3)/Lxyz_loc(2))**dim1
            elseif (nygrid==1) then
              npar_loc_x=(npar_loc*Lxyz_loc(1)/Lxyz_loc(3))**dim1
              npar_loc_y=1
              npar_loc_z=(npar_loc*Lxyz_loc(3)/Lxyz_loc(2))**dim1
            elseif (nzgrid==1) then
              npar_loc_x=(npar_loc*Lxyz_loc(1)/Lxyz_loc(2))**dim1
              npar_loc_y=(npar_loc*Lxyz_loc(2)/Lxyz_loc(1))**dim1
              npar_loc_z=1
            endif
          elseif (dimensionality==1) then
!  1-D
            if (nxgrid/=1) then
              npar_loc_x=npar_loc
              npar_loc_y=1
              npar_loc_z=1
            elseif (nygrid/=1) then
              npar_loc_x=1
              npar_loc_y=npar_loc
              npar_loc_z=1
            elseif (nzgrid/=1) then
              npar_loc_x=1
              npar_loc_y=1
              npar_loc_z=npar_loc
            endif
          endif
!  Distance between particles.
          dx_par=Lxyz_loc(1)/npar_loc_x
          dy_par=Lxyz_loc(2)/npar_loc_y
          dz_par=Lxyz_loc(3)/npar_loc_z
!  Place first particle.
          fp(1,ixp) = x(l1) ; fp(1,iyp) = y(m1) ; fp(1,izp) = z(n1)
          if (nxgrid/=1) fp(1,ixp) = xyz0_loc(1)+dx_par/2
          if (nygrid/=1) fp(1,iyp) = xyz0_loc(2)+dy_par/2
          if (nzgrid/=1) fp(1,izp) = xyz0_loc(3)+dz_par/2
!  Place all other particles iteratively.
          if (dimensionality==3) then
!  3-D
            do k=2,npar_loc
              fp(k,ixp)=fp(k-1,ixp)+dx_par
              fp(k,iyp)=fp(k-1,iyp)
              fp(k,izp)=fp(k-1,izp)
              if (fp(k,ixp)>xyz1_loc(1)) then
                fp(k,ixp)=fp(1,ixp)
                fp(k,iyp)=fp(k,iyp)+dy_par
              endif
              if (fp(k,iyp)>xyz1_loc(2)) then
                fp(k,iyp)=fp(1,iyp)
                fp(k,izp)=fp(k,izp)+dz_par
              endif
            enddo
          elseif (dimensionality==2) then
!  2-D
            if (nxgrid==1) then
              do k=2,npar_loc
                fp(k,ixp)=fp(k-1,ixp)
                fp(k,iyp)=fp(k-1,iyp)+dy_par
                fp(k,izp)=fp(k-1,izp)
                if (fp(k,iyp)>xyz1_loc(2)) then
                  fp(k,iyp)=fp(1,iyp)
                  fp(k,izp)=fp(k,izp)+dz_par
                endif
              enddo
            elseif (nygrid==1) then
              do k=2,npar_loc
                fp(k,ixp)=fp(k-1,ixp)+dx_par
                fp(k,iyp)=fp(k-1,iyp)
                fp(k,izp)=fp(k-1,izp)
                if (fp(k,ixp)>xyz1_loc(1)) then
                  fp(k,ixp)=fp(1,ixp)
                  fp(k,izp)=fp(k,izp)+dz_par
                endif
              enddo
            elseif (nzgrid==1) then
              do k=2,npar_loc
                fp(k,ixp)=fp(k-1,ixp)+dx_par
                fp(k,iyp)=fp(k-1,iyp)
                fp(k,izp)=fp(k-1,izp)
                if (fp(k,ixp)>xyz1_loc(1)) then
                  fp(k,ixp)=fp(1,ixp)
                  fp(k,iyp)=fp(k,iyp)+dy_par
                endif
              enddo
            endif
          elseif (dimensionality==1) then
!  1-D
            if (nxgrid/=1) then
              do k=2,npar_loc
                fp(k,ixp)=fp(k-1,ixp)+dx_par
                fp(k,iyp)=fp(k-1,iyp)
                fp(k,izp)=fp(k-1,izp)
              enddo
            elseif (nygrid/=1) then
              do k=2,npar_loc
                fp(k,ixp)=fp(k-1,ixp)
                fp(k,iyp)=fp(k-1,iyp)+dy_par
                fp(k,izp)=fp(k-1,izp)
              enddo
            elseif (nzgrid/=1) then
              do k=2,npar_loc
                fp(k,ixp)=fp(k-1,ixp)
                fp(k,iyp)=fp(k-1,iyp)
                fp(k,izp)=fp(k-1,izp)+dz_par
              enddo
            endif
          else
!  0-D
            fp(2:npar_loc,ixp)=fp(1,ixp)
            fp(2:npar_loc,iyp)=fp(1,iyp)
            fp(2:npar_loc,izp)=fp(1,izp)
          endif
          lequidistant=.true.
!
!  Shift particle locations slightly so that a mode appears.
!
        case ('shift')
          if (lroot) print*, 'init_particles: shift particle positions'
          if (.not. lequidistant) &
            call fatal_error('init_particles','must place particles equidistantly before shifting')
          k2_xxp=kx_xxp**2+ky_xxp**2+kz_xxp**2
          if (k2_xxp==0.0) &
            call fatal_error('init_particles','kx_xxp=ky_xxp=kz_xxp=0.0 is not allowed')
          do k=1,npar_loc
            fp(k,ixp) = fp(k,ixp) - kx_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp)+ky_xxp*fp(k,iyp)+kz_xxp*fp(k,izp))
            fp(k,iyp) = fp(k,iyp) - ky_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp)+ky_xxp*fp(k,iyp)+kz_xxp*fp(k,izp))
            fp(k,izp) = fp(k,izp) - kz_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp)+ky_xxp*fp(k,iyp)+kz_xxp*fp(k,izp))
          enddo
!
!  Shift to egg crate mode 2d, cos(x)cos(z)
!
        case ('cosxcosz')
          if (lroot) print*, 'init_particles: shift particle positions'
          if (.not. lequidistant) &
            call fatal_error('init_particles','must place particles equidistantly before shifting')
          k2_xxp=kx_xxp**2+kz_xxp**2
          if (k2_xxp==0.0) &
            call fatal_error('init_particles','kx_xxp=ky_xxp=kz_xxp=0.0 is not allowed')
          do k=1,npar_loc
            fp(k,ixp) = fp(k,ixp) - kx_xxp/k2_xxp*amplxxp*sin(kx_xxp*fp(k,ixp))*cos(kz_xxp*fp(k,izp))
            fp(k,izp) = fp(k,izp) - kz_xxp/k2_xxp*amplxxp*sin(kx_xxp*fp(k,ixp))*cos(kz_xxp*fp(k,izp))
          enddo
!
!  Shift to egg crate mode 2d, sin(x)sin(z)
!
        case ('sinxsinz')
          if (lroot) print*, 'init_particles: shift particle positions'
          if (.not. lequidistant) &
            call fatal_error('init_particles','must place particles equidistantly before shifting')
          k2_xxp=kx_xxp**2+kz_xxp**2
          if (k2_xxp==0.0) &
            call fatal_error('init_particles','kx_xxp=ky_xxp=kz_xxp=0.0 is not allowed')
          do k=1,npar_loc
          fp(k,ixp) = fp(k,ixp) + kx_xxp/k2_xxp*amplxxp*cos(kx_xxp*fp(k,ixp))*sin(kz_xxp*fp(k,izp))
          fp(k,izp) = fp(k,izp) + kz_xxp/k2_xxp*amplxxp*cos(kx_xxp*fp(k,ixp))*sin(kz_xxp*fp(k,izp))
          enddo
!
        case ('gaussian-z')
          if (lroot) print*, 'init_particles: Gaussian particle positions'
          do k=1,npar_loc
            do while (.true.)
              if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
              if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
              call random_number_wrapper(r)
              call random_number_wrapper(p)
              if (nprocz==2) then
                if (lfirst_proc_z) fp(k,izp)=-abs(zp0*sqrt(-2*alog(r))*cos(2*pi*p))
                if (llast_proc_z) fp(k,izp)=abs(zp0*sqrt(-2*alog(r))*cos(2*pi*p))
              else
                fp(k,izp)= zp0*sqrt(-2*alog(r))*cos(2*pi*p)
              endif
              if ((fp(k,izp)>=xyz0(3)).and.(fp(k,izp)<=xyz1(3))) exit
            enddo
          enddo
          if (nxgrid/=1) &
              fp(1:npar_loc,ixp)=xyz0_par(1)+fp(1:npar_loc,ixp)*Lxyz_par(1)
          if (nygrid/=1) &
              fp(1:npar_loc,iyp)=xyz0_par(2)+fp(1:npar_loc,iyp)*Lxyz_par(2)
!
        case ('gaussian-z-pure')
          if (lroot) print*, 'init_particles: Gaussian particle positions'
          do k=1,npar_loc
            call random_number_wrapper(r)
            call random_number_wrapper(p)
            if (nprocz==2) then
              if (lfirst_proc_z) fp(k,izp)=-abs(zp0*sqrt(-2*alog(r))*cos(2*pi*p))
              if (llast_proc_z) fp(k,izp)=abs(zp0*sqrt(-2*alog(r))*cos(2*pi*p))
            else
              fp(k,izp)= zp0*sqrt(-2*alog(r))*cos(2*pi*p)
            endif
          enddo
!
        case ('gaussian-r')
          if (lroot) print*, 'init_particles: Gaussian particle positions'
          do k=1,npar_loc
            call random_number_wrapper(r)
            call random_number_wrapper(p)
            call random_number_wrapper(q)
            fp(k,ixp)= xp0*sqrt(-2*alog(r))*cos(2*pi*p)*cos(2*pi*q)
            fp(k,iyp)= yp0*sqrt(-2*alog(r))*cos(2*pi*p)*sin(2*pi*q)
          enddo
!
        case default
          call fatal_error('init_particles','no such initxxp: '//trim(initxxp(j)))
!
        endselect
!
      enddo !  do j=1,ninit
!
!  Interface for user's own initial condition for position
!
      if (linitial_condition) call initial_condition_xxp(f,fp)
!
!  Particles are not allowed to be present in non-existing dimensions.
!  This would give huge problems with interpolation later.
!
      if (nxgrid==1) fp(1:npar_loc,ixp)=x(nghost+1)
      if (nygrid==1) fp(1:npar_loc,iyp)=y(nghost+1)
      if (nzgrid==1) fp(1:npar_loc,izp)=z(nghost+1)
!
!  Redistribute particles among processors (now that positions are determined).
!
      call boundconds_particles(fp,ipar)
!
!  Map particle positions on the grid. Doing this here is necessary for
!  particle velocity initial conditions that depend on the local particle
!  density, but it breaks conformity with not using particle blocks. To get
!  the same random numbers as when not using particle blocks, one should
!  not sort the particles until the end of the subroutine.
!
      if (learly_particle_map) then
        call map_nearest_grid(fp,ineargrid)
        call sort_particles_iblock(fp,ineargrid,ipar)
        call map_xxp_grid(f,fp,ineargrid)
        if (any(initvvp=='random') .or. &     ! Over reaction, but can
            any(initvvp=='random-x') .or. &   ! cause a lot of frustration
            any(initvvp=='random-y') .or. &   ! when comparing block results
            any(initvvp=='random-z')) then    ! with normal results
          call fatal_error('init_particles','for random particle velocity initial conditions'// &
                           'use learly_particle_map=F')
        endif
      endif
!
!  Initial particle velocity.
!
      do j=1,ninit
!
        select case (initvvp(j))
!
        case ('nothing')
          if (lroot.and.j==1) print*, 'init_particles: No particle velocity set'
!
        case ('zero')
          if (lroot) print*, 'init_particles: Zero particle velocity'
          fp(1:npar_loc,ivpx:ivpz)=0.0
!
        case ('zero-shear')
          if (lroot) print*, 'init_particles: Zero particle velocity'
          fp(1:npar_loc,ivpy)=-Sshear*fp(1:npar_loc,ixp)
          fp(1:npar_loc,ivpx)=0.0
          fp(1:npar_loc,ivpz)=0.0
!
        case ('constant')
          if (lroot) print*, 'init_particles: Constant particle velocity'
          if (lroot) print*, 'init_particles: vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
          if (lcylindrical_coords) then
            fp(1:npar_loc,ivpx)=vpx0*cos(fp(k,iyp))+vpy0*sin(fp(k,iyp))
            fp(1:npar_loc,ivpy)=vpy0*cos(fp(k,iyp))-vpx0*sin(fp(k,iyp))
            fp(1:npar_loc,ivpz)=vpz0
          else
            fp(1:npar_loc,ivpx)=vpx0
            fp(1:npar_loc,ivpy)=vpy0
            fp(1:npar_loc,ivpz)=vpz0
          endif
!
        case ('constant-1')
          if (lroot) print*, 'init_particles: Particle 1 velocity vx,vy,vz=', vpx1, vpy1, vpz1
          do k=1,npar_loc
            if (ipar(k)==1) then
              fp(k,ivpx)=vpx1
              fp(k,ivpy)=vpy1
              fp(k,ivpz)=vpz1
            endif
          enddo
!
        case ('constant-2')
          if (lroot) print*, 'init_particles: Particle 2 velocity vx,vy,vz=', vpx2, vpy2, vpz2
          do k=1,npar_loc
            if (ipar(k)==2) then
              fp(k,ivpx)=vpx2
              fp(k,ivpy)=vpy2
              fp(k,ivpz)=vpz2
            endif
          enddo
!
        case ('constant-3')
          if (lroot) print*, 'init_particles: Particle 3 velocity vx,vy,vz=', vpx3, vpy3, vpz3
          do k=1,npar_loc
            if (ipar(k)==3) then
              fp(k,ivpx)=vpx3
              fp(k,ivpy)=vpy3
              fp(k,ivpz)=vpz3
            endif
          enddo
!
        case ('sinwave-phase')
          if (lroot) print*, 'init_particles: sinwave-phase'
          if (lroot) print*, 'init_particles: vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
          do k=1,npar_loc
            fp(k,ivpx)=fp(k,ivpx)+vpx0*sin(kx_vpx*fp(k,ixp)+ky_vpx*fp(k,iyp)+kz_vpx*fp(k,izp)+phase_vpx)
            fp(k,ivpy)=fp(k,ivpy)+vpy0*sin(kx_vpy*fp(k,ixp)+ky_vpy*fp(k,iyp)+kz_vpy*fp(k,izp)+phase_vpy)
            fp(k,ivpz)=fp(k,ivpz)+vpz0*sin(kx_vpz*fp(k,ixp)+ky_vpz*fp(k,iyp)+kz_vpz*fp(k,izp)+phase_vpz)
          enddo
!
        case ('coswave-phase')
          if (lroot) print*, 'init_particles: coswave-phase'
          if (lroot) print*, 'init_particles: vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
          do k=1,npar_loc
            fp(k,ivpx)=fp(k,ivpx)+vpx0*cos(kx_vpx*fp(k,ixp)+ky_vpx*fp(k,iyp)+kz_vpx*fp(k,izp)+phase_vpx)
            fp(k,ivpy)=fp(k,ivpy)+vpy0*cos(kx_vpy*fp(k,ixp)+ky_vpy*fp(k,iyp)+kz_vpy*fp(k,izp)+phase_vpy)
            fp(k,ivpz)=fp(k,ivpz)+vpz0*cos(kx_vpz*fp(k,ixp)+ky_vpz*fp(k,iyp)+kz_vpz*fp(k,izp)+phase_vpz)
          enddo
!
        case ('random')
          if (lroot) print*, 'init_particles: Random particle velocities; delta_vp0=', delta_vp0
          do k=1,npar_loc
            call random_number_wrapper(r)
            fp(k,ivpx) = fp(k,ivpx) + delta_vp0*(2*r-1)
            call random_number_wrapper(r)
            fp(k,ivpy) = fp(k,ivpy) + delta_vp0*(2*r-1)
            call random_number_wrapper(r)
            fp(k,ivpz) = fp(k,ivpz) + delta_vp0*(2*r-1)
          enddo
!
        case ('random-x')
          if (lroot) print*, 'init_particles: Random particle x-velocity; delta_vp0=', delta_vp0
          do k=1,npar_loc
            call random_number_wrapper(r)
            fp(k,ivpx) = fp(k,ivpx) + delta_vp0*(2*r-1)
          enddo
!
        case ('random-y')
          if (lroot) print*, 'init_particles: Random particle y-velocity; delta_vp0=', delta_vp0
          do k=1,npar_loc
            call random_number_wrapper(r)
            fp(k,ivpy) = fp(k,ivpy) + delta_vp0*(2*r-1)
          enddo
!
        case ('random-z')
          if (lroot) print*, 'init_particles: Random particle z-velocity; delta_vp0=', delta_vp0
          do k=1,npar_loc
            call random_number_wrapper(r)
            fp(k,ivpz) = fp(k,ivpz) + delta_vp0*(2*r-1)
          enddo
!
        case ('average-to-zero')
          call mpireduce_sum(sum(fp(1:npar_loc,ivpx)),vpx_sum)
          call mpireduce_sum(sum(fp(1:npar_loc,ivpy)),vpy_sum)
          call mpireduce_sum(sum(fp(1:npar_loc,ivpz)),vpz_sum)
          call mpibcast_real(vpx_sum)
          call mpibcast_real(vpy_sum)
          call mpibcast_real(vpz_sum)
          fp(1:npar_loc,ivpx)=fp(1:npar_loc,ivpx)-vpx_sum/npar
          fp(1:npar_loc,ivpy)=fp(1:npar_loc,ivpy)-vpy_sum/npar
          fp(1:npar_loc,ivpz)=fp(1:npar_loc,ivpz)-vpz_sum/npar
!
        case ('jeans-wave-dustpar-x')
        ! assumes rhs_poisson_const=1 !
          do k=1,npar_loc
            fp(k,ivpx) = fp(k,ivpx) - amplxxp*(sqrt(1+4*1.0*1.0*tausp**2)-1)/ &
                (2*kx_xxp*1.0*tausp)*sin(kx_xxp*(fp(k,ixp)))
          enddo
!
        case ('dragforce_equilibrium')
!
!  Equilibrium between drag forces on dust and gas and other forces
!  (from Nakagawa, Sekiya, & Hayashi 1986).
!
          if (lroot) then
            print*, 'init_particles: drag equilibrium'
            print*, 'init_particles: beta_glnrho_global=', beta_glnrho_global
          endif
          cs=sqrt(cs20)
          if (.not.learly_particle_map) call fatal_error('init_particles', &
              'must have learly_particle_map=T for dragforce_equilibrium')
          if (ldensity_nolog) then
            call fill_blocks_with_bricks(f,fb,mfarray,irho)
          else
            call fill_blocks_with_bricks(f,fb,mfarray,ilnrho)
          endif
!
          eps = 0.0
          if (ldragforce_equi_global_eps) eps = eps_dtog
!
          if (lroot) print*, 'init_particles: average dust-to-gas ratio=', eps
!  Set gas velocity field.
          do l=l1,l2; do m=m1,m2; do n=n1,n2
!  Take either global or local dust-to-gas ratio.
            if (.not. ldragforce_equi_global_eps .and. ldragforce_gas_par) eps = f(l,m,n,irhop) / get_gas_density(f, l, m, n)
!
            f(l,m,n,iux) = f(l,m,n,iux) - beta_glnrho_global(1)*eps*Omega*tausp/ &
                           ((1.0+eps)**2+(Omega*tausp)**2)*cs
            f(l,m,n,iuy) = f(l,m,n,iuy) + beta_glnrho_global(1)*(1+eps+(Omega*tausp)**2)/ &
                           (2*((1.0+eps)**2+(Omega*tausp)**2))*cs
!
          enddo; enddo; enddo
!  Set particle velocity field.
          do k=1,npar_loc
!  Take either global or local dust-to-gas ratio.
            if (.not. ldragforce_equi_global_eps .and. ldragforce_gas_par) then
              ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
              ib=inearblock(k)
              eps = fb(ix0,iy0,iz0,irhop,ib) / get_gas_density_strat(fb,ix0,iy0,iz0,ib)
            endif
!
            fp(k,ivpx) = fp(k,ivpx) + beta_glnrho_global(1)*Omega*tausp/ &
                         ((1.0+eps)**2+(Omega*tausp)**2)*cs
            fp(k,ivpy) = fp(k,ivpy) + beta_glnrho_global(1)*(1+eps)/ &
                         (2*((1.0+eps)**2+(Omega*tausp)**2))*cs
          enddo
!
        case ('dragforce_equi_dust')
!
!  Equilibrium between drag force and Coriolis force on the dust.
!
          if (lroot) then
            print*, 'init_particles: drag equilibrium dust'
            print*, 'init_particles: beta_dPdr_dust=', beta_dPdr_dust
          endif
!  Set particle velocity field.
          cs=sqrt(cs20)
          do k=1,npar_loc
            fp(k,ivpx) = fp(k,ivpx) + 1/gamma*beta_dPdr_dust/(Omega*tausp+1/(Omega*tausp))*cs
            fp(k,ivpy) = fp(k,ivpy) - 1/gamma*beta_dPdr_dust*Omega*tausp*0.5/ &
                (Omega*tausp+1/(Omega*tausp))*cs
          enddo
!
       case ('Keplerian','keplerian')
!
!  Keplerian velocity based on gravr.
!
          if (lroot) then
            print*, 'init_particles: Keplerian velocity'
            if (lshear) call fatal_error("init_particles", &
                 "Keplerian initial condition is for global disks, not shearing boxes")
          endif
          do k=1,npar_loc
            if (lcartesian_coords) then
              rad=sqrt(fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2)
              OO=sqrt(gravr)*rad**(-1.5)
              fp(k,ivpx) = -OO*fp(k,iyp)
              fp(k,ivpy) =  OO*fp(k,ixp)
              fp(k,ivpz) =  0.0
            elseif (lcylindrical_coords) then
              rad=fp(k,ixp)
              OO=sqrt(gravr)*rad**(-1.5)
              fp(k,ivpx) =  0.0
              fp(k,ivpy) =  OO*rad
              fp(k,ivpz) =  0.0
            elseif (lspherical_coords) then
              rad=fp(k,ixp)*sin(fp(k,iyp))
              OO=sqrt(gravr)*rad**(-1.5)
              fp(k,ivpx) =  0.0
              fp(k,ivpy) =  0.0
              fp(k,ivpz) =  OO*rad
            endif
          enddo
!
!  Explosion.
!
       case ('explosion')
         do k=1,npar_loc
           rad=sqrt(fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2)
           fp(k,ivpx) = delta_vp0*fp(k,ixp)/rp_ext
           fp(k,ivpy) = delta_vp0*fp(k,iyp)/rp_ext
           fp(k,ivpz) = delta_vp0*fp(k,izp)/rp_ext
         enddo
!
        case default
          call fatal_error('init_particles','no such initvvp: '//trim(initvvp(j)))
!
        endselect
!
      enddo ! do j=1,ninit
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_vvp(f,fp)
!
!  Map particle velocity on the grid.
!
      call map_vvp_grid(f,fp,ineargrid)
!
!  Sort particles and map them on the grid. This is where it is done when
!  not using particle blocks, so better use learly_particle_map=F if trying
!  to emulate those results.
!
      if (.not. learly_particle_map) then
        call map_nearest_grid(fp,ineargrid)
        call sort_particles_iblock(fp,ineargrid,ipar)
        call map_xxp_grid(f,fp,ineargrid)
      endif
!
!  Distribute particles in blocks so that each CPU has approximately
!  the same number of particles.
!
      call load_balance_particles(f,fp,ipar)
      call map_nearest_grid(fp,ineargrid)
      call sort_particles_iblock(fp,ineargrid,ipar)
      call map_xxp_grid(f,fp,ineargrid)
!
    endsubroutine init_particles
!***********************************************************************
    subroutine insert_lost_particles(f,fp,ineargrid)
!
!  14-oct-12/dhruba: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (inout) :: fp,ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine insert_lost_particles
!***********************************************************************
    subroutine insert_particles(f,fp,ineargrid)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine insert_particles
!***********************************************************************
    subroutine particles_dragforce_stiff(f,fp,ineargrid)
!
!  10-june-11/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine particles_dragforce_stiff
!***********************************************************************
    subroutine pencil_criteria_particles()
!
!  All pencils that the Particles module depends on are specified here.
!
!  20-04-06/anders: coded
!
      if (idiag_npm/=0 .or. idiag_np2m/=0 .or. idiag_npmax/=0 .or. &
          idiag_npmin/=0 .or. idiag_npmx/=0 .or. idiag_npmy/=0 .or. &
          idiag_npmz/=0) lpenc_diagnos(i_np)=.true.
      if (idiag_rhopm/=0 .or. idiag_rhoprms/=0 .or. idiag_rhop2m/=0 .or. &
          idiag_rhopmax/=0 .or. idiag_rhopmin/=0 .or. idiag_rhopmphi/=0 .or. &
          idiag_rhopmx/=0 .or. idiag_rhopmy/=0 .or. idiag_rhopmz/=0) &
          lpenc_diagnos(i_rhop)=.true.
      if (idiag_dedragp/=0 .or. idiag_decollp/=0) then
        lpenc_diagnos(i_TT1)=.true.
        lpenc_diagnos(i_rho1)=.true.
      endif
      if (idiag_epspmx/=0 .or. idiag_epspmy/=0 .or. idiag_epspmz/=0 .or. &
          idiag_epspmin/=0 .or. idiag_epspmax/=0) &
          lpenc_diagnos(i_epsp)=.true.
      if (idiag_rpvpxmx /= 0 .or. idiag_rpvpymx /= 0 .or. idiag_rpvpzmx /= 0 .or. &
          idiag_rpvpx2mx /= 0 .or. idiag_rpvpy2mx /= 0 .or. idiag_rpvpz2mx /= 0) then
        lpenc_diagnos(i_rhop) = .true.
        lpenc_diagnos(i_uup) = .true.
      endif
      if (idiag_rhopmxz/=0 .or. idiag_rhopmxy/=0 .or. idiag_rhopmr/=0 .or. &
          idiag_sigmap/=0) lpenc_diagnos2d(i_rhop)=.true.
!
    endsubroutine pencil_criteria_particles
!***********************************************************************
    subroutine pencil_interdep_particles(lpencil_in)
!
!  Interdependency among pencils provided by the Particles module
!  is specified here.
!
!  16-feb-06/anders: dummy
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_rhop) .and. irhop==0) lpencil_in(i_np)=.true.
!
      if (lpencil_in(i_uup) .and. iuup == 0) &
          call fatal_error("pencil_interdep_particles", "p%uup is requested but not calculated")
!
      if (lpencil_in(i_epsp)) then
        lpencil_in(i_rhop)=.true.
        lpencil_in(i_rho1)=.true.
      endif
!
      if (lpencil_in(i_grhop)) then
        if (irhop /= 0) then
          if (nprocx /= 1.and.(.not.lcommunicate_rhop)) &
            call fatal_error("pencil_interdep_particles", &
                             "Switch on lcommunicate_rhop=T in particles_run_pars")
        else
          if (nprocx /= 1.and.(.not.lcommunicate_np)) &
            call fatal_error("pencil_interdep_particles", &
                             "Switch on lcommunicate_np=T in particles_run_pars")
        endif
      endif
!
    endsubroutine pencil_interdep_particles
!***********************************************************************
    subroutine calc_pencils_particles(f,p)
!
!  Calculate particle pencils.
!
!  16-feb-06/anders: dummy
!
      use Sub, only: grad
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      if (lpencil(i_np)) p%np=f(l1:l2,m,n,inp)
!
      if (lpencil(i_rhop)) then
        if (irhop/=0) then
          p%rhop=f(l1:l2,m,n,irhop)
        else
          p%rhop=rhop_swarm*f(l1:l2,m,n,inp)
        endif
      endif
!
      if (lpencil(i_grhop)) then
        if (irhop/=0) then
          call grad(f,irhop,p%grhop)
        else
          call grad(f,inp,p%grhop)
          p%grhop=rhop_swarm*p%grhop
        endif
      endif
!
      if (lpencil(i_epsp)) p%epsp=p%rhop*p%rho1
!
      if (lpencil(i_uup) .and. iuup > 0) p%uup = f(l1:l2,m,n,iupx:iupz)
!
    endsubroutine calc_pencils_particles
!***********************************************************************
    subroutine dxxp_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of dust particle position.
!
!  02-jan-05/anders: coded
!
      use General, only: random_number_wrapper, random_seed_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      logical :: lheader, lfirstcall=.true.
!
      intent (in) :: f, fp, ineargrid
      intent (inout) :: df, dfp
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
!
!  Identify module and boundary conditions.
!
      if (lheader) print*,'dxxp_dt: Calculate dxxp_dt'
      if (lheader) then
        print*, 'dxxp_dt: Particles boundary condition bcpx=', bcpx
        print*, 'dxxp_dt: Particles boundary condition bcpy=', bcpy
        print*, 'dxxp_dt: Particles boundary condition bcpz=', bcpz
        print*, 'dxxp_dt: Set rate of change of particle position equal to particle velocity.'
      endif
!
!  The rate of change of a particle's position is the particle's velocity.
!  If pointmasses are used do the evolution in that module, for better
!  conservation of the Jacobi constant.
!
      if (.not.lpointmasses) then
        if (lcartesian_coords) then
!
          if (nxgrid/=1) dfp(1:npar_loc,ixp) = dfp(1:npar_loc,ixp) + fp(1:npar_loc,ivpx)
          if (nygrid/=1) dfp(1:npar_loc,iyp) = dfp(1:npar_loc,iyp) + fp(1:npar_loc,ivpy)
          if (nzgrid/=1) dfp(1:npar_loc,izp) = dfp(1:npar_loc,izp) + fp(1:npar_loc,ivpz)
!
        elseif (lcylindrical_coords) then
!
          if (nxgrid/=1) dfp(1:npar_loc,ixp) = dfp(1:npar_loc,ixp) + fp(1:npar_loc,ivpx)
          if (nygrid/=1) &
               dfp(1:npar_loc,iyp) = dfp(1:npar_loc,iyp) + &
               fp(1:npar_loc,ivpy)/max(fp(1:npar_loc,ixp),tini)
          if (nzgrid/=1) dfp(1:npar_loc,izp) = dfp(1:npar_loc,izp) + fp(1:npar_loc,ivpz)
!
        elseif (lspherical_coords) then
!
          if (nxgrid/=1) dfp(1:npar_loc,ixp) = dfp(1:npar_loc,ixp) + fp(1:npar_loc,ivpx)
          if (nygrid/=1) &
               dfp(1:npar_loc,iyp) = dfp(1:npar_loc,iyp) + &
               fp(1:npar_loc,ivpy)/max(fp(1:npar_loc,ixp),tini)
          if (nzgrid/=1) &
               dfp(1:npar_loc,izp) = dfp(1:npar_loc,izp) + &
               fp(1:npar_loc,ivpz)/(max(fp(1:npar_loc,ixp),tini)*&
                                                      sin(fp(1:npar_loc,iyp)))
!
        endif
      endif
!
!  With shear there is an extra term due to the background shear flow.
!
      if (lshear.and.nygrid/=1) dfp(1:npar_loc,iyp) = &
          dfp(1:npar_loc,iyp) - qshear*Omega*fp(1:npar_loc,ixp)
!
      if (lfirstcall) lfirstcall=.false.
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dxxp_dt
!***********************************************************************
    subroutine dvvp_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of dust particle velocity.
!
!  29-dec-04/anders: coded
!
      use Diagnostics
      use EquationOfState, only: cs20
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension(3) :: ggp
      real :: Omega2, rr, vv, OO2
      integer :: k, npar_found
      logical :: lheader, lfirstcall=.true.
!
      intent (in) :: f, fp, ineargrid
      intent (inout) :: df, dfp
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
      if (lheader) print*,'dvvp_dt: Calculate dvvp_dt'
!
!  Add Coriolis force from rotating coordinate frame.
!
      if (Omega/=0.) then
        if (lcoriolis_force_par) then
          if (lheader) print*,'dvvp_dt: Add Coriolis force; Omega=', Omega
          Omega2=2*Omega
          if (.not.lspherical_coords) then
            dfp(1:npar_loc,ivpx) = dfp(1:npar_loc,ivpx) + Omega2*fp(1:npar_loc,ivpy)
            dfp(1:npar_loc,ivpy) = dfp(1:npar_loc,ivpy) - Omega2*fp(1:npar_loc,ivpx)
          else
            call not_implemented('dvvp_dt','Coriolis force on the particles for spherical coordinates')
          endif
        endif
!
! Centrifugal force
!
        if (lcentrifugal_force_par) then
          if (lheader) print*,'dvvp_dt: Add Centrifugal force; Omega=', Omega
          if (lcartesian_coords) then
!
            dfp(1:npar_loc,ivpx) = dfp(1:npar_loc,ivpx) + Omega**2*fp(1:npar_loc,ixp)
!
            dfp(1:npar_loc,ivpy) = dfp(1:npar_loc,ivpy) + Omega**2*fp(1:npar_loc,iyp)
!
          elseif (lcylindrical_coords) then
            dfp(1:npar_loc,ivpx) = dfp(1:npar_loc,ivpx) + Omega**2*fp(1:npar_loc,ixp)
          else
            call not_implemented('dvvp_dt','Centrifugal force on the particles for spherical coordinates')
          endif
        endif
!
!  With shear there is an extra term due to the background shear flow.
!
        if (lshear) dfp(1:npar_loc,ivpy) = dfp(1:npar_loc,ivpy) + qshear*Omega*fp(1:npar_loc,ivpx)
      endif
!
!  Add constant background pressure gradient beta=alpha*H0/r0, where alpha
!  comes from a global pressure gradient P = P0*(r/r0)^alpha.
!  (the term must be added to the dust equation of motion when measuring
!  velocities relative to the shear flow modified by the global pressure grad.)
!
      if (beta_dPdr_dust/=0.0 .and. t>=tstart_dragforce_par) &
        dfp(1:npar_loc,ivpx) = dfp(1:npar_loc,ivpx) + 1/gamma*cs20*beta_dPdr_dust_scaled
!
!  Gravity on the particles.
!
      if (t>=tstart_grav_par) then
!
!  Gravity in the x-direction.
!
        select case (gravx_profile)
!
          case ('zero','')
            if (lheader) print*, 'dvvp_dt: No gravity in x-direction.'
!
          case ('linear')
            if (lheader) print*, 'dvvp_dt: Linear gravity field in x-direction.'
            dfp(1:npar_loc,ivpx)=dfp(1:npar_loc,ivpx) - nu_epicycle2*fp(1:npar_loc,ixp)
!
          case ('plain')
            if (lheader) print*, 'dvvp_dt: Plain gravity field in x-direction.'
            dfp(1:npar_loc,ivpx)=dfp(1:npar_loc,ivpx) - gravx
!
          case ('sinusoidal')
            if (lheader) &
                print*, 'dvvp_dt: Sinusoidal gravity field in x-direction.'
            dfp(1:npar_loc,ivpx)=dfp(1:npar_loc,ivpx) - gravx*sin(kx_gg*fp(1:npar_loc,ixp))
!
          case ('z2')
            if (lheader) print *, 'dvvp_dt: g_x = gravx * z^2, gravx = ', gravx
            dfp(1:npar_loc,ivpx) = dfp(1:npar_loc,ivpx) + gravx * fp(1:npar_loc,izp)**2
!
          case default
            call fatal_error('dvvp_dt','no such gravx_profile: '//trim(gravx_profile))
!
        endselect
!
!  Gravity in the z-direction.
!
        select case (gravz_profile)
!
          case ('zero','')
            if (lheader) print*, 'dvvp_dt: No gravity in z-direction.'
!
          case ('const','plain')
            if (lheader) print*, 'dvvp_dt: Constant gravity field in z-direction.'
            dfp(1:npar_loc,ivpz)=dfp(1:npar_loc,ivpz) + gravz
!
          case ('linear')
            if (lheader) print*, 'dvvp_dt: Linear gravity field in z-direction.'
            dfp(1:npar_loc,ivpz)=dfp(1:npar_loc,ivpz) - nu_epicycle2*fp(1:npar_loc,izp)
!
          case ('sinusoidal')
            if (lheader) print*, 'dvvp_dt: Sinusoidal gravity field in z-direction.'
            dfp(1:npar_loc,ivpz)=dfp(1:npar_loc,ivpz) - gravz*sin(kz_gg*fp(1:npar_loc,izp))
!
          case default
            call fatal_error('dvvp_dt','no such gravz_profile: '//trim(gravz_profile))
!
        endselect
!
!  Radial gravity.
!
        select case (gravr_profile)
!
        case ('zero','')
          if (lheader) print*, 'dvvp_dt: No radial gravity'
!
        case ('newtonian-central','newtonian')
          if (lpointmasses) call fatal_error('dvvp_dt','You are using massive particles. '// &
              'The N-body code should take care of the stellar-like '// &
              'gravity on the dust. Switch off the gravr_profile="newtonian" on particles_init')
          if (lheader) print*, 'dvvp_dt: Newtonian gravity from a fixed central object'
          do k=1,npar_loc
            if (lcartesian_coords) then
              if (lcylindrical_gravity_par) then
                rr=sqrt(fp(k,ixp)**2+fp(k,iyp)**2+gravsmooth2)
              else
                rr=sqrt(fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2+gravsmooth2)
              endif
              OO2=rr**(-3.)*gravr
              ggp(1) = -fp(k,ixp)*OO2
              ggp(2) = -fp(k,iyp)*OO2
              if (lcylindrical_gravity_par) then
                ggp(3) = 0.
              else
                ggp(3) = -fp(k,izp)*OO2
              endif
              dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + ggp
            elseif (lcylindrical_coords) then
              if (lcylindrical_gravity_par) then
                rr=sqrt(fp(k,ixp)**2+gravsmooth2)
              else
                rr=sqrt(fp(k,ixp)**2+fp(k,izp)**2+gravsmooth2)
              endif
              OO2=rr**(-3.)*gravr
              ggp(1) = -fp(k,ixp)*OO2
              ggp(2) = 0.0
              if (lcylindrical_gravity_par) then
                ggp(3) = 0.
              else
                ggp(3) = -fp(k,izp)*OO2
              endif
              dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + ggp
            elseif (lspherical_coords) then
              rr=sqrt(fp(k,ixp)**2+gravsmooth2)
              OO2=rr**(-3)*gravr
              ggp(1) = -fp(k,ixp)*OO2
              ggp(2) = 0.0; ggp(3) = 0.0
              if (lcylindrical_gravity_par) call fatal_error("dvvp_dt",&
                   "No cylindrical gravity in spherical coordinates.")
              dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + ggp
            endif
!  Limit time-step if particles close to gravity source.
            if (ldt_grav_par.and.(lfirst.and.ldt)) then
              if (lcartesian_coords) then
                vv=sqrt(fp(k,ivpx)**2+fp(k,ivpy)**2+fp(k,ivpz)**2)
              elseif (lcylindrical_coords) then
                vv=sqrt(fp(k,ivpx)**2+fp(k,ivpz)**2)
              elseif (lspherical_coords) then
                vv=abs(fp(k,ivpx))
              endif
              dt1_max(ineargrid(k,1)-nghost)=max(dt1_max(ineargrid(k,1)-nghost),vv/rr/cdtpgrav)
            endif
          enddo
!
        case default
          call fatal_error('dvvp_dt','no such gravr_profile: '//trim(gravr_profile))
!
        endselect
!
      endif
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_nparmin/=0) call max_name(-npar_loc,idiag_nparmin,lneg=.true.)
        call max_name(+npar_loc,idiag_nparmax)
        call sum_par_name(fp(1:npar_loc,ixp),idiag_xpm)
        call sum_par_name(fp(1:npar_loc,iyp),idiag_ypm)
        call sum_par_name(fp(1:npar_loc,izp),idiag_zpm)
        if (idiag_xp2m/=0) call sum_par_name(fp(1:npar_loc,ixp)**2,idiag_xp2m)
        if (idiag_yp2m/=0) call sum_par_name(fp(1:npar_loc,iyp)**2,idiag_yp2m)
        if (idiag_zp2m/=0) call sum_par_name(fp(1:npar_loc,izp)**2,idiag_zp2m)
        if (idiag_rpm/=0)  call sum_par_name(sqrt(fp(1:npar_loc,ixp)**2+ &
            fp(1:npar_loc,iyp)**2+fp(1:npar_loc,izp)**2),idiag_rpm)
        if (idiag_rp2m/=0) call sum_par_name(fp(1:npar_loc,ixp)**2+ &
            fp(1:npar_loc,iyp)**2+fp(1:npar_loc,izp)**2,idiag_rp2m)
        call sum_par_name(fp(1:npar_loc,ivpx),idiag_vpxm)
        call sum_par_name(fp(1:npar_loc,ivpy),idiag_vpym)
        call sum_par_name(fp(1:npar_loc,ivpz),idiag_vpzm)
        if (idiag_lpxm/=0) call sum_par_name( &
            fp(1:npar_loc,iyp)*fp(1:npar_loc,ivpz)- &
            fp(1:npar_loc,izp)*fp(1:npar_loc,ivpy),idiag_lpxm)
        if (idiag_lpym/=0) call sum_par_name( &
            fp(1:npar_loc,izp)*fp(1:npar_loc,ivpx)- &
            fp(1:npar_loc,ixp)*fp(1:npar_loc,ivpz),idiag_lpym)
        if (idiag_lpzm/=0) call sum_par_name( &
            fp(1:npar_loc,ixp)*fp(1:npar_loc,ivpy)- &
            fp(1:npar_loc,iyp)*fp(1:npar_loc,ivpx),idiag_lpzm)
        if (idiag_lpx2m/=0) call sum_par_name( &
            (fp(1:npar_loc,iyp)*fp(1:npar_loc,ivpz)- &
            fp(1:npar_loc,izp)*fp(1:npar_loc,ivpy))**2,idiag_lpx2m)
        if (idiag_lpy2m/=0) call sum_par_name( &
            (fp(1:npar_loc,izp)*fp(1:npar_loc,ivpx)- &
            fp(1:npar_loc,ixp)*fp(1:npar_loc,ivpz))**2,idiag_lpy2m)
        if (idiag_lpz2m/=0) call sum_par_name((fp(1:npar_loc,ixp)*fp(1:npar_loc,ivpy)- &
            fp(1:npar_loc,iyp)*fp(1:npar_loc,ivpx))**2,idiag_lpz2m)
        if (idiag_vpx2m/=0) call sum_par_name(fp(1:npar_loc,ivpx)**2,idiag_vpx2m)
        if (idiag_vpy2m/=0) call sum_par_name(fp(1:npar_loc,ivpy)**2,idiag_vpy2m)
        if (idiag_vpz2m/=0) call sum_par_name(fp(1:npar_loc,ivpz)**2,idiag_vpz2m)
        if (idiag_ekinp/=0) then
          if (lparticles_density) then
            call sum_par_name(0.5*fp(1:npar_loc,irhopswarm)* &
                sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2),idiag_ekinp)
          else
            if (lcartesian_coords.and.(all(lequidist))) then
              call sum_par_name(0.5*rhop_swarm*npar_per_cell* &
                   sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2),idiag_ekinp)
            else
              call sum_par_name(0.5*mp_swarm*sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2),idiag_ekinp)
            endif
          endif
        endif
        if (idiag_epotpm/=0) call sum_par_name( &
            -gravr/sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)),idiag_epotpm)
        if (idiag_vpmax/=0) call max_par_name(sqrt(sum(fp(1:npar_loc,ivpx:ivpz)**2,2)),idiag_vpmax)
        call max_par_name(fp(1:npar_loc,ivpx),idiag_vpxmax)
        call max_par_name(fp(1:npar_loc,ivpy),idiag_vpymax)
        call max_par_name(fp(1:npar_loc,ivpz),idiag_vpzmax)
        if (idiag_eccpxm/=0) call sum_par_name( &
            sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2)*fp(1:npar_loc,ixp)/gravr- &
            sum(fp(1:npar_loc,ixp:izp)*fp(1:npar_loc,ivpx:ivpz),dim=2)* &
            fp(1:npar_loc,ivpx)/gravr-fp(1:npar_loc,ixp)/ &
            sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)),idiag_eccpxm)
        if (idiag_eccpym/=0) call sum_par_name( &
            sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2)*fp(1:npar_loc,iyp)/gravr- &
            sum(fp(1:npar_loc,ixp:izp)*fp(1:npar_loc,ivpx:ivpz),dim=2)* &
            fp(1:npar_loc,ivpy)/gravr-fp(1:npar_loc,iyp)/ &
            sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)),idiag_eccpym)
        if (idiag_eccpzm/=0) call sum_par_name( &
            sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2)*fp(1:npar_loc,izp)/gravr- &
            sum(fp(1:npar_loc,ixp:izp)*fp(1:npar_loc,ivpx:ivpz),dim=2)* &
            fp(1:npar_loc,ivpz)/gravr-fp(1:npar_loc,izp)/ &
            sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)),idiag_eccpzm)
        if (idiag_eccpx2m/=0) call sum_par_name(( &
            sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2)*fp(1:npar_loc,ixp)/gravr- &
            sum(fp(1:npar_loc,ixp:izp)*fp(1:npar_loc,ivpx:ivpz),dim=2)* &
            fp(1:npar_loc,ivpx)/gravr-fp(1:npar_loc,ixp)/ &
            sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)))**2,idiag_eccpx2m)
        if (idiag_eccpy2m/=0) call sum_par_name(( &
            sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2)*fp(1:npar_loc,iyp)/gravr- &
            sum(fp(1:npar_loc,ixp:izp)*fp(1:npar_loc,ivpx:ivpz),dim=2)* &
            fp(1:npar_loc,ivpy)/gravr-fp(1:npar_loc,iyp)/ &
            sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)))**2,idiag_eccpy2m)
        if (idiag_eccpz2m/=0) call sum_par_name(( &
            sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2)*fp(1:npar_loc,izp)/gravr- &
            sum(fp(1:npar_loc,ixp:izp)*fp(1:npar_loc,ivpx:ivpz),dim=2)* &
            fp(1:npar_loc,ivpz)/gravr-fp(1:npar_loc,izp)/ &
            sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)))**2,idiag_eccpz2m)
        if (idiag_mpt/=0) then
          if (lparticles_density) then
            call integrate_par_name((/fp(1:npar_loc,irhopswarm)/),idiag_mpt)
          elseif (lparticles_radius.and.lparticles_number) then
            call integrate_par_name((/four_pi_rhopmat_over_three* &
                fp(1:npar_loc,iap)**3*fp(1:npar_loc,inpswarm)/),idiag_mpt)
          endif
        endif
        if (idiag_npargone/=0) then
          call count_particles(ipar,npar_found)
          call save_name(real(npar-npar_found),idiag_npargone)
        endif
        if (idiag_deshearbcsm/=0) call sum_name(energy_gain_shear_bcs/npar,idiag_deshearbcsm)
      endif
!
      if (lfirstcall) lfirstcall=.false.
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dvvp_dt
!***********************************************************************
    subroutine dxxp_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle position (called from main pencil loop).
!
!  25-apr-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dxxp_dt_pencil
!***********************************************************************
    subroutine dvvp_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of dust particle velocity (called from main pencil loop).
!
!  25-apr-06/anders: coded
!
      use Particles_spin, only: calc_liftforce
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)

      call calc_diagnostics_particles(fp,p,ineargrid)

    endsubroutine dvvp_dt_pencil
!***********************************************************************
    subroutine calc_diagnostics_particles(fp,p,ineargrid)

      use Diagnostics

      real, dimension (mpar_loc,mparray) :: fp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
!  Diagnostic output
!
      if (ldiagnos) then
        call sum_mn_name(p%np,idiag_npm)
        if (idiag_np2m/=0) call sum_mn_name(p%np**2,idiag_np2m)
        call max_mn_name(p%np,idiag_npmax)
        if (idiag_npmin/=0)    call max_mn_name(-p%np,idiag_npmin,lneg=.true.)
        call sum_mn_name(p%rhop,idiag_rhopm)
        if (idiag_rhop2m/=0 )  call sum_mn_name(p%rhop**2,idiag_rhop2m)
        if (idiag_rhoprms/=0)  call sum_mn_name(p%rhop**2,idiag_rhoprms,lsqrt=.true.)
        call max_mn_name(p%rhop,idiag_rhopmax)
        if (idiag_rhopmin/=0)  call max_mn_name(-p%rhop,idiag_rhopmin,lneg=.true.)
        call max_mn_name(p%epsp,idiag_epspmax)
        if (idiag_epspmin/=0)  call max_mn_name(-p%epsp,idiag_epspmin,lneg=.true.)
      endif
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1
!
      if (l1davgfirst) then
        call yzsum_mn_name_x(p%np,idiag_npmx)
        call xzsum_mn_name_y(p%np,idiag_npmy)
        call xysum_mn_name_z(p%np,idiag_npmz)
        call yzsum_mn_name_x(p%rhop,idiag_rhopmx)
        call xzsum_mn_name_y(p%rhop,idiag_rhopmy)
        call xysum_mn_name_z(p%rhop,idiag_rhopmz)
        call yzsum_mn_name_x(p%epsp,idiag_epspmx)
        call xzsum_mn_name_y(p%epsp,idiag_epspmy)
        call xysum_mn_name_z(p%epsp,idiag_epspmz)
        call phizsum_mn_name_r(p%rhop,idiag_rhopmr)
!
        if (idiag_rpvpxmx /= 0) call yzsum_mn_name_x(p%rhop * p%uup(:,1), idiag_rpvpxmx)
        if (idiag_rpvpymx /= 0) call yzsum_mn_name_x(p%rhop * p%uup(:,2), idiag_rpvpymx)
        if (idiag_rpvpzmx /= 0) call yzsum_mn_name_x(p%rhop * p%uup(:,3), idiag_rpvpzmx)
        if (idiag_rpvpx2mx /= 0) call yzsum_mn_name_x(p%rhop * p%uup(:,1)**2, idiag_rpvpx2mx)
        if (idiag_rpvpy2mx /= 0) call yzsum_mn_name_x(p%rhop * p%uup(:,2)**2, idiag_rpvpy2mx)
        if (idiag_rpvpz2mx /= 0) call yzsum_mn_name_x(p%rhop * p%uup(:,3)**2, idiag_rpvpz2mx)
      endif
!
      if (l2davgfirst) then
        call phisum_mn_name_rz(p%rhop,idiag_rhopmphi)
        call zsum_mn_name_xy(p%rhop,idiag_rhopmxy)
        call ysum_mn_name_xz(p%rhop,idiag_rhopmxz)
        call zsum_mn_name_xy(p%rhop,idiag_sigmap,lint=.true.)
      endif
!
    endsubroutine calc_diagnostics_particles
!***********************************************************************
    subroutine dxxp_dt_blocks(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle position in blocks.
!
!  25-apr-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      integer :: k, ix0, iy0, iz0, iblock
      real :: dt1_advpx, dt1_advpy, dt1_advpz
!
!  Contribution of dust particles to time step.
!
      do iblock=0,nblock_loc-1
        if (lfirst.and.ldt.and.ldt_adv_par) then
          if (npar_iblock(iblock)/=0) then
            do k=k1_iblock(iblock),k2_iblock(iblock)
              ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
              dt1_advpx=abs(fp(k,ivpx))*dx1b(ix0,iblock)
              if (lshear .and. .not. lshearadvection_as_shift) then
                dt1_advpy=(-qshear*Omega*fp(k,ixp)+abs(fp(k,ivpy)))*dy1b(iy0,iblock)
              else
                dt1_advpy=abs(fp(k,ivpy))*dy1b(iy0,iblock)
              endif
              dt1_advpz=abs(fp(k,ivpz))*dz1b(iz0,iblock)
              dt1_max(ix0+nghost-l1)=max(dt1_max(ix0+nghost-l1), &
                 sqrt(dt1_advpx**2+dt1_advpy**2+dt1_advpz**2)/cdtp)
            enddo
          endif
        endif
      enddo
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(dfp)
!
    endsubroutine dxxp_dt_blocks
!***********************************************************************
    subroutine dvvp_dt_blocks(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle velocity in blocks.
!
!  25-apr-06/anders: coded
!
      use Diagnostics
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (mxb,myb,mzb) :: dt1_drag, dt1_drag_gas, dt1_drag_dust
      real, dimension (3) :: dragforce, uup
      real :: rho1_point, tausp1_par, rhop_swarm_par
      real :: weight, weight_x, weight_y, weight_z
      integer :: k, l, ix0, iy0, iz0, ib, iblock, lb, mb, nb
      integer :: ixx, iyy, izz, ixx0, iyy0, izz0, ixx1, iyy1, izz1
      logical :: lsink
!
      intent (inout) :: f, df, dfp, fp, ineargrid
!
!  Loop over adopted blocks.
!
      do iblock=0,nblock_loc-1
        ib=iblock
!
!  Drag force on particles and on gas.
!
        if (ldragforce_dust_par .and. t>=tstart_dragforce_par) then
          if (headtt) print*, 'dvvp_dt_blocks: Add drag force; tausp=', tausp
!
          if (npar_iblock(iblock)/=0) then
!
            if (lfirst.and.ldt) then
              dt1_drag_dust=0.0
              if (ldragforce_gas_par) dt1_drag_gas=0.0
            endif
!
!  Loop over all particles in considered block.
!
            do k=k1_iblock(iblock),k2_iblock(iblock)
!
!  Exclude sink particles from the drag calculations
!
              lsink=.false.
              if (lparticles_sink) then
                if (fp(k,iaps)>0.0) lsink=.true.
              endif
              if (.not.lsink) then
                ix0=ineargrid(k,1)
                iy0=ineargrid(k,2)
                iz0=ineargrid(k,3)
!
!  The interpolated gas velocity must be calculated here.
!
                if (lhydro) then
                  if (lparticlemesh_cic) then
                    call interpolate_linear(f,iux,iuz,fp(k,ixp:izp),uup,ineargrid(k,:),iblock,ipar(k))
                  elseif (lparticlemesh_tsc) then
                    if (linterpolate_spline) then
                      call interpolate_quadratic_spline(f,iux,iuz, &
                          fp(k,ixp:izp),uup,ineargrid(k,:),iblock,ipar(k))
                    else
                      call interpolate_quadratic(f,iux,iuz, &
                          fp(k,ixp:izp),uup,ineargrid(k,:),iblock,ipar(k))
                    endif
                  else
                    uup=fb(ix0,iy0,iz0,iux:iuz,iblock)
                  endif
                else
                  uup=0.0
                endif
!
!  Get the friction time. For the case of |uup| ~> cs, the Epstein drag law
!  is dependent on the relative mach number, hence the need to feed uup as
!  an optional argument to get_frictiontime.
!
                call get_frictiontime(f,fp,ineargrid,k,tausp1_par,iblock,uup)
!
!  Calculate and add drag force.
!
                dragforce = -tausp1_par*(fp(k,ivpx:ivpz)-uup)
!
                dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + dragforce
!
!  Back-reaction friction force from particles on gas. Three methods are
!  implemented for assigning a particle to the mesh (see Hockney & Eastwood):
!
!    0. NGP (Nearest Grid Point)
!       The entire effect of the particle goes to the nearest grid point.
!    1. CIC (Cloud In Cell)
!       The particle has a region of influence with the size of a grid cell.
!       This is equivalent to a first order (spline) interpolation scheme.
!    2. TSC (Triangular Shaped Cloud)
!       The particle is spread over a length of two grid cells, but with
!       a density that falls linearly outwards.
!       This is equivalent to a second order spline interpolation scheme.
!
                if (ldragforce_gas_par) then
!
!  Cloud In Cell (CIC) scheme.
!
                  if (lparticlemesh_cic) then
                    ixx0=ix0; iyy0=iy0; izz0=iz0
                    ixx1=ix0; iyy1=iy0; izz1=iz0
!
!  Particle influences the 8 surrounding grid points. The reference point is
!  the grid point at the lower left corner.
!
                    if ( (xb(ix0,ib)>fp(k,ixp)) .and. nxgrid/=1) ixx0=ixx0-1
                    if ( (yb(iy0,ib)>fp(k,iyp)) .and. nygrid/=1) iyy0=iyy0-1
                    if ( (zb(iz0,ib)>fp(k,izp)) .and. nzgrid/=1) izz0=izz0-1
                    if (nxgrid/=1) ixx1=ixx0+1
                    if (nygrid/=1) iyy1=iyy0+1
                    if (nzgrid/=1) izz1=izz0+1
                    do izz=izz0,izz1; do iyy=iyy0,iyy1; do ixx=ixx0,ixx1
                      weight=1.0
                      if (nxgrid/=1) weight=weight*( 1.0-abs(fp(k,ixp)-xb(ixx,ib))*dx1b(ixx,ib) )
                      if (nygrid/=1) weight=weight*( 1.0-abs(fp(k,iyp)-yb(iyy,ib))*dy1b(iyy,ib) )
                      if (nzgrid/=1) weight=weight*( 1.0-abs(fp(k,izp)-zb(izz,ib))*dz1b(izz,ib) )
                      rho1_point = 1.0 / get_gas_density_strat(fb,ixx,iyy,izz,ib)
!  Add friction force to grid point.
                      call get_rhopswarm(mp_swarm,fp,k,ixx,iyy,izz,ib,rhop_swarm_par)
                      if (.not.lcompensate_sedimentation) then
                        dfb(ixx,iyy,izz,iux:iuz,ib)=dfb(ixx,iyy,izz,iux:iuz,ib)- &
                           rhop_swarm_par*rho1_point*dragforce*weight
                      else
                        dfb(ixx,iyy,izz,iux:iuz,ib)=dfb(ixx,iyy,izz,iux:iuz,ib)- &
                           rhop_swarm_par*rho1_point*dragforce*weight*compensate_sedimentation
                      endif
                    enddo; enddo; enddo
!
!  Triangular Shaped Cloud (TSC) scheme.
!
                  elseif (lparticlemesh_tsc) then
!
!  Particle influences the 27 surrounding grid points, but has a density that
!  decreases with the distance from the particle centre.
!
                    if (nxgrid/=1) then
                      ixx0=ix0-1; ixx1=ix0+1
                    else
                      ixx0=ix0  ; ixx1=ix0
                    endif
                    if (nygrid/=1) then
                      iyy0=iy0-1; iyy1=iy0+1
                    else
                      iyy0=iy0  ; iyy1=iy0
                    endif
                    if (nzgrid/=1) then
                      izz0=iz0-1; izz1=iz0+1
                    else
                      izz0=iz0  ; izz1=iz0
                    endif
!
!  The nearest grid point is influenced differently than the left and right
!  neighbours are. A particle that is situated exactly on a grid point gives
!  3/4 contribution to that grid point and 1/8 to each of the neighbours.
!
                    do izz=izz0,izz1; do iyy=iyy0,iyy1; do ixx=ixx0,ixx1
                      if ( ((ixx-ix0)==-1) .or. ((ixx-ix0)==+1) ) then
                        weight_x=1.125-1.5* abs(fp(k,ixp)-xb(ixx,ib))*dx1b(ixx,ib) + &
                                       0.5*(abs(fp(k,ixp)-xb(ixx,ib))*dx1b(ixx,ib))**2
                      else
                        if (nxgrid/=1) weight_x=0.75-((fp(k,ixp)-xb(ixx,ib))*dx1b(ixx,ib))**2
                      endif
                      if ( ((iyy-iy0)==-1) .or. ((iyy-iy0)==+1) ) then
                        weight_y=1.125-1.5* abs(fp(k,iyp)-yb(iyy,ib))*dy1b(iyy,ib) + &
                                       0.5*(abs(fp(k,iyp)-yb(iyy,ib))*dy1b(iyy,ib))**2
                      else
                        if (nygrid/=1) weight_y=0.75-((fp(k,iyp)-yb(iyy,ib))*dy1b(iyy,ib))**2
                      endif
                      if ( ((izz-iz0)==-1) .or. ((izz-iz0)==+1) ) then
                        weight_z=1.125-1.5* abs(fp(k,izp)-zb(izz,ib))*dz1b(izz,ib) + &
                                       0.5*(abs(fp(k,izp)-zb(izz,ib))*dz1b(izz,ib))**2
                      else
                        if (nzgrid/=1) weight_z=0.75-((fp(k,izp)-zb(izz,ib))*dz1b(izz,ib))**2
                      endif
!
                      weight=1.0
!
                      if (nxgrid/=1) weight=weight*weight_x
                      if (nygrid/=1) weight=weight*weight_y
                      if (nzgrid/=1) weight=weight*weight_z
                      rho1_point = 1.0 / get_gas_density_strat(fb,ixx,iyy,izz,ib)
!  Add friction force to grid point.
                      call get_rhopswarm(mp_swarm,fp,k,ixx,iyy,izz,ib,rhop_swarm_par)
                      if (.not.lcompensate_sedimentation) then
                        dfb(ixx,iyy,izz,iux:iuz,ib)=dfb(ixx,iyy,izz,iux:iuz,ib) - &
                          rhop_swarm_par*rho1_point*dragforce*weight
                      else
                        dfb(ixx,iyy,izz,iux:iuz,ib)=dfb(ixx,iyy,izz,iux:iuz,ib) - &
                          rhop_swarm_par*rho1_point*dragforce*weight*compensate_sedimentation
                      endif
                   enddo; enddo; enddo
                  else
!
!  Nearest Grid Point (NGP) scheme.
!
                    rho1_point = 1.0 / get_gas_density_strat(fb,ix0,iy0,iz0,ib)
                    !WL: Why is this l being defined?
                    l=ineargrid(k,1)
                    call get_rhopswarm(mp_swarm,fp,k,ix0,iy0,iz0,ib, &
                        rhop_swarm_par)
                    if (.not.lcompensate_sedimentation) then
                      dfb(ix0,iy0,iz0,iux:iuz,ib) = dfb(ix0,iy0,iz0,iux:iuz,ib) - &
                        rhop_swarm_par*rho1_point*dragforce
                    else
                      dfb(ix0,iy0,iz0,iux:iuz,ib) = dfb(ix0,iy0,iz0,iux:iuz,ib) - &
                        rhop_swarm_par*rho1_point*dragforce*compensate_sedimentation
                    endif
                  endif
                endif
!
!  The minimum friction time of particles in a grid cell sets the local friction
!  time-step when there is only drag force on the dust,
!    dt1_drag = max(1/tausp)
!
!  With drag force on the gas as well, the maximum time-step is set as
!    dt1_drag = Sum_k[eps_k/tau_k]
!
                if (lfirst.and.ldt) then
                  dt1_drag_dust(ix0,iy0,iz0)= &
                       max(dt1_drag_dust(ix0,iy0,iz0),tausp1_par)
                  if (ldragforce_gas_par) then
                    rho1_point = 1.0 / get_gas_density_strat(fb,ix0,iy0,iz0,ib)
                    if (fb(ix0,iy0,iz0,inp,ib)/=0.0) &
                        dt1_drag_gas(ix0,iy0,iz0)=dt1_drag_gas(ix0,iy0,iz0)+ &
                        rhop_swarm_par*rho1_point*tausp1_par
                  endif
                endif
              endif
            enddo
          endif
!
!  Contribution of friction force to time-step. Dust and gas inverse friction
!  time-steps are added up to give a valid expression even when the two are
!  of similar magnitude.
!
          if (lfirst.and.ldt) then
            if (ldragforce_gas_par) then
              dt1_drag=dt1_drag_dust+dt1_drag_gas
            else
              dt1_drag=dt1_drag_dust
            endif
            dt1_drag=dt1_drag/cdtp
            do nb=n1b,n2b; do mb=m1b,m2b; do lb=l1b,l2b
              if (fb(lb,mb,nb,inp,iblock)/=0.0) &
                dt1_max(lb-nghostb)=max(dt1_max(lb-nghostb),dt1_drag(lb,mb,nb))
            enddo; enddo; enddo
          endif
        else
!
!  No particles in this block.
!
          if (lfirst.and.ldt) dt1_drag=0.0
        endif
      enddo
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_nparbmax/=0) call max_name(maxval(npar_iblock(0:nblock_loc-1)),idiag_nparbmax)
        if (idiag_nblockmin/=0) call max_name(-nblock_loc,idiag_nblockmin,lneg=.true.)
        call max_name(nblock_loc,idiag_nblockmax)
        if (idiag_dtdragp /= 0) call max_name(maxval(dt1_drag), idiag_dtdragp, l_dt=.true.)
      endif
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
!
    endsubroutine dvvp_dt_blocks
!***********************************************************************
    subroutine remove_particles_sink_simple(f,fp,dfp,ineargrid)
!
!  Subroutine for taking particles out of the simulation due to their proximity
!  to a sink particle or sink point.
!
!  25-sep-08/anders: coded
!
      use Solid_Cells
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      real :: rp
      integer :: k
!
      if (lsinkpoint) then
        k=1
        do while (k<=npar_loc)
          rp=sqrt((fp(k,ixp)-xsinkpoint)**2+(fp(k,iyp)-ysinkpoint)**2+(fp(k,izp)-zsinkpoint)**2)
          if (rp<rsinkpoint) then
            call remove_particle(fp,ipar,k,dfp,ineargrid)
          else
            k=k+1
          endif
        enddo
      endif
!
!  Remove particles if they are within a solid geometry
!
      if (lsolid_cells) then
        k=1
        do while (k<=npar_loc)
          if (in_solid_cell(fp(k,ixp:izp),fp(k,iap))) then
            call remove_particle(fp,ipar,k,dfp,ineargrid)
          else
            k=k+1
          endif
        enddo
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine remove_particles_sink_simple
!***********************************************************************
    subroutine create_particles_sink_simple(f,fp,dfp,ineargrid)
!
!  Subroutine for creating new sink particles or sink points.
!
!  Just a dummy routine for now.
!
!  25-sep-08/anders: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine create_particles_sink_simple
!***********************************************************************
    subroutine get_frictiontime(f,fp,ineargrid,k,tausp1_par,iblock,uup,nochange_opt)
!
!  Calculate the friction time.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      real, dimension (3), intent(in) :: uup
      integer :: k
      real :: tausp1_par
      integer :: iblock
      logical, optional :: nochange_opt
!
      real :: tmp, epsp, OO
      integer :: ix0, iy0, iz0, jspec
      logical :: nochange=.false.
!
      if (present(nochange_opt)) then
        if (nochange_opt) then
          nochange=.true.
        else
          nochange=.false.
        endif
      else
        nochange=.false.
      endif
!
      ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
!
!  Epstein drag law.
!
      if (ldraglaw_epstein) then
        if (iap/=0) then
          if (fp(k,iap)/=0.0) tausp1_par = tausp01 / fp(k,iap)
        else
          if (npar_species>1) then
            jspec=npar_species*(ipar(k)-1)/npar+1
            tmp=tausp1_species(jspec)
          else
            tmp=tausp1
          endif
!
!  Scale friction time with local density. inearblock(k)=iblock
!
          if (ldraglaw_variable_density) then
            if (ldensity_nolog) then
              tausp1_par=tmp*fb(ix0,iy0,iz0,irho,inearblock(k))
            else
              tausp1_par=tmp*exp(fb(ix0,iy0,iz0,ilnrho,inearblock(k)))
            endif
!
!  Discriminate between constant tau and special case for
!  1/tau=omega when omega is not constant (as, for instance,
!  global Keplerian disks, for which omega=rad**(-3/2)
!
          elseif (ldraglaw_variable) then
            if (lcartesian_coords) then
              OO=(fp(k,ixp)**2 + fp(k,iyp)**2)**(-0.75)
            elseif (lcylindrical_coords) then
              OO=fp(k,ixp)**(-1.5)
            elseif (lspherical_coords) then
              OO=(fp(k,ixp)*sin(fp(k,iyp)))**(-1.5)
            else
              call fatal_error("get_frictiontime","no valid coord system")
              OO=0.
            endif
            tausp1_par=tmp*OO
          else
!
!  Constant friction time
!
            tausp1_par=tmp
          endif
        endif
      else if (ldraglaw_eps_stk_transonic) then
!
         call calc_draglaw_parameters(fp,k,uup,ix0,iy0,iz0,iblock,tausp1_par,lstokes=.true.)
!
      endif
!
!  Change friction time artificially.
!
      if (.not. nochange) then
!
!  Increase friction time linearly with dust density where the dust-to-gas
!  ratio is higher than a chosen value. Supposed to mimick decreased cooling
!  when the gas follows the dust.
!
        if (epsp_friction_increase/=0.0) then
          epsp = fb(ix0,iy0,iz0,irhop,iblock) / get_gas_density_strat(fb,ix0,iy0,iz0,iblock)
          if (epsp>epsp_friction_increase) tausp1_par=tausp1_par/(epsp/epsp_friction_increase)
        endif
!
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine get_frictiontime
!***********************************************************************
    subroutine calc_draglaw_parameters(fp,k,uup,ix0,iy0,iz0,ib,tausp1_par,lstokes)
!
      use EquationOfState, only: rho0,cs0,cs20
!
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension(3) :: uup, duu
      type (pencil_case) :: p
      real :: tausp1_par,tmp,tmp1
      integer :: k, ix0, iy0, iz0, jspec, ib
      real :: kd,fd,mach,mach2,fac,OO,rho,cs2
      real :: knudsen,reynolds,lambda
      real :: inv_particle_radius,kn_crit
      logical, optional :: lstokes
      logical, save :: lfirstcall
!
      intent(in) :: k,uup,ix0,iy0,iz0,ib
      intent(out) :: tausp1_par
!
!  Epstein drag away from the limit of subsonic particle motion. The drag
!  force is given by (Schaaf 1963)
!
!       Feps=-pi*a**2 * rhog * |Delta(u)| * Delta(u) &                    (1)
!        *[(1+1/m**2+1/(4*m**4))*erf(m) + (1/m+1/(2*m**3)*exp(-m**2)/sqrt(pi)]
!
!  where Delta(u) is the relative dust-to-gas velocity (vector)
!  and m=|Delta(u)|/cs (scalar) is the relative mach number of the flow
!
!  As erf is too cumbersome a function to implement numerically, an interpolation
!  between the limits of
!
!     subsonic:    Feps=-sqrt(128*pi)/3*a**2*rhog*cs*Delta(u)             (2)
!     supersonic:  Feps=-pi*a**2*rhog*|Delta(u)|*Delta(u)                 (3)
!
!  is used, leading to an expression that can be used for arbitrary velocities
!  as derived by Kwok (1975).
!
!     transonic:  Feps=-sqrt(128*pi)/3*a**2*rhog*cs*fd*Delta(u)           (4)
!
!  where fd=sqrt(1 + 9*pi/128*m**2)                                       (5)
!
!  The force Feps is divided by the mass of the particle mp=4/3*pi*a**3*rhopmat
!  to yield the acceleration feps=Feps/mp
!
!         feps = -sqrt(8/pi)*rhog*cs*fd*Delta(u)/[a*rhopmat]              (6)
!
!  Epstein drag ceases to work when the particle diameter becomes comparable
!  to the mean free path (lambda) of the gas molecules. In this case, the force
!  is given by Stokes friction in the viscous case (low dust Reynolds number)
!
!      Fsto=-6*pi*a*mu_kin*Delta(u)                                       (7)
!
!  where mu_kin is the kinematic viscosity of the gas
!
!      mu_kin=1/3*rhog*vth*lambda                                         (8)
!
!  and vth=sqrt(8/pi)*cs is the mean thermal velocity of the gas. For high dust
!  Reynolds numbers the viscosity if uninmportant and the drag force of the tur-
!  bulent flow past the particle is given by Newtonian friction
!
!     Fnew=-1.3*pi*a**2*rhog*|Delta(u)|*Delta(u)
!
!  The two cases are once again connected by an interpolating factor
!
!     F'sto=-6*pi*a*kd*mu_kin*Delta(u)
!
!  where kd is a factor that contains the Reynolds number of the flow over the
!  particle (defined in the code, some lines below).
!
!  The following interpolation then works for flows of arbitrary Knudsen, Mach and Reynolds
!  numbers
!
!     Fdrag = [Kn'/(Kn'+1)]**2 * Feps +  [1/(Kn'+1)]**2 * F'sto
!
!  Where Kn'=3*Kn is the critical Knudsen number where the viscous (Stokes) drag and the subsonic
!  Epstein drag are equal.
!
!  (The discussion above was taken from Paardekooper 2006, Woite & Helling 2003 and Kwok 1975)
!
!  In the 2D case, the density rhog is to be replaced by
!
!     rhog=Sigmag/[sqrt(2*pi)H]
!         =Sigmag*Omega/[sqrt(2*pi)*cs]
!
!  which removes the dependence of (6) on cs. We are left with
!
!         feps = -2/pi*sigmag*Omega*fd*Delta(u)/[a*rhopmat]
!
!  the constant terms are tausp1. The same follows for Stokes drag
!
!  Friction time for different species
!
      if (npar_species==1) then
        tmp=tausp
        tmp1=tausp1
      else
        jspec=npar_species*(ipar(k)-1)/npar+1
        tmp=tausp_species(jspec)
        tmp1=tausp1_species(jspec)
      endif
!
!  Relative velocity
!
      duu=fp(k,ivpx:ivpz)-uup
!
!  Sound speed, Mach number and the correction fd to flows of arbitrary mach number
!
      if (llocal_iso) then
        cs2 = cs20/fp(k,ixp)**cs2_powerlaw
      else
        call not_implemented("calc_draglaw_parameters","for other than local isothermal")
      endif
!
      mach=sqrt((duu(1)**2+duu(2)**2+duu(3)**2)/cs2)
      fd=sqrt(1+(9.0*pi/128)*mach**2)
!
!  Angular frequency
!
      if (nzgrid==1) then
!  then omega is needed
        if (ldraglaw_variable) then
          !these omegas assume GM=1
          if (lcartesian_coords) then
            OO=(fp(k,ixp)**2 + fp(k,iyp)**2)**(-0.75)
          elseif (lcylindrical_coords) then
            OO=fp(k,ixp)**(-1.5)
          elseif (lspherical_coords) then
            OO=(fp(k,ixp)*sin(fp(k,iyp)))**(-1.5)
          else
            call fatal_error("calc_draglaw_parameters","no valid coord system")
            OO=0.
          endif
        else
          OO=nu_epicycle
        endif
      endif
!
!  Possibility to include the transition from Epstein to Stokes drag
!
      if (present(lstokes)) then
!
        if (lfirstcall) print*, 'calc_draglaw_parameters: Epstein-Stokes transonic drag law'
!
!  For Stokes drag, the mean free path is needed
!
!   lambda = 1/rhog*(mu/sigma_coll)_H2
!
!  were mu_H2 is the mean molecular weight of the hydrogen molecule (3.9e-24 g),
!  and sigma_coll its cross section (2e-15 cm^2).
!  Assume that (mu/sigma_coll) is the input parameter mean_free_path_gas
!
        if (mean_free_path_gas==0) then
          call fatal_error("calc_draglaw_parameters","You want to use Stokes drag"// &
              "but forgot to set 'mean_free_path_gas' in the *.in files")
        endif
!
        if (ldensity_nolog) then
          rho=fb(ix0,iy0,iz0,irho,ib)
        else
          rho=exp(fb(ix0,iy0,iz0,irho,ib))
        endif
!
        if (nzgrid==1) then
          !the sqrt(2pi) factor is inside the mean_free_path_gas constant
          lambda=mean_free_path_gas * sqrt(cs2)*rho0/(rho*OO*cs0)
        else
          lambda=mean_free_path_gas * rho0/rho
        endif
!
!  The Knudsen number is the ratio of the mean free path to the particle
!  radius, 2s. To keep consistency with the formulation evolving for radius,
!  tausp1 is C/(s*rhopmat) where C is 2/pi for 2d runs and sqrt(8/pi) for 3D
!  runs (because of the sqrt(2*pi) factor coming from the substitution
!  Sigma=rho/(sqrt(2*pi)*H). 's' is the particle radius
!
        if (luse_tau_ap) then
          ! use tausp as the radius (in meter) to make life easier
          inv_particle_radius=tmp1
        else
          if (nzgrid==1) then
            inv_particle_radius=0.5*pi*tmp1     !rhopmat=1, particle_radius in meters
          else
            inv_particle_radius=sqrt(pi/8)*tmp1 !rhopmat=1, particle_radius in meters
          endif
        endif
!
        knudsen=0.5*lambda*inv_particle_radius
!
!  The Stokes drag depends non-linearly on
!
!    Re = 2*s*rho_g*|delta(v)|/mu_kin
!
        reynolds=3*sqrt(pi/8)*mach/knudsen
!
!  the Reynolds number of the flow over the particle. It can parameterized by
!
        if (reynolds<=500) then
          kd=1.0+0.15*reynolds**0.687
        elseif ((reynolds>500).and.(reynolds<=1500)) then
          kd=3.96e-6*reynolds**2.4
        elseif (reynolds>1500) then
          kd=0.11*reynolds
        else
          call fatal_error_local("calc_draglaw_parameters", "'reynolds' seems to be NaN")
          kd=0.
        endif
!
!  And we finally have the Stokes correction to intermediate Knudsen numbers
!  kn_crit is the critical knudsen number where viscous (low reynolds)
!  Stokes and subsonic Epstein friction are equal (Woitke & Helling, 2003)
!
        kn_crit=3*knudsen
        fac=kn_crit/(kn_crit+1)**2 * (kn_crit*fd + kd)
!
      else
!
!  Only use Epstein drag
!
        if (lfirstcall) print*,'get_frictiontime: Epstein transonic drag law'
!
        fac=fd
!
      endif
!
! Calculate tausp1_par for 2d and 3d cases with and without particle_radius
! as a dynamical variable
!
!  Normalize to make tausp1 not dependent on cs0 or rho0
!  Bad because it comes at the expense of evil divisions
!
      if (nzgrid==1) then
        if (luse_tau_ap) then
          tausp1_par=tmp1*2*pi_1*OO*rho*fac/(rho0*rhopmat)
        else
          tausp1_par=tmp1*OO*rho*fac/ rho0
        endif
      else
        if (luse_tau_ap) then
          tausp1_par=tmp1*sqrt(8*pi_1*cs2)*rho*fac/(rho0*cs0)
        else
          tausp1_par=tmp1*sqrt(cs2)*rho*fac/(rho0*cs0)
        endif
      endif
!
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine calc_draglaw_parameters
!***********************************************************************
    subroutine read_particles_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_init_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_init_pars
!***********************************************************************
    subroutine write_particles_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_init_pars)
!
    endsubroutine write_particles_init_pars
!***********************************************************************
    subroutine read_particles_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_run_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_run_pars
!***********************************************************************
    subroutine write_particles_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_run_pars)
!
    endsubroutine write_particles_run_pars
!***********************************************************************
    subroutine powersnap_particles(f)
!
!  Calculate power spectra of dust particle variables.
!
!  01-jan-06/anders: coded
!
      use Power_spectrum, only: power_1d
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (lpar_spec) call power_1d(f,'p',0,irhop)
!
    endsubroutine powersnap_particles
!***********************************************************************
    subroutine rprint_particles(lreset,lwrite)
!
!  Read and register print parameters relevant for particles.
!
!  29-dec-04/anders: coded
!
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname,inamez,inamey,inamex,inamexy,inamexz,inamer,inamerz
      logical :: lwr
!
!  Write information to index.pro.
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset.
!
      if (lreset) then
        idiag_xpm=0; idiag_ypm=0; idiag_zpm=0
        idiag_xp2m=0; idiag_yp2m=0; idiag_zp2m=0; idiag_rpm=0; idiag_rp2m=0
        idiag_vpxm=0; idiag_vpym=0; idiag_vpzm=0
        idiag_vpx2m=0; idiag_vpy2m=0; idiag_vpz2m=0; idiag_ekinp=0
        idiag_vpxmax=0; idiag_vpymax=0; idiag_vpzmax=0; idiag_vpmax=0
        idiag_lpxm=0; idiag_lpym=0; idiag_lpzm=0
        idiag_lpx2m=0; idiag_lpy2m=0; idiag_lpz2m=0
        idiag_npm=0; idiag_np2m=0; idiag_npmax=0; idiag_npmin=0
        idiag_dtdragp=0; idiag_dedragp=0
        idiag_rhopm=0; idiag_rhoprms=0; idiag_rhop2m=0; idiag_rhopmax=0
        idiag_rhopmin=0; idiag_decollp=0; idiag_rhopmphi=0
        idiag_epspmin=0; idiag_epspmax=0
        idiag_nparmin=0; idiag_nparmax=0; idiag_nmigmax=0; idiag_mpt=0
        idiag_npmx=0; idiag_npmy=0; idiag_npmz=0; idiag_epotpm=0
        idiag_rhopmx=0; idiag_rhopmy=0; idiag_rhopmz=0
        idiag_epspmx=0; idiag_epspmy=0; idiag_epspmz=0
        idiag_rpvpxmx = 0; idiag_rpvpymx = 0; idiag_rpvpzmx = 0
        idiag_rpvpx2mx = 0; idiag_rpvpy2mx = 0; idiag_rpvpz2mx = 0
        idiag_rhopmxy=0; idiag_rhopmxz=0; idiag_rhopmr=0; idiag_sigmap=0
        idiag_dvpx2m=0; idiag_dvpy2m=0; idiag_dvpz2m=0
        idiag_dvpmax=0; idiag_dvpm=0; idiag_nparbmax=0
        idiag_nblockmin=0; idiag_nblockmax=0
        idiag_eccpxm=0; idiag_eccpym=0; idiag_eccpzm=0
        idiag_eccpx2m=0; idiag_eccpy2m=0; idiag_eccpz2m=0
        idiag_npargone=0; idiag_deshearbcsm=0
      endif
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot .and. ip<14) print*,'rprint_particles: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'nparmin',idiag_nparmin)
        call parse_name(iname,cname(iname),cform(iname),'nparmax',idiag_nparmax)
        call parse_name(iname,cname(iname),cform(iname),'xpm',idiag_xpm)
        call parse_name(iname,cname(iname),cform(iname),'ypm',idiag_ypm)
        call parse_name(iname,cname(iname),cform(iname),'zpm',idiag_zpm)
        call parse_name(iname,cname(iname),cform(iname),'xp2m',idiag_xp2m)
        call parse_name(iname,cname(iname),cform(iname),'yp2m',idiag_yp2m)
        call parse_name(iname,cname(iname),cform(iname),'zp2m',idiag_zp2m)
        call parse_name(iname,cname(iname),cform(iname),'rpm',idiag_rpm)
        call parse_name(iname,cname(iname),cform(iname),'rp2m',idiag_rp2m)
        call parse_name(iname,cname(iname),cform(iname),'vpxm',idiag_vpxm)
        call parse_name(iname,cname(iname),cform(iname),'vpym',idiag_vpym)
        call parse_name(iname,cname(iname),cform(iname),'vpzm',idiag_vpzm)
        call parse_name(iname,cname(iname),cform(iname),'vpx2m',idiag_vpx2m)
        call parse_name(iname,cname(iname),cform(iname),'vpy2m',idiag_vpy2m)
        call parse_name(iname,cname(iname),cform(iname),'vpz2m',idiag_vpz2m)
        call parse_name(iname,cname(iname),cform(iname),'ekinp',idiag_ekinp)
        call parse_name(iname,cname(iname),cform(iname),'vpxmax',idiag_vpxmax)
        call parse_name(iname,cname(iname),cform(iname),'vpymax',idiag_vpymax)
        call parse_name(iname,cname(iname),cform(iname),'vpzmax',idiag_vpzmax)
        call parse_name(iname,cname(iname),cform(iname),'vpmax',idiag_vpmax)
        call parse_name(iname,cname(iname),cform(iname),'lpxm',idiag_lpxm)
        call parse_name(iname,cname(iname),cform(iname),'lpym',idiag_lpym)
        call parse_name(iname,cname(iname),cform(iname),'lpzm',idiag_lpzm)
        call parse_name(iname,cname(iname),cform(iname),'lpx2m',idiag_lpx2m)
        call parse_name(iname,cname(iname),cform(iname),'lpy2m',idiag_lpy2m)
        call parse_name(iname,cname(iname),cform(iname),'lpz2m',idiag_lpz2m)
        call parse_name(iname,cname(iname),cform(iname),'eccpxm',idiag_eccpxm)
        call parse_name(iname,cname(iname),cform(iname),'eccpym',idiag_eccpym)
        call parse_name(iname,cname(iname),cform(iname),'eccpzm',idiag_eccpzm)
        call parse_name(iname,cname(iname),cform(iname),'eccpx2m',idiag_eccpx2m)
        call parse_name(iname,cname(iname),cform(iname),'eccpy2m',idiag_eccpy2m)
        call parse_name(iname,cname(iname),cform(iname),'eccpz2m',idiag_eccpz2m)
        call parse_name(iname,cname(iname),cform(iname),'dtdragp',idiag_dtdragp)
        call parse_name(iname,cname(iname),cform(iname),'npm',idiag_npm)
        call parse_name(iname,cname(iname),cform(iname),'np2m',idiag_np2m)
        call parse_name(iname,cname(iname),cform(iname),'npmax',idiag_npmax)
        call parse_name(iname,cname(iname),cform(iname),'npmin',idiag_npmin)
        call parse_name(iname,cname(iname),cform(iname),'rhopm',idiag_rhopm)
        call parse_name(iname,cname(iname),cform(iname),'rhoprms',idiag_rhoprms)
        call parse_name(iname,cname(iname),cform(iname),'rhop2m',idiag_rhop2m)
        call parse_name(iname,cname(iname),cform(iname),'rhopmin',idiag_rhopmin)
        call parse_name(iname,cname(iname),cform(iname),'rhopmax',idiag_rhopmax)
        call parse_name(iname,cname(iname),cform(iname),'epspmin',idiag_epspmin)
        call parse_name(iname,cname(iname),cform(iname),'epspmax',idiag_epspmax)
        call parse_name(iname,cname(iname),cform(iname),'rhopmphi',idiag_rhopmphi)
        call parse_name(iname,cname(iname),cform(iname),'nmigmax',idiag_nmigmax)
        call parse_name(iname,cname(iname),cform(iname),'mpt',idiag_mpt)
        call parse_name(iname,cname(iname),cform(iname),'dvpx2m',idiag_dvpx2m)
        call parse_name(iname,cname(iname),cform(iname),'dvpy2m',idiag_dvpy2m)
        call parse_name(iname,cname(iname),cform(iname),'dvpz2m',idiag_dvpz2m)
        call parse_name(iname,cname(iname),cform(iname),'dvpm',idiag_dvpm)
        call parse_name(iname,cname(iname),cform(iname),'dvpmax',idiag_dvpmax)
        call parse_name(iname,cname(iname),cform(iname),'dedragp',idiag_dedragp)
        call parse_name(iname,cname(iname),cform(iname),'decollp',idiag_decollp)
        call parse_name(iname,cname(iname),cform(iname),'epotpm',idiag_epotpm)
        call parse_name(iname,cname(iname),cform(iname),'npargone',idiag_npargone)
        call parse_name(iname,cname(iname),cform(iname),'nparbmax',idiag_nparbmax)
        call parse_name(iname,cname(iname),cform(iname),'nblockmin',idiag_nblockmin)
        call parse_name(iname,cname(iname),cform(iname),'nblockmax',idiag_nblockmax)
        call parse_name(iname,cname(iname),cform(iname),'deshearbcsm',idiag_deshearbcsm)
      enddo
!
!  Check for those quantities for which we want x-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'npmx',idiag_npmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'rhopmx',idiag_rhopmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'epspmx',idiag_epspmx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), 'rpvpxmx', idiag_rpvpxmx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), 'rpvpymx', idiag_rpvpymx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), 'rpvpzmx', idiag_rpvpzmx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), 'rpvpx2mx', idiag_rpvpx2mx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), 'rpvpy2mx', idiag_rpvpy2mx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), 'rpvpz2mx', idiag_rpvpz2mx)
      enddo
!
!  Check for those quantities for which we want y-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'npmy',idiag_npmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'rhopmy',idiag_npmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'epspmy',idiag_epspmy)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'npmz',idiag_npmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rhopmz',idiag_rhopmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'epspmz',idiag_epspmz)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'rhopmxy',idiag_rhopmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'sigmap',idiag_sigmap)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamexz=1,nnamexz
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'rhopmxz',idiag_rhopmxz)
      enddo
!
!  Check for those quantities for which we want phiz-averages.
!
      do inamer=1,nnamer
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'rhopmr',idiag_rhopmr)
      enddo
!
!  Check for those quantities for which we want phi-averages.
!
      do inamerz=1,nnamerz
        call parse_name(inamerz,cnamerz(inamerz),cformrz(inamerz),'rhopmphi',idiag_rhopmphi)
      enddo
!
    endsubroutine rprint_particles
!***********************************************************************
    real function get_gas_density_strat(fb, ix, iy, iz, ib) result(rho)
!
!  Reads the gas density at one location in a block.
!
!  04-mar-13/ccyang: coded.
!
      use EquationOfState, only: get_stratz
!
      real, dimension(mxb,myb,mzb,mfarray,0:nblockmax-1), intent(in) :: fb
      integer, intent(in) :: ix, iy, iz, ib
!
      real, dimension(1) :: rho0
!
      stratz: if (lstratz) then
        call get_stratz((/zb(iz,ib)/), rho0z=rho0)
        rho = rho0(1) * (1.0 + fb(ix, iy, iz, irho, ib))
      elseif (ldensity_nolog) then stratz
        rho = fb(ix, iy, iz, irho, ib)
      else stratz
        rho = exp(fb(ix, iy, iz, ilnrho, ib))
      endif stratz
!
    endfunction get_gas_density_strat
!***********************************************************************
    subroutine periodic_boundcond_on_aux(f)
!
! dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine periodic_boundcond_on_aux
!***********************************************************************
endmodule Particles
