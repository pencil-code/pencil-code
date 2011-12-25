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
!
! PENCILS PROVIDED np; rhop; epsp; rhop_swarm
!
!***************************************************************
module Particles
!
  use Cdata
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
  use Particles_sub
  use Particles_radius
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'particles.h'
!
  real, dimension (npar_species) :: tausp_species=0.0, tausp1_species=0.0
  real :: xp0=0.0, yp0=0.0, zp0=0.0, vpx0=0.0, vpy0=0.0, vpz0=0.0
  real :: Lx0=0.0, Ly0=0.0, Lz0=0.0
  real :: delta_vp0=1.0, tausp=0.0, tausp1=0.0, eps_dtog=0.01
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
  real :: pdlaw=0.0
  real :: xsinkpoint=0.0, ysinkpoint=0.0, zsinkpoint=0.0, rsinkpoint=0.0
  logical :: ldragforce_dust_par=.false., ldragforce_gas_par=.false.
  logical :: lpar_spec=.false., learly_particle_map=.true.
  logical :: ldragforce_equi_global_eps=.false.
  logical :: ldraglaw_epstein=.true.
  logical :: ldt_grav_par=.false., ldt_adv_par=.true.
  logical :: lsinkpoint=.false., lglobalrandom=.false.
  logical :: lcoriolis_force_par=.true., lcentrifugal_force_par=.false.
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
      lparticlemesh_cic, lparticlemesh_tsc, linterpolate_spline, &
      tstart_dragforce_par, tstart_grav_par, tausp_species, &
      learly_particle_map, epsp_friction_increase, lmigration_real_check, &
      ldraglaw_epstein, lcheck_exact_frontier, pdlaw, ldt_grav_par, &
      lsinkpoint,xsinkpoint, ysinkpoint, zsinkpoint, rsinkpoint, &
      lcoriolis_force_par, lcentrifugal_force_par, ldt_adv_par, Lx0, Ly0, &
      Lz0, lglobalrandom, linsert_particles_continuously, &
      lrandom_particle_pencils, lnocalc_np, lnocalc_rhop, it1_loadbalance, &
      np_const, rhop_const, lrandom_particle_blocks, lreblock_particles_run, &
      lbrick_partition
!
  namelist /particles_run_pars/ &
      bcpx, bcpy, bcpz, tausp, dsnap_par_minor, beta_dPdr_dust, &
      ldragforce_gas_par, ldragforce_dust_par, mp_swarm, np_swarm, &
      eps_dtog, cdtp, cdtpgrav, lpar_spec, linterp_reality_check, &
      nu_epicycle, gravx_profile, gravz_profile, gravr_profile, gravx, gravz, &
      gravr, gravsmooth, kx_gg, kz_gg, lmigration_redo, tstart_dragforce_par, &
      tstart_grav_par, lparticlemesh_cic, lparticlemesh_tsc, &
      epsp_friction_increase, learly_particle_map, lmigration_real_check, &
      ldraglaw_epstein, lcheck_exact_frontier, ldt_grav_par, lsinkpoint, &
      xsinkpoint, ysinkpoint, zsinkpoint, rsinkpoint, lcoriolis_force_par, &
      lcentrifugal_force_par, ldt_adv_par, linsert_particles_continuously, &
      lrandom_particle_pencils, lnocalc_np, lnocalc_rhop, it1_loadbalance, &
      np_const, rhop_const, lrandom_particle_blocks, lreblock_particles_run, &
      lbrick_partition
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
  integer :: idiag_rhoptilm=0, idiag_dtdragp=0
  integer :: idiag_nparmin=0, idiag_nparmax=0, idiag_npargone=0
  integer :: idiag_rhopm=0, idiag_rhoprms=0, idiag_rhop2m=0, idiag_rhopmax=0
  integer :: idiag_rhopmin=0, idiag_decollp=0, idiag_rhopmphi=0
  integer :: idiag_epspmin=0, idiag_epspmax=0
  integer :: idiag_npmx=0, idiag_npmy=0, idiag_npmz=0
  integer :: idiag_rhopmx=0, idiag_rhopmy=0, idiag_rhopmz=0
  integer :: idiag_epspmx=0, idiag_epspmy=0, idiag_epspmz=0
  integer :: idiag_mpt=0, idiag_dedragp=0, idiag_rhopmxy=0, idiag_rhopmr=0
  integer :: idiag_dvpx2m=0, idiag_dvpy2m=0, idiag_dvpz2m=0
  integer :: idiag_dvpm=0, idiag_dvpmax=0, idiag_epotpm=0
  integer :: idiag_rhopmxz=0, idiag_nparbmax=0
  integer :: idiag_eccpxm=0, idiag_eccpym=0, idiag_eccpzm=0
  integer :: idiag_eccpx2m=0, idiag_eccpy2m=0, idiag_eccpz2m=0
!
  contains
!***********************************************************************
    subroutine register_particles()
!
!  Set up indices for access to the fp and dfp arrays
!
!  29-dec-04/anders: coded
!
      use FArrayManager
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Indices for particle position.
!
      ixp=npvar+1
      iyp=npvar+2
      izp=npvar+3
!
!  Indices for particle velocity.
!
      ivpx=npvar+4
      ivpy=npvar+5
      ivpz=npvar+6
!
!  Increase npvar accordingly.
!
      npvar=npvar+6
!
!  Set indices for auxiliary variables
!
      call farray_register_auxiliary('np',inp)
      call farray_register_auxiliary('rhop',irhop)
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_particles','npvar > mpvar')
      endif
!
    endsubroutine register_particles
!***********************************************************************
    subroutine initialize_particles(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  29-dec-04/anders: coded
!
      use EquationOfState, only: cs0,rho0
      use FArrayManager
      use SharedVariables, only: put_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      real :: rhom
      integer :: ierr, jspec
!
!  This module is incompatible with normal domain decomposition.
!
      if (.not.lparticles_blocks) then
        if (lroot) then
          print*, 'initialize_particles: must use PARTICLES =                   PARTICLES_DUST'
          print*, '                      for normal domain decomposition'
        endif
        call fatal_error('initialize_particles','')
      endif
!
!  The inverse stopping time is needed for drag force and collisional cooling.
!
      if (tausp/=0.0) tausp1=1/tausp
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
          print*, 'initialize_particles: '// &
              'Number of particle species = ', npar_species
          print*, 'initialize_particles: tausp_species = ', tausp_species
        endif
!
!  Must have set tausp_species for drag force.
!
        if (ldragforce_dust_par .or. ldragforce_gas_par) then
          if (any(tausp_species==0.0)) then
            if (lroot) print*, &
                'initialize_particles: drag force must have tausp_species/=0 !'
                call fatal_error('initialize_particles','')
          endif
!
!  Inverse friction time is needed for drag force.
!
          do jspec=1,npar_species
            if (tausp_species(jspec)/=0.0) &
                tausp1_species(jspec)=1/tausp_species(jspec)
          enddo
        endif
      else
        tausp_species(1)=tausp
        if (tausp_species(1)/=0.0) &
            tausp1_species(1)=1/tausp_species(1)
      endif
!
!  Share friction time (but only if Epstein drag regime!).
!
      if (ldraglaw_epstein) then
        call put_shared_variable( 'tausp_species', tausp_species)
        call put_shared_variable('tausp1_species',tausp1_species)
      endif
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
        if ((lgravz .and. lgravz_gas) .or. gravz_profile=='linear' ) then
          rhom=sqrt(2*pi)*1.0*1.0/Lz  ! rhom = Sigma/Lz, Sigma=sqrt(2*pi)*H*rho1
        else
          rhom=1.0
        endif
        if (rhop_swarm==0.0) &
             rhop_swarm = eps_dtog*rhom/(real(npar)/(nxgrid*nygrid*nzgrid))
        if (mp_swarm==0.0) &
             mp_swarm   = eps_dtog*rhom*box_volume/(real(npar))
        if (lroot) print*, 'initialize_particles: '// &
            'dust-to-gas ratio eps_dtog=', eps_dtog
      endif
!
      if (lroot) then
        print*, 'initialize_particles: '// &
            'mass per constituent particle mpmat=', mpmat
        print*, 'initialize_particles: '// &
            'mass per superparticle mp_swarm =', mp_swarm
        print*, 'initialize_particles: '// &
            'number density per superparticle np_swarm=', np_swarm
        print*, 'initialize_particles: '// &
            'mass density per superparticle rhop_swarm=', rhop_swarm
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
!  Share Keplerian gravity.
!
      call put_shared_variable('gravr',gravr,ierr)
      if (ierr/=0) call fatal_error('initialize_particles', &
          'there was a problem when sharing gravr')
!
!  Gas density is needed for back-reaction friction force.
!
      if (ldragforce_gas_par .and. .not. ldensity) then
        if (lroot) then
          print*, 'initialize_particles: friction force on gas only works '
          print*, '                      together with gas density module!'
        endif
        call fatal_error('initialize_particles','')
      endif
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
!  Drag force on gas right now assumed rhop_swarm is the same for all particles.
!
      if (ldragforce_gas_par.and.(lparticles_radius.or.lparticles_number)) then
        if (lroot) print*, 'initialize_particles: drag force on gas is '// &
            'not yet implemented for variable particle radius or number'
        call fatal_error('initialize_particles','')
      endif
!
!  Write constants to disk.
!
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro',position='append')
          write (1,*) 'np_swarm=', np_swarm
          write (1,*) 'mpmat=', mpmat
          write (1,*) 'mp_swarm=', mp_swarm
        close (1)
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles
!***********************************************************************
    subroutine init_particles(f,fp,ineargrid)
!
!  Initial positions and velocities of dust particles.
!
!  29-dec-04/anders: coded
!
      use EquationOfState, only: gamma, beta_glnrho_global, cs20
      use General, only: random_number_wrapper
      use Mpicomm, only: mpireduce_sum
      use Sub
      use InitialCondition, only: initial_condition_xxp,&
                                  initial_condition_vvp
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: Lxyz_par, xyz0_par, xyz1_par
      real :: vpx_sum, vpy_sum, vpz_sum
      real :: r, p, q, px, py, pz, eps, cs, k2_xxp, rp2
      real :: dim1, npar_loc_x, npar_loc_y, npar_loc_z, dx_par, dy_par, dz_par
      real :: rad,rad_scl,phi,tmp,OO,xx0,yy0,r2
      integer :: l, j, k, ix0, iy0, iz0, ib
      logical :: lequidistant=.false.
!
      intent (out) :: f, fp, ineargrid
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
          if (lroot) &
              print*, 'init_particles: All particles at x,y,z=', xp0, yp0, zp0
          fp(1:npar_loc,ixp)=xp0
          fp(1:npar_loc,iyp)=yp0
          fp(1:npar_loc,izp)=zp0
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
        case ('random-hole')
          if (lroot) print*, 'init_particles: Random particle positions '// &
              'with inner hole'
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
          if (lroot) print*, 'init_particles: Random particle positions '// &
               'within a box'
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
       case ('random-cylindrical','random-cyl')
          if (lroot) print*, 'init_particles: Random particle '//&
               'cylindrical positions with power-law pdlaw=',pdlaw
!
          do k=1,npar_loc
!
! Start the particles obeying a power law pdlaw
!
            tmp=2-pdlaw
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
              call fatal_error('init_particles','random-cylindrical '// &
                  'not implemented for spherical coordinates')
            endif
!
            if (nzgrid/=1) call random_number_wrapper(fp(k,izp))
            if (nzgrid/=1) &
                fp(k,izp)=xyz0_par(3)+fp(k,izp)*Lxyz_par(3)
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
          if (.not. lequidistant) then
            if (lroot) print*, 'init_particles: must place particles equidistantly before shifting!'
            call fatal_error('init_particles','')
          endif
          k2_xxp=kx_xxp**2+ky_xxp**2+kz_xxp**2
          if (k2_xxp==0.0) then
            if (lroot) print*, &
                'init_particles: kx_xxp=ky_xxp=kz_xxp=0.0 is not allowed!'
            call fatal_error('init_particles','')
          endif
          do k=1,npar_loc
            fp(k,ixp) = fp(k,ixp) - kx_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp)+ky_xxp*fp(k,iyp)+kz_xxp*fp(k,izp))
            fp(k,iyp) = fp(k,iyp) - ky_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp)+ky_xxp*fp(k,iyp)+kz_xxp*fp(k,izp))
            fp(k,izp) = fp(k,izp) - kz_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp)+ky_xxp*fp(k,iyp)+kz_xxp*fp(k,izp))
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
          if (lroot) &
              print*, 'init_particles: No such such value for initxxp: ', &
              trim(initxxp(j))
          call fatal_error('init_particles','')
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
          if (lroot) print*, 'init_particles: for random particle velocity '// &
              'initial conditions you should use learly_particle_map=F'
          call fatal_error('init_particles','')
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
          if (lroot) &
              print*, 'init_particles: vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
          if (lcylindrical_coords) then
            fp(1:npar_loc,ivpx)&
                =vpx0*cos(fp(k,iyp))+vpy0*sin(fp(k,iyp))
            fp(1:npar_loc,ivpy)&
                =vpy0*cos(fp(k,iyp))-vpx0*sin(fp(k,iyp))
            fp(1:npar_loc,ivpz)=vpz0
          else
            fp(1:npar_loc,ivpx)=vpx0
            fp(1:npar_loc,ivpy)=vpy0
            fp(1:npar_loc,ivpz)=vpz0
          endif
!
        case ('sinwave-phase')
          if (lroot) print*, 'init_particles: sinwave-phase'
          if (lroot) &
              print*, 'init_particles: vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
          do k=1,npar_loc
            fp(k,ivpx)=fp(k,ivpx)+vpx0*sin(kx_vpx*fp(k,ixp)+ky_vpx*fp(k,iyp)+kz_vpx*fp(k,izp)+phase_vpx)
            fp(k,ivpy)=fp(k,ivpy)+vpy0*sin(kx_vpy*fp(k,ixp)+ky_vpy*fp(k,iyp)+kz_vpy*fp(k,izp)+phase_vpy)
            fp(k,ivpz)=fp(k,ivpz)+vpz0*sin(kx_vpz*fp(k,ixp)+ky_vpz*fp(k,iyp)+kz_vpz*fp(k,izp)+phase_vpz)
          enddo
!
        case ('coswave-phase')
          if (lroot) print*, 'init_particles: coswave-phase'
          if (lroot) &
              print*, 'init_particles: vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
          do k=1,npar_loc
            fp(k,ivpx)=fp(k,ivpx)+vpx0*cos(kx_vpx*fp(k,ixp)+ky_vpx*fp(k,iyp)+kz_vpx*fp(k,izp)+phase_vpx)
            fp(k,ivpy)=fp(k,ivpy)+vpy0*cos(kx_vpy*fp(k,ixp)+ky_vpy*fp(k,iyp)+kz_vpy*fp(k,izp)+phase_vpy)
            fp(k,ivpz)=fp(k,ivpz)+vpz0*cos(kx_vpz*fp(k,ixp)+ky_vpz*fp(k,iyp)+kz_vpz*fp(k,izp)+phase_vpz)
          enddo
!
        case ('random')
          if (lroot) print*, 'init_particles: Random particle velocities; '// &
              'delta_vp0=', delta_vp0
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
          if (lroot) print*, 'init_particles: Random particle x-velocity; '// &
              'delta_vp0=', delta_vp0
          do k=1,npar_loc
            call random_number_wrapper(r)
            fp(k,ivpx) = fp(k,ivpx) + delta_vp0*(2*r-1)
          enddo
!
        case ('random-y')
          if (lroot) print*, 'init_particles: Random particle y-velocity; '// &
              'delta_vp0=', delta_vp0
          do k=1,npar_loc
            call random_number_wrapper(r)
            fp(k,ivpy) = fp(k,ivpy) + delta_vp0*(2*r-1)
          enddo
!
        case ('random-z')
          if (lroot) print*, 'init_particles: Random particle z-velocity; '// &
              'delta_vp0=', delta_vp0
          do k=1,npar_loc
            call random_number_wrapper(r)
            fp(k,ivpz) = fp(k,ivpz) + delta_vp0*(2*r-1)
          enddo
!
        case ('average-to-zero')
          call mpireduce_sum(sum(fp(1:npar_loc,ivpx)),vpx_sum)
          fp(1:npar_loc,ivpx)=fp(1:npar_loc,ivpx)-vpx_sum/npar
          call mpireduce_sum(sum(fp(1:npar_loc,ivpy)),vpy_sum)
          fp(1:npar_loc,ivpy)=fp(1:npar_loc,ivpy)-vpy_sum/npar
          call mpireduce_sum(sum(fp(1:npar_loc,ivpz)),vpz_sum)
          fp(1:npar_loc,ivpz)=fp(1:npar_loc,ivpz)-vpz_sum/npar
!
        case ('jeans-wave-dustpar-x')
        ! assumes rhs_poisson_const=1 !
          do k=1,npar_loc
            fp(k,ivpx) = fp(k,ivpx) - amplxxp* &
                (sqrt(1+4*1.0*1.0*tausp**2)-1)/ &
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
          call fill_blocks_with_bricks(f,fb,mfarray,ilnrho)
!  Calculate average dust-to-gas ratio in box.
          if (ldensity_nolog) then
            eps = sum(f(l1:l2,m1:m2,n1:n2,irhop))/ &
                sum(f(l1:l2,m1:m2,n1:n2,irho))
          else
            eps = sum(f(l1:l2,m1:m2,n1:n2,irhop))/ &
                sum(exp(f(l1:l2,m1:m2,n1:n2,ilnrho)))
          endif
!
          if (lroot) print*, 'init_particles: average dust-to-gas ratio=', eps
!  Set gas velocity field.
          do l=l1,l2; do m=m1,m2; do n=n1,n2
!  Take either global or local dust-to-gas ratio.
            if (.not. ldragforce_equi_global_eps) then
              if (ldensity_nolog) then
                eps=f(l,m,n,irhop)/f(l,m,n,irho)
              else
                eps=f(l,m,n,irhop)/exp(f(l,m,n,ilnrho))
              endif
            endif
!
            f(l,m,n,iux) = f(l,m,n,iux) - &
                beta_glnrho_global(1)*eps*Omega*tausp/ &
                ((1.0+eps)**2+(Omega*tausp)**2)*cs
            f(l,m,n,iuy) = f(l,m,n,iuy) + &
                beta_glnrho_global(1)*(1+eps+(Omega*tausp)**2)/ &
                (2*((1.0+eps)**2+(Omega*tausp)**2))*cs
!
          enddo; enddo; enddo
!  Set particle velocity field.
          do k=1,npar_loc
!  Take either global or local dust-to-gas ratio.
            if (.not. ldragforce_equi_global_eps) then
              ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
              ib=inearblock(k)
              if (ldensity_nolog) then
                eps=fb(ix0,iy0,iz0,irhop,ib)/fb(ix0,iy0,iz0,irho,ib)
              else
                eps=fb(ix0,iy0,iz0,irhop,ib)/exp(fb(ix0,iy0,iz0,ilnrho,ib))
              endif
            endif
!
            fp(k,ivpx) = fp(k,ivpx) + &
                beta_glnrho_global(1)*Omega*tausp/ &
                ((1.0+eps)**2+(Omega*tausp)**2)*cs
            fp(k,ivpy) = fp(k,ivpy) + &
                beta_glnrho_global(1)*(1+eps)/ &
                (2*((1.0+eps)**2+(Omega*tausp)**2))*cs
!
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
            fp(k,ivpx) = fp(k,ivpx) + &
                1/gamma*beta_dPdr_dust/ &
                (Omega*tausp+1/(Omega*tausp))*cs
            fp(k,ivpy) = fp(k,ivpy) - &
                1/gamma*beta_dPdr_dust*Omega*tausp*0.5/ &
                (Omega*tausp+1/(Omega*tausp))*cs
          enddo
!
       case ('Keplerian','keplerian')
!
!  Keplerian velocity based on gravr.
!
          if (lroot) then
            print*, 'init_particles: Keplerian velocity'
            if (lspherical_coords) call fatal_error('init_particles', &
                 'Keplerian particle initial condition: '// &
                 'not implemented for spherical coordinates')
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
          if (lroot) &
              print*, 'init_particles: No such such value for initvvp: ', &
              trim(initvvp(j))
          call fatal_error('','')
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
    subroutine insert_particles(f,fp,ineargrid)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
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
      real, dimension (mpar_loc,mpvar) :: fp
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
      if (idiag_rhopmxy/=0 .or. idiag_rhopmr/=0 .or. idiag_rhopmxz/=0) &
          lpenc_diagnos2d(i_rhop)=.true.
      if (idiag_dedragp/=0 .or. idiag_decollp/=0) then
        lpenc_diagnos(i_TT1)=.true.
        lpenc_diagnos(i_rho1)=.true.
      endif
      if (idiag_epspmx/=0 .or. idiag_epspmy/=0 .or. idiag_epspmz/=0 .or. &
          idiag_epspmin/=0 .or. idiag_epspmax/=0) &
          lpenc_diagnos(i_epsp)=.true.
      if (idiag_rhopmxy/=0 .or. idiag_rhopmxz/=0) lpenc_diagnos2d(i_rhop)=.true.
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
      if (lpencil_in(i_epsp)) then
        lpencil_in(i_rhop)=.true.
        lpencil_in(i_rho1)=.true.
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
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      if (lpencil(i_np)) p%np=f(l1:l2,m,n,inp)
!
      if (lpencil(i_rhop)) then
        if (irhop/=0) then
          p%rhop=f(l1:l2,m,n,irhop)
        else
          call get_rhopswarm(mp_swarm,m,n,p%rhop_swarm)
          p%rhop=p%rhop_swarm*f(l1:l2,m,n,inp)
        endif
      endif
!
      if (lpencil(i_epsp)) p%epsp=p%rhop*p%rho1
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
      real, dimension (mpar_loc,mpvar) :: fp, dfp
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
      endif
!
      if (lheader) print*, 'dxxp_dt: Set rate of change of particle '// &
          'position equal to particle velocity.'
!
!  The rate of change of a particle's position is the particle's velocity.
!
      if (lcartesian_coords) then
!
        if (nxgrid/=1) &
             dfp(1:npar_loc,ixp) = dfp(1:npar_loc,ixp) + fp(1:npar_loc,ivpx)
        if (nygrid/=1) &
             dfp(1:npar_loc,iyp) = dfp(1:npar_loc,iyp) + fp(1:npar_loc,ivpy)
        if (nzgrid/=1) &
             dfp(1:npar_loc,izp) = dfp(1:npar_loc,izp) + fp(1:npar_loc,ivpz)
!
      elseif (lcylindrical_coords) then
!
        if (nxgrid/=1) &
             dfp(1:npar_loc,ixp) = dfp(1:npar_loc,ixp) + fp(1:npar_loc,ivpx)
        if (nygrid/=1) &
             dfp(1:npar_loc,iyp) = dfp(1:npar_loc,iyp) + &
             fp(1:npar_loc,ivpy)/max(fp(1:npar_loc,ixp),tini)
        if (nzgrid/=1) &
             dfp(1:npar_loc,izp) = dfp(1:npar_loc,izp) + fp(1:npar_loc,ivpz)
!
      elseif (lspherical_coords) then
!
        if (nxgrid/=1) &
             dfp(1:npar_loc,ixp) = dfp(1:npar_loc,ixp) + fp(1:npar_loc,ivpx)
        if (nygrid/=1) &
             dfp(1:npar_loc,iyp) = dfp(1:npar_loc,iyp) + &
             fp(1:npar_loc,ivpy)/max(fp(1:npar_loc,ixp),tini)
        if (nzgrid/=1) &
             dfp(1:npar_loc,izp) = dfp(1:npar_loc,izp) + &
             fp(1:npar_loc,ivpz)/(max(fp(1:npar_loc,ixp),tini)*&
                                                    sin(fp(1:npar_loc,iyp)))
!
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
      use EquationOfState, only: cs20, gamma
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension(3) :: ggp
      real :: Omega2, rsph, vsph, OO2
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
            print*,'dvvp_dt: Coriolis force on the particles is '
            print*,'not yet implemented for spherical coordinates.'
            call fatal_error('dvvp_dt','')
          endif
        endif
!
! Centrifugal force
!
        if (lcentrifugal_force_par) then
          if (lheader) print*,'dvvp_dt: Add Centrifugal force; Omega=', Omega
          if (lcartesian_coords) then
!
            dfp(1:npar_loc,ivpx) = dfp(1:npar_loc,ivpx) + &
                Omega**2*fp(1:npar_loc,ixp)
!
            dfp(1:npar_loc,ivpy) = dfp(1:npar_loc,ivpy) + &
                Omega**2*fp(1:npar_loc,iyp)
!
          elseif (lcylindrical_coords) then
            dfp(1:npar_loc,ivpx) = &
                dfp(1:npar_loc,ivpx) + Omega**2*fp(1:npar_loc,ixp)
          else
            print*,'dvvp_dt: Centrifugal force on the particles is '
            print*,'not implemented for spherical coordinates.'
            call fatal_error('dvvp_dt','')
          endif
        endif
!
!  With shear there is an extra term due to the background shear flow.
!
        if (lshear) dfp(1:npar_loc,ivpy) = &
            dfp(1:npar_loc,ivpy) + qshear*Omega*fp(1:npar_loc,ivpx)
      endif
!
!  Add constant background pressure gradient beta=alpha*H0/r0, where alpha
!  comes from a global pressure gradient P = P0*(r/r0)^alpha.
!  (the term must be added to the dust equation of motion when measuring
!  velocities relative to the shear flow modified by the global pressure grad.)
!
      if (beta_dPdr_dust/=0.0 .and. t>=tstart_dragforce_par) then
        dfp(1:npar_loc,ivpx) = &
            dfp(1:npar_loc,ivpx) + 1/gamma*cs20*beta_dPdr_dust_scaled
      endif
!
!  Gravity on the particles.
!
      if (t>=tstart_grav_par) then
!
!  Gravity in the x-direction.
!
        select case (gravx_profile)
!
          case ('')
            if (lheader) print*, 'dvvp_dt: No gravity in x-direction.'
!
          case ('zero')
            if (lheader) print*, 'dvvp_dt: No gravity in x-direction.'
!
          case ('linear')
            if (lheader) print*, 'dvvp_dt: Linear gravity field in x-direction.'
            dfp(1:npar_loc,ivpx)=dfp(1:npar_loc,ivpx) - &
                nu_epicycle2*fp(1:npar_loc,ixp)
!
          case ('plain')
            if (lheader) print*, 'dvvp_dt: Plain gravity field in x-direction.'
            dfp(1:npar_loc,ivpx)=dfp(1:npar_loc,ivpx) - gravx
!
          case ('sinusoidal')
            if (lheader) &
                print*, 'dvvp_dt: Sinusoidal gravity field in x-direction.'
            dfp(1:npar_loc,ivpx)=dfp(1:npar_loc,ivpx) - &
                gravx*sin(kx_gg*fp(1:npar_loc,ixp))
!
          case default
            call fatal_error('dvvp_dt','chosen gravx_profile is not valid!')
!
        endselect
!
!  Gravity in the z-direction.
!
        select case (gravz_profile)
!
          case ('')
            if (lheader) print*, 'dvvp_dt: No gravity in z-direction.'
!
          case ('zero')
            if (lheader) print*, 'dvvp_dt: No gravity in z-direction.'
!
          case ('linear')
            if (lheader) print*, 'dvvp_dt: Linear gravity field in z-direction.'
            dfp(1:npar_loc,ivpz)=dfp(1:npar_loc,ivpz) - &
                nu_epicycle2*fp(1:npar_loc,izp)
!
          case ('sinusoidal')
            if (lheader) &
                print*, 'dvvp_dt: Sinusoidal gravity field in z-direction.'
            dfp(1:npar_loc,ivpz)=dfp(1:npar_loc,ivpz) - &
                gravz*sin(kz_gg*fp(1:npar_loc,izp))
!
          case default
            call fatal_error('dvvp_dt','chosen gravz_profile is not valid!')
!
        endselect
!
!  Radial gravity.
!
        select case (gravr_profile)
!
        case ('')
          if (lheader) print*, 'dvvp_dt: No radial gravity'
!
        case ('zero')
          if (lheader) print*, 'dvvp_dt: No radial gravity'
!
        case ('newtonian-central','newtonian')
          if (lparticles_nbody) &
              call fatal_error('dvvp_dt','You are using massive particles. '//&
              'The N-body code should take care of the stellar-like '// &
              'gravity on the dust. Switch off the '// &
              'gravr_profile=''newtonian'' on particles_init')
          if (lheader) &
               print*, 'dvvp_dt: Newtonian gravity from a fixed central object'
          do k=1,npar_loc
            if (lcartesian_coords) then
              rsph=sqrt(fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2+gravsmooth2)
              OO2=rsph**(-3)*gravr
              ggp(1) = -fp(k,ixp)*OO2
              ggp(2) = -fp(k,iyp)*OO2
              ggp(3) = -fp(k,izp)*OO2
              dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + ggp
            elseif (lcylindrical_coords) then
              rsph=sqrt(fp(k,ixp)**2+fp(k,izp)**2+gravsmooth2)
              OO2=rsph**(-3)*gravr
              ggp(1) = -fp(k,ixp)*OO2
              ggp(2) = 0.0
              ggp(3) = -fp(k,izp)*OO2
              dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + ggp
            elseif (lspherical_coords) then
              rsph=sqrt(fp(k,ixp)**2+gravsmooth2)
              OO2=rsph**(-3)*gravr
              ggp(1) = -fp(k,ixp)*OO2
              ggp(2) = 0.0; ggp(3) = 0.0
              dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + ggp
            endif
!  Limit time-step if particles close to gravity source.
            if (ldt_grav_par.and.(lfirst.and.ldt)) then
              if (lcartesian_coords) then
                vsph=sqrt(fp(k,ivpx)**2+fp(k,ivpy)**2+fp(k,ivpz)**2)
              elseif (lcylindrical_coords) then
                vsph=sqrt(fp(k,ivpx)**2+fp(k,ivpz)**2)
              elseif (lspherical_coords) then
                vsph=abs(fp(k,ivpx))
              endif
              dt1_max(ineargrid(k,1))= &
                  max(dt1_max(ineargrid(k,1)),vsph/rsph/cdtpgrav)
            endif
          enddo
!
        case default
          call fatal_error('dvvp_dt','chosen gravr_profile is not valid!')
!
        endselect
!
      endif
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_nparmin/=0) call max_name(-npar_loc,idiag_nparmin,lneg=.true.)
        if (idiag_nparmax/=0) call max_name(+npar_loc,idiag_nparmax)
        if (idiag_xpm/=0)  call sum_par_name(fp(1:npar_loc,ixp),idiag_xpm)
        if (idiag_ypm/=0)  call sum_par_name(fp(1:npar_loc,iyp),idiag_ypm)
        if (idiag_zpm/=0)  call sum_par_name(fp(1:npar_loc,izp),idiag_zpm)
        if (idiag_xp2m/=0) call sum_par_name(fp(1:npar_loc,ixp)**2,idiag_xp2m)
        if (idiag_yp2m/=0) call sum_par_name(fp(1:npar_loc,iyp)**2,idiag_yp2m)
        if (idiag_zp2m/=0) call sum_par_name(fp(1:npar_loc,izp)**2,idiag_zp2m)
        if (idiag_rpm/=0)  call sum_par_name(sqrt(fp(1:npar_loc,ixp)**2+ &
            fp(1:npar_loc,iyp)**2+fp(1:npar_loc,izp)**2),idiag_rpm)
        if (idiag_rp2m/=0) call sum_par_name(fp(1:npar_loc,ixp)**2+ &
            fp(1:npar_loc,iyp)**2+fp(1:npar_loc,izp)**2,idiag_rp2m)
        if (idiag_vpxm/=0) call sum_par_name(fp(1:npar_loc,ivpx),idiag_vpxm)
        if (idiag_vpym/=0) call sum_par_name(fp(1:npar_loc,ivpy),idiag_vpym)
        if (idiag_vpzm/=0) call sum_par_name(fp(1:npar_loc,ivpz),idiag_vpzm)
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
        if (idiag_lpz2m/=0) call sum_par_name( &
            (fp(1:npar_loc,ixp)*fp(1:npar_loc,ivpy)- &
            fp(1:npar_loc,iyp)*fp(1:npar_loc,ivpx))**2,idiag_lpz2m)
        if (idiag_vpx2m/=0) &
            call sum_par_name(fp(1:npar_loc,ivpx)**2,idiag_vpx2m)
        if (idiag_vpy2m/=0) &
            call sum_par_name(fp(1:npar_loc,ivpy)**2,idiag_vpy2m)
        if (idiag_vpz2m/=0) &
            call sum_par_name(fp(1:npar_loc,ivpz)**2,idiag_vpz2m)
        if (idiag_ekinp/=0) then
          if (lcartesian_coords.and.(all(lequidist))) then
            call sum_par_name(0.5*rhop_swarm*npar_per_cell* &
                 sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2),idiag_ekinp)
          else
            call sum_par_name(0.5*mp_swarm* &
                 sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2),idiag_ekinp)
          endif
        endif
        if (idiag_epotpm/=0) call sum_par_name( &
            -gravr/sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)),idiag_epotpm)
        if (idiag_vpmax/=0) call max_par_name( &
            sqrt(sum(fp(1:npar_loc,ivpx:ivpz)**2,2)),idiag_vpmax)
        if (idiag_vpxmax/=0) call max_par_name(fp(1:npar_loc,ivpx),idiag_vpxmax)
        if (idiag_vpymax/=0) call max_par_name(fp(1:npar_loc,ivpy),idiag_vpymax)
        if (idiag_vpzmax/=0) call max_par_name(fp(1:npar_loc,ivpz),idiag_vpzmax)
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
        if (idiag_rhoptilm/=0) then
          do k=1,npar_loc
            if (lparticles_number) np_swarm=fp(k,inpswarm)
            call sum_par_name((/four_pi_rhopmat_over_three* &
                fp(k,iap)**3*np_swarm/),idiag_rhoptilm)
          enddo
        endif
        if (idiag_mpt/=0) then
          do k=1,npar_loc
            if (lparticles_number) np_swarm=fp(k,inpswarm)
            call integrate_par_name( &
                (/four_pi_rhopmat_over_three*fp(k,iap)**3*np_swarm/),idiag_mpt)
          enddo
        endif
        if (idiag_npargone/=0) then
          call count_particles(ipar,npar_found)
          if (idiag_npargone/=0) &
              call save_name(float(npar-npar_found),idiag_npargone)
        endif
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
      real, dimension (mpar_loc,mpvar) :: fp, dfp
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
      use Cparam, only: lparticles_spin
      use Diagnostics
      use EquationOfState, only: cs20, gamma
      use Particles_spin, only: calc_liftforce
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_npm/=0)      call sum_mn_name(p%np,idiag_npm)
        if (idiag_np2m/=0)     call sum_mn_name(p%np**2,idiag_np2m)
        if (idiag_npmax/=0)    call max_mn_name(p%np,idiag_npmax)
        if (idiag_npmin/=0)    call max_mn_name(-p%np,idiag_npmin,lneg=.true.)
        if (idiag_rhopm/=0)    call sum_mn_name(p%rhop,idiag_rhopm)
        if (idiag_rhop2m/=0 )  call sum_mn_name(p%rhop**2,idiag_rhop2m)
        if (idiag_rhoprms/=0)  call sum_mn_name(p%rhop**2,idiag_rhoprms,lsqrt=.true.)
        if (idiag_rhopmax/=0)  call max_mn_name(p%rhop,idiag_rhopmax)
        if (idiag_rhopmin/=0)  call max_mn_name(-p%rhop,idiag_rhopmin,lneg=.true.)
        if (idiag_epspmax/=0)  call max_mn_name(p%epsp,idiag_epspmax)
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
        if (idiag_rhopmr/=0)  call phizsum_mn_name_r(p%rhop,idiag_rhopmr)
      endif
!
      if (l2davgfirst) then
        if (idiag_rhopmphi/=0) call phisum_mn_name_rz(p%rhop,idiag_rhopmphi)
        if (idiag_rhopmxy/=0)  call zsum_mn_name_xy(p%rhop,idiag_rhopmxy)
        if (idiag_rhopmxz/=0)  call ysum_mn_name_xz(p%rhop,idiag_rhopmxz)
      endif
!
    endsubroutine dvvp_dt_pencil
!***********************************************************************
    subroutine dxxp_dt_blocks(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle position in blocks.
!
!  25-apr-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
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
              if (lshear) then
                dt1_advpy=(-qshear*Omega*fp(k,ixp)+abs(fp(k,ivpy)))&
                           *dy1b(iy0,iblock)
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
      use EquationOfState, only: cs20, gamma
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (mxb,myb,mzb) :: dt1_drag, dt1_drag_gas, dt1_drag_dust
      real, dimension (3) :: dragforce, uup
      real :: rho1_point, tausp1_par, rhop_swarm_pt
      real :: weight, weight_x, weight_y, weight_z
      integer :: k, l, ix0, iy0, iz0, ib, iblock, lb, mb, nb
      integer :: ixx, iyy, izz, ixx0, iyy0, izz0, ixx1, iyy1, izz1
      logical :: lnbody
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
!  Exclude the massive nbody particles from the drag calculations
!
              lnbody=(lparticles_nbody.and.any(ipar(k)==ipar_nbody))
              if (.not.lnbody) then
                ix0=ineargrid(k,1)
                iy0=ineargrid(k,2)
                iz0=ineargrid(k,3)
!
!  The interpolated gas velocity must be calculated here.
!
                if (lhydro) then
                  if (lparticlemesh_cic) then
                    call interpolate_linear(f,iux,iuz, &
                        fp(k,ixp:izp),uup,ineargrid(k,:),iblock,ipar(k))
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
                call get_frictiontime(f,fp,ineargrid,k,tausp1_par,iblock)
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
                      if (nxgrid/=1) weight=weight* &
                          ( 1.0-abs(fp(k,ixp)-xb(ixx,ib))*dx1b(ixx,ib) )
                      if (nygrid/=1) weight=weight* &
                          ( 1.0-abs(fp(k,iyp)-yb(iyy,ib))*dy1b(iyy,ib) )
                      if (nzgrid/=1) weight=weight* &
                          ( 1.0-abs(fp(k,izp)-zb(izz,ib))*dz1b(izz,ib) )
                      if (ldensity_nolog) then
                        rho1_point=1/fb(ixx,iyy,izz,ilnrho,ib)
                      else
                        rho1_point=exp(-fb(ixx,iyy,izz,ilnrho,ib))
                      endif
!  Add friction force to grid point.
                      call get_rhopswarm(mp_swarm,ixx,iyy,izz,ib,rhop_swarm_pt)
                      dfb(ixx,iyy,izz,iux:iuz,ib)=dfb(ixx,iyy,izz,iux:iuz,ib)- &
                          rhop_swarm_pt*rho1_point*dragforce*weight
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
                        if (nxgrid/=1) &
                             weight_x=0.75-((fp(k,ixp)-xb(ixx,ib))*dx1b(ixx,ib))**2
                      endif
                      if ( ((iyy-iy0)==-1) .or. ((iyy-iy0)==+1) ) then
                        weight_y=1.125-1.5* abs(fp(k,iyp)-yb(iyy,ib))*dy1b(iyy,ib) + &
                                       0.5*(abs(fp(k,iyp)-yb(iyy,ib))*dy1b(iyy,ib))**2
                      else
                        if (nygrid/=1) &
                             weight_y=0.75-((fp(k,iyp)-yb(iyy,ib))*dy1b(iyy,ib))**2
                      endif
                      if ( ((izz-iz0)==-1) .or. ((izz-iz0)==+1) ) then
                        weight_z=1.125-1.5* abs(fp(k,izp)-zb(izz,ib))*dz1b(izz,ib) + &
                                       0.5*(abs(fp(k,izp)-zb(izz,ib))*dz1b(izz,ib))**2
                      else
                        if (nzgrid/=1) &
                             weight_z=0.75-((fp(k,izp)-zb(izz,ib))*dz1b(izz,ib))**2
                      endif
!
                      weight=1.0
!
                      if (nxgrid/=1) weight=weight*weight_x
                      if (nygrid/=1) weight=weight*weight_y
                      if (nzgrid/=1) weight=weight*weight_z
                      if (ldensity_nolog) then
                        rho1_point=1/fb(ixx,iyy,izz,irho,ib)
                      else
                        rho1_point=exp(-fb(ixx,iyy,izz,ilnrho,ib))
                      endif
!  Add friction force to grid point.
                      call get_rhopswarm(mp_swarm,ixx,iyy,izz,ib,rhop_swarm_pt)
                      dfb(ixx,iyy,izz,iux:iuz,ib)=dfb(ixx,iyy,izz,iux:iuz,ib) -&
                          rhop_swarm_pt*rho1_point*dragforce*weight
                    enddo; enddo; enddo
                  else
!
!  Nearest Grid Point (NGP) scheme.
!
                    if (ldensity_nolog) then
                      rho1_point=1/fb(ix0,iy0,iz0,ilnrho,ib)
                    else
                      rho1_point=exp(-fb(ix0,iy0,iz0,ilnrho,ib))
                    endif
                    !WL: Why is this l being defined?
                    l=ineargrid(k,1)
                    call get_rhopswarm(mp_swarm,ix0,iy0,iz0,ib,rhop_swarm_pt)
                    dfb(ix0,iy0,iz0,iux:iuz,ib) = dfb(ix0,iy0,iz0,iux:iuz,ib) -&
                        rhop_swarm_pt*rho1_point*dragforce
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
                    if (ldensity_nolog) then
                      rho1_point=1/fb(ix0,iy0,iz0,ilnrho,ib)
                    else
                      rho1_point=exp(-fb(ix0,iy0,iz0,ilnrho,ib))
                    endif
                    if (fb(ix0,iy0,iz0,inp,ib)/=0.0) &
                        dt1_drag_gas(ix0,iy0,iz0)=dt1_drag_gas(ix0,iy0,iz0)+ &
                        (fb(ix0,iy0,iz0,irhop,ib)*rho1_point)/ &
                        fb(ix0,iy0,iz0,inp,ib)*tausp1_par
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
              if (fb(lb,mb,nb,inp,iblock)/=0.0) then
                dt1_max(lb-nghostb)=max(dt1_max(lb-nghostb),dt1_drag(lb,mb,nb))
              endif
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
        if (idiag_nparbmax/=0) &
            call max_name(maxval(npar_iblock(0:nblock_loc-1)),idiag_nparbmax)
      endif
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
!
    endsubroutine dvvp_dt_blocks
!***********************************************************************
    subroutine remove_particles_sink(f,fp,dfp,ineargrid)
!
!  Subroutine for taking particles out of the simulation due to their proximity
!  to a sink particle or sink point.
!
!  25-sep-08/anders: coded
!
      use Solid_Cells
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mpvar) :: fp, dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      real :: rp
      integer :: k
!
      if (lsinkpoint) then
        k=1
        do while (k<=npar_loc)
          rp=sqrt((fp(k,ixp)-xsinkpoint)**2+(fp(k,iyp)-ysinkpoint)**2+ &
             (fp(k,izp)-zsinkpoint)**2)
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
    endsubroutine remove_particles_sink
!***********************************************************************
    subroutine create_sink_particles(f,fp,dfp,ineargrid)
!
!  Subroutine for creating new sink particles or sink points.
!
!  Just a dummy routine for now.
!
!  25-sep-08/anders: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mpvar) :: fp, dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine create_sink_particles
!***********************************************************************
    subroutine get_frictiontime(f,fp,ineargrid,k,tausp1_par,iblock,nochange_opt)
!
!  Calculate the friction time.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer :: k
      real :: tausp1_par
      integer :: iblock
      logical, optional :: nochange_opt
!
      real :: tmp, epsp
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
        if (npar_species>1) then
          jspec=npar_species*(ipar(k)-1)/npar+1
          tmp=tausp1_species(jspec)
        else
          tmp=tausp1
        endif
        tausp1_par=tmp
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
          if (ldensity_nolog) then
            epsp=fb(ix0,iy0,iz0,irhop,iblock)/fb(ix0,iy0,iz0,irho,iblock)
          else
            epsp=fb(ix0,iy0,iz0,irhop,iblock)/exp(fb(ix0,iy0,iz0,ilnrho,iblock))
          endif
          if (epsp>epsp_friction_increase) &
              tausp1_par=tausp1_par/(epsp/epsp_friction_increase)
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
    subroutine read_particles_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_init_pars
!***********************************************************************
    subroutine write_particles_init_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_init_pars)
!
    endsubroutine write_particles_init_pars
!***********************************************************************
    subroutine read_particles_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_run_pars
!***********************************************************************
    subroutine write_particles_run_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_run_pars)
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
      if (lwr) then
        write(3,*) 'ixp=', ixp
        write(3,*) 'iyp=', iyp
        write(3,*) 'izp=', izp
        write(3,*) 'ivpx=', ivpx
        write(3,*) 'ivpy=', ivpy
        write(3,*) 'ivpz=', ivpz
        write(3,*) 'inp=', inp
        write(3,*) 'irhop=', irhop
      endif
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
        idiag_rhoptilm=0; idiag_dtdragp=0; idiag_dedragp=0
        idiag_rhopm=0; idiag_rhoprms=0; idiag_rhop2m=0; idiag_rhopmax=0
        idiag_rhopmin=0; idiag_decollp=0; idiag_rhopmphi=0
        idiag_epspmin=0; idiag_epspmax=0
        idiag_nparmin=0; idiag_nparmax=0; idiag_nmigmax=0; idiag_mpt=0
        idiag_npmx=0; idiag_npmy=0; idiag_npmz=0; idiag_epotpm=0
        idiag_rhopmx=0; idiag_rhopmy=0; idiag_rhopmz=0
        idiag_epspmx=0; idiag_epspmy=0; idiag_epspmz=0
        idiag_rhopmxy=0; idiag_rhopmxz=0; idiag_rhopmr=0
        idiag_dvpx2m=0; idiag_dvpy2m=0; idiag_dvpz2m=0
        idiag_dvpmax=0; idiag_dvpm=0; idiag_nparbmax=0
        idiag_eccpxm=0; idiag_eccpym=0; idiag_eccpzm=0
        idiag_eccpx2m=0; idiag_eccpy2m=0; idiag_eccpz2m=0
        idiag_npargone=0
      endif
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot .and. ip<14) print*,'rprint_particles: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'nparmin',idiag_nparmin)
        call parse_name(iname,cname(iname),cform(iname),'nparmax',idiag_nparmax)
        call parse_name(iname,cname(iname),cform(iname),'nparbmax',idiag_nparbmax)
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
        call parse_name(iname,cname(iname),cform(iname), &
            'rhoptilm',idiag_rhoptilm)
        call parse_name(iname,cname(iname),cform(iname), &
            'dedragp',idiag_dedragp)
        call parse_name(iname,cname(iname),cform(iname), &
            'decollp',idiag_decollp)
        call parse_name(iname,cname(iname),cform(iname), &
            'epotpm',idiag_epotpm)
        call parse_name(iname,cname(iname),cform(iname), &
            'npargone',idiag_npargone)
      enddo
!
!  Check for those quantities for which we want x-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'npmx',idiag_npmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'rhopmx',idiag_rhopmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'epspmx',idiag_epspmx)
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
endmodule Particles
