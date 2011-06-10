! $Id$
!
!  This modules takes care of instantaneous collisions between
!  superparticles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_collisions = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Particles_collisions
!
  use Cparam
  use Cdata
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
  use Particles_sub
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'particles_collisions.h'
!
  real, pointer, dimension (:) :: tausp_species, tausp1_species
  real, pointer :: gravr
  real :: lambda_mfp_single=1.0, coeff_restitution=1.0
  real :: energy_gain_inelastic=0.0
  integer :: ncoll_max_par=-1, npart_max_par=-1, npart_max_par_2=0
  logical :: lcollision_random_angle=.false., lcollision_big_ball=.false.
  logical :: lshear_in_vp=.false., lkeplerian_flat=.false.
  logical :: ltauc_from_tauf=.false., lapproaching_collisions=.false.
  logical :: lstop_at_first_collision=.false.
  character (len=labellen) :: icoll='big-ball'
!
  integer :: idiag_ncollpm=0, idiag_npartpm=0, idiag_decollpm=0
!
  namelist /particles_coll_run_pars/ &
      lambda_mfp_single, coeff_restitution, icoll, lshear_in_vp, &
      ncoll_max_par, npart_max_par, lkeplerian_flat, ltauc_from_tauf, &
      lapproaching_collisions, lstop_at_first_collision
!
  contains
!***********************************************************************
    subroutine initialize_particles_collisions(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  07-oct-08/anders: coded
!
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, intent(in) :: lstarting
!
      select case (icoll)
      case ('random-angle'); lcollision_random_angle=.true.
      case ('big-ball');     lcollision_big_ball=.true.
      case default
        if (lroot) print*, 'No such value for icoll: ', trim(icoll)
        call fatal_error('initialize_particles_collisions','')
      endselect
!
!  Allocate neighbour array necessary for identifying collisions.
!
      if (.not.allocated(kneighbour)) allocate(kneighbour(mpar_loc))
!
      if (lparticles_blocks) then
        if (npart_max_par/=-1 .and. (.not.lrandom_particle_blocks)) then
          if (lroot) then
            print*, 'initialize_particles_collisions: '// &
                'with npart_max_par/=-1 one should set'
            print*, '    lrandom_particle_blocks=T in particles_run_pars'
            print*, '    to avoid artificial sub domains in the collisions'
          endif
          call fatal_error('initialize_particles_collisions','')
        endif
      else
        if (npart_max_par/=-1 .and. (.not.lrandom_particle_pencils)) then
          if (lroot) then
            print*, 'initialize_particles_collisions: '// &
                'with npart_max_par/=-1 one should set'
            print*, '    lrandom_particle_pencils=T in particles_run_pars'
            print*, '    to avoid artificial sub domains in the collisions'
          endif
          call fatal_error('initialize_particles_collisions','')
        endif
      endif
!
!  The maximum number of collision partners must be an even number.
!
      if (npart_max_par/=-1) then
        if (mod(npart_max_par,2)/=0) &
            call fatal_error('initialize_particles_collisions', &
            'npart_par_max must be an even number')
        npart_max_par_2=npart_max_par/2
      endif
!
!  Get radial gravity from gravity module.
!
      if (lkeplerian_flat) call get_shared_variable('gravr',gravr)
!
!  Get friction time from dust particle module.
!
      if (ltauc_from_tauf) then
        call get_shared_variable( 'tausp_species', tausp_species)
        call get_shared_variable('tausp1_species',tausp1_species)
      endif
!
!  Friction time must be set when calculating collision time from friction
!  time.
!
      if (ltauc_from_tauf) then
        if (npar_species==1 .and. tausp1_species(1)==0.0) call fatal_error( &
            'initialize_particles_collisions', 'tausp must be set when '// &
            'calculating collision time from friction time')
      endif
!
      if (.not.(lcartesian_coords.and.(all(lequidist)))) call fatal_error( &
           'initialize_particles_collisions','collisions only implemented '// &
           'for Cartesian equidistant grids.')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_collisions
!***********************************************************************
    subroutine particles_collisions_timestep(fp,ineargrid)
!
!  Time-step contribution from particle collisions.
!
!  30-nov-10/anders: dummy
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine particles_collisions_timestep
!***********************************************************************
    subroutine particles_collisions_pencils(fp,ineargrid)
!
!  Calculate collisions between superparticles by comparing the collision
!  time-scale to the time-step. A random number is used to determine
!  whether two superparticles collide in this time-step.
!
!  Collisions change the velocity vectors of the colliding particles
!  instantaneously.
!
!  23-mar-09/anders: coded
!
      use Diagnostics
      use EquationOfState, only: cs0, rho0
      use General
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: xpj, xpk, vpj, vpk
      real :: deltavjk, tau_coll1, prob, r
      real :: espec, asemi, omega_orbit
      real :: tausp_j, tausp_k, tausp1_j, tausp1_k
      integer, dimension (nx) :: np_pencil
      integer :: l, j, k, npart_point, ncoll, ncoll_par, npart_par
      integer :: jspec, kspec
!
      intent (in) :: ineargrid
      intent (inout) :: fp
!
!  Reset collision counter.
!
      ncoll=0
!
!  Pencil loop.
!
      do imn=1,ny*nz
!
!  Create list of shepherd and neighbour particles for each grid cell in the
!  current pencil.
!
!  Note: with npart_max_par>0, it is safest to also set
!  lrandom_particle_pencils=T. Otherwise there is a risk that a particle
!  always interacts with the same subset of other particles.
!
        call shepherd_neighbour_pencil(fp,ineargrid,kshepherd,kneighbour)
!
!  Calculate number of particles per grid point. This is only needed in order
!  to limit the number of collision partners. Note that f(:,:,:,inp) is not
!  up-to-date at this time.
!
        if (npart_max_par/=-1) then
          np_pencil=0
          if (npar_imn(imn)/=0) then
            do k=k1_imn(imn),k2_imn(imn)
              np_pencil(ineargrid(k,1)-nghost)= &
                  np_pencil(ineargrid(k,1)-nghost)+1
            enddo
          endif
        endif
!
        do l=l1,l2
          k=kshepherd(l-nghost)
          if (k>0) then
            if (npart_max_par/=-1) npart_point=np_pencil(l-nghost)-1
            do while (k/=0)
              j=k
              if (ltauc_from_tauf) then
                kspec=npar_species*(ipar(k)-1)/npar+1
                tausp_k=tausp_species(kspec)
                tausp1_k=tausp1_species(kspec)
              endif
              npart_par=0
              ncoll_par=0
              do while (.true.)
                j=kneighbour(j)
                if (j==0) then
                  if (npart_max_par/=-1 .and. npart_max_par<npart_point) then
                    j=kshepherd(l-nghost)
                  else
                    exit
                  endif
                endif
                if (ltauc_from_tauf) then
                  jspec=npar_species*(ipar(j)-1)/npar+1
                  tausp_j=tausp_species(jspec)
                  tausp1_j=tausp1_species(jspec)
                endif
!
!  Calculate the relative speed of particles j and k.
!
                xpk=fp(k,ixp:izp)
                vpk=fp(k,ivpx:ivpz)
                if (lshear .and. lshear_in_vp) vpk(2)=vpk(2)-qshear*Omega*xpk(1)
                xpj=fp(j,ixp:izp)
                vpj=fp(j,ivpx:ivpz)
                if (lshear .and. lshear_in_vp) vpj(2)=vpj(2)-qshear*Omega*xpj(1)
!
!  Only consider collisions between particles approaching each other?
!
                if ((.not.lapproaching_collisions) .or. &
                    sum((vpk-vpj)*(xpk-xpj))<0.0) then
!
!  For Keplerian particle discs, where the scale height of the particles is
!  given by their velocity dispersion Hp~vrms/OmegaK, we can use the 2-D
!  approach suggested by Lithwick & Chiang (2007). The collision time scale is
!
!    tcol=1/(n*sigma*vrms)=lambda/vrms
!
!  The particle density is n~Sigma/Hp, giving
!
!    tcol=1/(Sigma*sigma*OmegaK)=lambda0/(Omega*dx)
!
!  Here we used that lambda0=1/(npswarm*sigma) has been calculated as if
!  the particle scale height was dx.
!
                  if (lkeplerian_flat) then
                    espec=sum(vpk**2)/2-gravr/sqrt(sum(xpk**2))
                    asemi=-gravr/(2*espec)
                    omega_orbit=sqrt(gravr/asemi**3)
                    deltavjk=omega_orbit*dx
                  else
                    deltavjk=sqrt(sum((vpk-vpj)**2))
                  endif
!
!  The time-scale for collisions between a representative particle from
!  superparticle k and the particle swarm in superparticle j is
!
!    tau_coll = 1/(n*sigma*dv) = lambda/dv
!
!  where lambda is the mean free path of a particle relative to a single
!  superparticle and sigma is the collisional cross section.
!
!  We can further write the collision time of a representative particle from
!  swarm k colliding with the particles of swarm j in terms of their friction
!  times.
!
!    tau_coll_k = 1/(nj*sigma_jk*dv_jk) = mj/(rhoj*sigma_jk*dv_jk)
!               = (4/3)*pi*rhopmat*aj^3/(rhoj*pi*(aj+ak)^2*dv_jk)
!               = (4/3)*rhopmat*aj^3/(rhoj*(aj+ak)^2*dv_jk)
!               = (4/3)*tau_fric_j*(rhog/rhoj)*(cs/dv_jk)*
!                     tau_fric_j^2/(tau_fric_j+tau_fric_k)^2
!
                  if (ltauc_from_tauf) then
                    if (lparticles_radius.or.lparticles_number) then
                      if (lroot) print*, 'particles_collisions_pencils: ', &
                          'not implemented for variable particle radius '// &
                          'or particle number'
                      call fatal_error('particles_collisions_pencils','')
                    endif
                    if (npar_species>1) then
                      tau_coll1=0.75*min(tausp1_j,tausp1_k)*deltavjk/cs0* &
                          rhop_swarm/rho0*(tausp_j+tausp_k)**2* &
                          min(tausp1_j,tausp1_k)**2
                    else
                      tau_coll1=3*tausp1_j*deltavjk/cs0*rhop_swarm/rho0
                    endif
                  else
!
!  If the mean free path is a user supplied constant, then we can readily
!  calculate the collision time-scale.
!
                    tau_coll1=deltavjk/lambda_mfp_single
                  endif
!
!  Increase collision rate artificially for fewer collisions.
!
                  if (npart_max_par/=-1 .and. npart_max_par<npart_point) then
                    tau_coll1=tau_coll1*npart_point/npart_max_par
                  endif
!
                  if (tau_coll1/=0.0) then
!
!  The probability for a collision in this time-step is dt/tau_coll.
!
                    prob=dt*tau_coll1
                    call random_number_wrapper(r)
                    if (r<=prob) then
                      call particle_collision(xpj,xpk,vpj,vpk,j,k)
                      if (lshear .and. lshear_in_vp) then
                        vpk(2)=vpk(2)+qshear*Omega*xpk(1)
                        vpj(2)=vpj(2)+qshear*Omega*xpj(1)
                      endif
                      fp(k,ivpx:ivpz)=vpk
                      fp(j,ivpx:ivpz)=vpj
                      ncoll=ncoll+1
                      ncoll_par=ncoll_par+1
                    endif
                  endif
                endif
                npart_par=npart_par+1
                if (ncoll_max_par/=-1 .and. ncoll_par==ncoll_max_par) exit
                if (npart_max_par/=-1 .and. npart_par==npart_max_par_2) exit
              enddo
              k=kneighbour(k)
!
!  Collision diagnostics. Since this subroutine is called in the last sub-
!  time-step, we can not use ldiagnos. Therefore we calculate collision
!  diagnostics in the preceding time-step. This has the side effect that
!  collision diagnostics are not particle normalized in it==1 or if it1==1.
!
              if (it==1 .or. mod(it,it1)==0) then
                if (idiag_ncollpm/=0) &
                    call sum_par_name((/float(ncoll_par)/),idiag_ncollpm)
                if (idiag_npartpm/=0) &
                    call sum_par_name((/float(npart_par)/),idiag_npartpm)
                if (idiag_decollpm/=0) &
                    call save_name(energy_gain_inelastic/npar,idiag_decollpm)
              endif
!
            enddo
          endif
        enddo
      enddo
!
!  We need to register the diagnostic type, even if there are no particles
!  at the local processor. In the same spirit we calculate diagnostics even
!  for it==1 on processors that do have particles (above). Otherwise a
!  processor starting with N>0 particles, but arriving at mod(it,it1)==0 with
!  zero particles, will not register collision diagnostics properly.
!
      if (it==1) then
        if (npar_loc==0) then
          if (idiag_ncollpm/=0) &
              call sum_par_name(fp(1:npar_loc,ixp),idiag_ncollpm)
          if (idiag_npartpm/=0) &
              call sum_par_name(fp(1:npar_loc,ixp),idiag_npartpm)
          if (idiag_decollpm/=0) call save_name(0.0,idiag_decollpm)
        endif
      endif
!
    endsubroutine particles_collisions_pencils
!***********************************************************************
    subroutine particles_collisions_blocks(fp,ineargrid)
!
!  Calculate collisions between superparticles by comparing the collision
!  time-scale to the time-step. A random number is used to determine
!  whether two superparticles collide in this time-step.
!
!  Collisions change the velocity vectors of the colliding particles
!  instantaneously.
!
!  23-mar-09/anders: coded
!
      use Diagnostics
      use EquationOfState, only: cs0, rho0
      use General
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: xpj, xpk, vpj, vpk
      real :: deltavjk, tau_coll1, prob, r
      real :: espec, asemi, omega_orbit
      real :: tausp_j, tausp_k, tausp1_j, tausp1_k
      integer, dimension (nxb,nyb,nzb) :: kshepherdb, np_block
      integer :: ix, iy, iz, iblock, ix0, iy0, iz0
      integer :: j, k, npart_point, ncoll, ncoll_par, npart_par
      integer :: jspec, kspec
!
      intent (in) :: ineargrid
      intent (inout) :: fp
!
!  Reset collision counter.
!
      ncoll=0
!
!  Block loop.
!
      do iblock=0,nblock_loc-1
!
!  Create list of shepherd and neighbour particles for each grid cell in the
!  current block.
!
!  Note: with npart_max_par>0, it is safest to also set
!  lrandom_particle_pencils=T. Otherwise there is a risk that a particle
!  always interacts with the same subset of other particles.
!
        call shepherd_neighbour_block(fp,ineargrid,kshepherdb,kneighbour,iblock)
!
!  Calculate number of particles per grid point. This is only needed in order
!  to limit the number of collision partners. Note that f(:,:,:,inp) is not
!  up-to-date at this time.
!
        if (npart_max_par/=-1) then
          np_block=0
          if (npar_iblock(iblock)/=0) then
            do k=k1_iblock(iblock),k2_iblock(iblock)
              ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
              np_block(ix0-nghostb,iy0-nghostb,iz0-nghostb)= &
                  np_block(ix0-nghostb,iy0-nghostb,iz0-nghostb)+1
            enddo
          endif
        endif
!
        do iz=n1b,n2b; do iy=m1b,m2b; do ix=l1b,l2b
          k=kshepherdb(ix-nghostb,iy-nghostb,iz-nghostb)
          if (k>0) then
            if (npart_max_par/=-1) &
                npart_point=np_block(ix-nghostb,iy-nghostb,iz-nghostb)-1
            do while (k/=0)
              j=k
              if (ltauc_from_tauf) then
                kspec=npar_species*(ipar(k)-1)/npar+1
                tausp_k=tausp_species(kspec)
                tausp1_k=tausp1_species(kspec)
              endif
              npart_par=0
              ncoll_par=0
              do while (.true.)
                j=kneighbour(j)
                if (j==0) then
                  if (npart_max_par/=-1 .and. npart_max_par<npart_point) then
                    j=kshepherdb(ix-nghostb,iy-nghostb,iz-nghostb)
                  else
                    exit
                  endif
                endif
                if (ltauc_from_tauf) then
                  jspec=npar_species*(ipar(j)-1)/npar+1
                  tausp_j=tausp_species(jspec)
                  tausp1_j=tausp1_species(jspec)
                endif
!
!  Calculate the relative speed of particles j and k.
!
                xpk=fp(k,ixp:izp)
                vpk=fp(k,ivpx:ivpz)
                if (lshear .and. lshear_in_vp) vpk(2)=vpk(2)-qshear*Omega*xpk(1)
                xpj=fp(j,ixp:izp)
                vpj=fp(j,ivpx:ivpz)
                if (lshear .and. lshear_in_vp) vpj(2)=vpj(2)-qshear*Omega*xpj(1)
!
!  Only consider collisions between particles approaching each other?
!
                if ((.not.lapproaching_collisions) .or. &
                    sum((vpk-vpj)*(xpk-xpj))<0.0) then
!
!  For Keplerian particle discs, where the scale height of the particles is
!  given by their velocity dispersion Hp~vrms/OmegaK, we can use the 2-D
!  approach suggested by Lithwick & Chiang (2007). The collision time scale is
!
!    tcol=1/(n*sigma*vrms)=lambda/vrms
!
!  The particle density is n~Sigma/Hp, giving
!
!    tcol=1/(Sigma*sigma*OmegaK)=lambda0/(Omega*dx)
!
!  Here we used that lambda0=1/(npswarm*sigma) has been calculated as if
!  the particle scale height was dx.
!
                  if (lkeplerian_flat) then
                    espec=sum(vpk**2)/2-gravr/sqrt(sum(xpk**2))
                    asemi=-gravr/(2*espec)
                    omega_orbit=sqrt(gravr/asemi**3)
                    deltavjk=omega_orbit*dx
                  else
                    deltavjk=sqrt(sum((vpk-vpj)**2))
                  endif
!
!  The time-scale for collisions between a representative particle from
!  superparticle k and the particle swarm in superparticle j is
!
!    tau_coll = 1/(n*sigma*dv) = lambda/dv
!
!  where lambda is the mean free path of a particle relative to a single
!  superparticle and sigma is the collisional cross section.
!
!  We can further write the collision time of a representative particle from
!  swarm k colliding with the particles of swarm j in terms of their friction
!  times.
!
!    tau_coll_k = 1/(nj*sigma_jk*dv_jk) = mj/(rhoj*sigma_jk*dv_jk)
!               = (4/3)*pi*rhopmat*aj^3/(rhoj*pi*(aj+ak)^2*dv_jk)
!               = (4/3)*rhopmat*aj^3/(rhoj*(aj+ak)^2*dv_jk)
!               = (4/3)*tau_fric_j*(rhog/rhoj)*(cs/dv_jk)*
!                     tau_fric_j^2/(tau_fric_j+tau_fric_k)^2
!
                  if (ltauc_from_tauf) then
                    if (lparticles_radius.or.lparticles_number) then
                      if (lroot) print*, 'particles_collisions_blocks: ', &
                          'not implemented for variable particle radius '// &
                          'or particle number'
                      call fatal_error('particles_collisions_blocks','')
                    endif
                    if (npar_species>1) then
                      tau_coll1=0.75*min(tausp1_j,tausp1_k)*deltavjk/cs0* &
                          rhop_swarm/rho0*(tausp_j+tausp_k)**2* &
                          min(tausp1_j,tausp1_k)**2
                    else
                      tau_coll1=3*tausp1_j*deltavjk/cs0*rhop_swarm/rho0
                    endif
                  else
!
!  If the mean free path is a user supplied constant, then we can readily
!  calculate the collision time-scale.
!
                    tau_coll1=deltavjk/lambda_mfp_single
                  endif
!
!  Increase collision rate artificially for fewer collisions.
!
                  if (npart_max_par/=-1 .and. npart_max_par<npart_point) then
                    tau_coll1=tau_coll1*npart_point/npart_max_par
                  endif
!
                  if (tau_coll1/=0.0) then
!
!  The probability for a collision in this time-step is dt/tau_coll.
!
                    prob=dt*tau_coll1
                    call random_number_wrapper(r)
                    if (r<=prob) then
                      call particle_collision(xpj,xpk,vpj,vpk,j,k)
                      if (lshear .and. lshear_in_vp) then
                        vpk(2)=vpk(2)+qshear*Omega*xpk(1)
                        vpj(2)=vpj(2)+qshear*Omega*xpj(1)
                      endif
                      fp(k,ivpx:ivpz)=vpk
                      fp(j,ivpx:ivpz)=vpj
                      ncoll=ncoll+1
                      ncoll_par=ncoll_par+1
                    endif
                  endif
                endif
                npart_par=npart_par+1
                if (ncoll_max_par/=-1 .and. ncoll_par==ncoll_max_par) exit
                if (npart_max_par/=-1 .and. npart_par==npart_max_par_2) exit
              enddo
              k=kneighbour(k)
!
!  Collision diagnostics. Since this subroutine is called in the last sub-
!  time-step, we can not use ldiagnos. Therefore we calculate collision
!  diagnostics in the preceding time-step. This has the side effect that
!  collision diagnostics are not particle normalized in it==1 or if it1==1.
!
              if (it==1 .or. mod(it,it1)==0) then
                if (idiag_ncollpm/=0) &
                    call sum_par_name((/float(ncoll_par)/),idiag_ncollpm)
                if (idiag_npartpm/=0) &
                    call sum_par_name((/float(npart_par)/),idiag_npartpm)
                if (idiag_decollpm/=0) &
                    call save_name(energy_gain_inelastic/npar,idiag_decollpm)
              endif
!
            enddo
          endif
        enddo; enddo; enddo
      enddo
!
!  We need to register the diagnostic type, even if there are no particles
!  at the local processor. In the same spirit we calculate diagnostics even
!  for it==1 on processors that do have particles (above). Otherwise a
!  processor starting with N>0 particles, but arriving at mod(it,it1)==0 with
!  zero particles, will not register collision diagnostics properly.
!
      if (it==1) then
        if (npar_loc==0) then
          if (idiag_ncollpm/=0) &
              call sum_par_name(fp(1:npar_loc,ixp),idiag_ncollpm)
          if (idiag_npartpm/=0) &
              call sum_par_name(fp(1:npar_loc,ixp),idiag_npartpm)
          if (idiag_decollpm/=0) call save_name(0.0,idiag_decollpm)
        endif
      endif
!
    endsubroutine particles_collisions_blocks
!***********************************************************************
    subroutine particle_collision(xpj,xpk,vpj,vpk,j,k)
!
!  Calculate collision between two particles.
!
!  13-nov-09/anders: coded
!
      use General
      use Mpicomm
      use Sub
!
      real, dimension (3) :: xpj, xpk, vpj, vpk
      integer :: j,k
!
      real, dimension (3) :: vvcm, vvkcm, vvkcmnew, vvkcm_normal, vvkcm_parall
      real, dimension (3) :: tmp1, tmp2, nvec
      real :: theta_rot, phi_rot
!
      intent (inout) :: xpj, xpk, vpj,vpk
      intent (in) ::j, k
!
      if (ip<=6) then
        print*, 'particle_collision: collision between'// &
            ' superparticles ', ipar(k), ' and ', ipar(j)
        print*, 'particle_collision: xj, xk='
        print*, '  ', xpj, sqrt(sum(xpj**2))
        print*, '  ', xpk, sqrt(sum(xpk**2))
        print*, 'particle_collision: **before**'
        print*, 'particle_collision: vj, vk, vj+vk='
        print*, '  ', vpj, sqrt(sum(vpj**2))
        print*, '  ', vpk, sqrt(sum(vpk**2))
        print*, '  ', vpj+vpk
        print*, 'particle_collision: ej, ek, ej+ek='
        print*, '  ', 0.5*(sum(vpj**2)), 0.5*(sum(vpk**2)), &
            0.5*(sum(vpj**2))+0.5*(sum(vpk**2))
        print*, 'particle_collision: lj, lk, lj+lk='
        call cross(xpj,vpj,tmp1)
        print*, '  ', tmp1
        call cross(xpk,vpk,tmp2)
        print*, '  ', tmp2
        print*, '  ', tmp1+tmp2
      endif
!
!  Choose random unit vector direction in COM frame. Here theta_rot is the
!  pole angle and phi_rot is the azimuthal angle. While we can choose phi_rot
!  randomly, theta_rot must be chosen with care so that equal areas on the
!  unit sphere are equally probable.
!    (see http://mathworld.wolfram.com/SpherePointPicking.html)
!
      if (lcollision_random_angle) then
        call random_number_wrapper(theta_rot)
        theta_rot=acos(2*theta_rot-1)
        call random_number_wrapper(phi_rot)
        phi_rot  =phi_rot*2*pi
        vvcm=0.5*(vpj+vpk)
        vvkcm=vpk-vvcm
        vvkcmnew(1)=+sin(theta_rot)*cos(phi_rot)
        vvkcmnew(2)=+sin(theta_rot)*sin(phi_rot)
        vvkcmnew(3)=+cos(theta_rot)
!
!  Multiply unit vector by length of velocity vector.
!
        vvkcmnew=vvkcmnew* &
            sqrt(vvkcm(1)**2+vvkcm(2)**2+vvkcm(3)**2)
!
!  Inelastic collision leads to energy loss of both particles.
!
        if (coeff_restitution/=1.0) then
          if (energy_gain_inelastic/=impossible) &
              energy_gain_inelastic=energy_gain_inelastic- &
              (1.0-coeff_restitution**2)*sum(vvkcm_normal**2)
          vvkcmnew=vvkcmnew*coeff_restitution
        endif
!
!  Change velocity vectors in normal frame.
!
        vpk=+vvkcmnew+vvcm
        vpj=-vvkcmnew+vvcm
!
      endif
!
!  Alternative method for collisions where each particle is considered to be
!  a big ball stretching to the surface of the other particle. We can then
!  solve the collision problem with perfect conservation of momentum, energy
!  and angular momentum.
!
      if (lcollision_big_ball) then
        nvec=xpj-xpk
        nvec=nvec/sqrt(sum(nvec**2))
        vvcm=0.5*(vpj+vpk)
        vvkcm=vpk-vvcm
        vvkcm_normal=nvec*(sum(vvkcm*nvec))
        vvkcm_parall=vvkcm-vvkcm_normal
        if (ip<=6) then
          print*, 'particle_collision: nvec, vvkcm, vvkcm_normal, vvkcm_parall='
          print*, '  ', nvec
          print*, '  ', vvkcm
          print*, '  ', vvkcm_normal
          print*, '  ', vvkcm_parall
        endif
!
!  Inelastic collision leads to energy loss of both particles.
!
        if (coeff_restitution/=1.0) then
          if (energy_gain_inelastic/=impossible) &
              energy_gain_inelastic=energy_gain_inelastic- &
              (1.0-coeff_restitution**2)*sum(vvkcm_normal**2)
          vvkcm_normal=vvkcm_normal*coeff_restitution
        endif
        vpk=vvcm+vvkcm_parall-vvkcm_normal
        vpj=vvcm-vvkcm_parall+vvkcm_normal
      endif
      if (ip<=6) then
        print*, 'particle_collision: **after**'
        print*, 'particle_collision: vj, vk, vj+vk='
        print*, '  ', vpj, sqrt(sum(vpj**2))
        print*, '  ', vpk, sqrt(sum(vpk**2))
        print*, '  ', vpj+vpk
        print*, 'particle_collision: ej, ek, ej+ek='
        print*, '  ', 0.5*(sum(vpj**2)), 0.5*(sum(vpk**2)), &
            0.5*(sum(vpj**2))+0.5*(sum(vpk**2))
        print*, 'particle_collision: lj, lk, lj+lk='
        call cross(xpj,vpj,tmp1)
        print*, '  ', tmp1
        call cross(xpk,vpk,tmp2)
        print*, '  ', tmp2
        print*, '  ', tmp1+tmp2
      endif
!
!  Stop after one collision (for testing purposes).
!
      if (lstop_at_first_collision) call fatal_error('particle_collision', &
          'stopping after first collision')
!
    endsubroutine particle_collision
!***********************************************************************
    subroutine read_particles_coll_run_pars(unit,iostat)
!
!  Read run parameters from run.in.
!
!  28-mar-09/anders: adapted
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,nml=particles_coll_run_pars,err=99,iostat=iostat)
      else
        read(unit,nml=particles_coll_run_pars,err=99)
      endif
!
99    return
!
    endsubroutine read_particles_coll_run_pars
!***********************************************************************
    subroutine write_particles_coll_run_pars(unit)
!
!  Write run parameters to param.nml.
!
!  28-mar-09/anders: adapted
!
      integer, intent(in) :: unit
!
      write(unit,NML=particles_coll_run_pars)
!
    endsubroutine write_particles_coll_run_pars
!*******************************************************************
    subroutine rprint_particles_collisions(lreset,lwrite)
!
!  Read and register diagnostic parameters.
!
!  28-mar-09/anders: adapted
!
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
!
      if (lreset) then
        idiag_ncollpm=0; idiag_npartpm=0; idiag_decollpm=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ncollpm',idiag_ncollpm)
        call parse_name(iname,cname(iname),cform(iname),'npartpm',idiag_npartpm)
        call parse_name(iname,cname(iname),cform(iname), &
            'decollpm',idiag_decollpm)
      enddo
!
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_particles_collisions
!***********************************************************************
endmodule Particles_collisions
