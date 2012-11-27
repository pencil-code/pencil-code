! $Id$
!
!  This modules takes care of instantaneous coagulation, shattering,
!  erosion, and bouncing of superparticles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_coagulation = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Particles_coagulation
!
  use Cdata
  use Cparam
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
  use Particles_sub
  use General, only: keep_compiler_quiet
!
  implicit none
!
  include 'particles_coagulation.h'
!
  real :: cdtpcoag=0.2, cdtpcoag1=5.0
  real :: kernel_cst=1.0, kernel_lin=1.0, kernel_pro=1.0
  real :: four_pi_rhopmat_over_three2=0.0
  real :: three_over_four_pi_rhopmat=0.0
  real :: GNewton=6.67428e-11, deltav_grav_floor=0.0
  real :: critical_mass_ratio_sticking=1.0
  real :: minimum_particle_mass=0.0, minimum_particle_radius=0.0
  real, pointer :: rhs_poisson_const
  logical :: lcoag_simultaneous=.false., lnoselfcollision=.true.
  logical :: lshear_in_vp=.true.
  logical :: lkernel_test=.false., lconstant_kernel_test=.false.
  logical :: llinear_kernel_test=.false., lproduct_kernel_test=.false.
  logical :: lgravitational_cross_section=.false.
  logical :: lzsomdullemond=.false.     ! use Zsom and Dullemond method
!
  real, dimension(:,:), allocatable :: r_ik_mat, cum_func_sec_ik
  real, dimension(:), allocatable :: r_i_tot, cum_func_first_i
  real :: r_total = 0.0
  real :: delta_r = 0.0
!
  integer :: idiag_ncoagpm=0, idiag_ncoagpartpm=0
!
  namelist /particles_coag_run_pars/ &
      cdtpcoag, lcoag_simultaneous, lshear_in_vp, lconstant_kernel_test, &
      kernel_cst, llinear_kernel_test, kernel_lin, lproduct_kernel_test, &
      kernel_pro, lnoselfcollision, lgravitational_cross_section, &
      GNewton, deltav_grav_floor, critical_mass_ratio_sticking, &
      minimum_particle_mass, minimum_particle_radius, lzsomdullemond
!
  contains
!***********************************************************************
    subroutine initialize_particles_coag(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-10/anders: coded
!  15-nov-12/KWJ: modified
!
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, intent(in) :: lstarting
!
!  Fatal error if Particle_radius module not used.
!
      if (.not.lparticles_radius) &
          call fatal_error('initialize_particles_coag', &
          'must use Particles_radius module for coagulation')
!
!  Allocate neighbour array necessary for identifying collisions.
!
      if (.not.allocated(kneighbour)) allocate(kneighbour(mpar_loc))
!
!  Get the gravity constant from the Selfgravity module.
!
      if (lselfgravity) then
        call get_shared_variable('rhs_poisson_const',rhs_poisson_const)
        GNewton=rhs_poisson_const/(4*pi)
      endif
!
!  Precalculate inverse of coagulation time-step parameter.
!
      cdtpcoag1=1/cdtpcoag
!
!  Short hand for any kernel test.
!
      lkernel_test=lconstant_kernel_test.or.llinear_kernel_test.or. &
          lproduct_kernel_test
!
!  Calculate squared and inverse volume factor.
!
      four_pi_rhopmat_over_three2=four_pi_rhopmat_over_three**2
      three_over_four_pi_rhopmat=1/four_pi_rhopmat_over_three
!
!  Define the minimum particle mass or radius that will be allowed. Small
!  particles may form by erosion and fragmentation.
!
      if (minimum_particle_mass/=0.0 .and. minimum_particle_radius/=0.0) then
        call fatal_error('initialize_particles_coag', &
            'not allowed to set both minimum_particle_mass and '// &
            'minimum_particle_radius')
      endif
      if (minimum_particle_mass/=0.0) then
        minimum_particle_radius = &
            (minimum_particle_mass*three_over_four_pi_rhopmat)**(1.0/3.0)
      elseif (minimum_particle_radius/=0.0) then
        minimum_particle_mass = &
            four_pi_rhopmat_over_three*minimum_particle_radius**3
      endif
!
      if (.not.(lcartesian_coords.and.(all(lequidist)))) call fatal_error( &
           'initialize_particles_coagulation', 'Coagulation only '// &
           'implemented for Cartesian equidistant grids.')
!
!  If using the Zsom-Dullemond Monte Carlo method
!
      if(lzsomdullemond.and.(.not.allocated(r_ik_mat))) then
        allocate(r_ik_mat(mpar_loc,mpar_loc))
        allocate(r_i_tot(mpar_loc))
        allocate(cum_func_first_i(mpar_loc))
        allocate(cum_func_sec_ik(mpar_loc,mpar_loc))
      end if
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_coag
!***********************************************************************
    subroutine particles_coagulation_timestep(fp,ineargrid)
!
!  Time-step contribution from particle coagulation.
!
!  30-nov-10/anders: coded
!  15-nov-12/KWJ: modified
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: xpj, xpk, vpj, vpk
      real :: deltavjk, dt1_coag_par, kernel
      real :: npswarmj, npswarmk
      integer :: j, k, l
!
!  If using the Zsom and Dullemond Monte Carlo method
!
      if(lzsomdullemond) then
        call particles_coag_timestep_zd(fp)
        return
      end if
!
      if (lfirst.and.ldt) then
!
!  Create list of shepherd and neighbour particles for each grid cell in the
!  current pencil.
!
        call shepherd_neighbour_pencil(fp,ineargrid,kshepherd,kneighbour)
!
!  Calculate coagulation time-step.
!
        do l=l1,l2
!
!  Start with shepherd particle.
!
          k=kshepherd(l-nghost)
          if (k>0) then
            do while (k/=0)
              if (lparticles_number) then
                npswarmk=fp(k,inpswarm)
              else
                npswarmk=rhop_swarm/(four_pi_rhopmat_over_three*fp(k,iap)**3)
              endif
              dt1_coag_par=0.0
              j=kshepherd(l-nghost)
!
!  Move through neighbours.
!
              do while (.true.)
                if (lparticles_number) then
                  npswarmj=fp(j,inpswarm)
                else
                  npswarmj=rhop_swarm/(four_pi_rhopmat_over_three*fp(j,iap)**3)
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
!  Special treatment for kernel tets.
!
                if (lkernel_test) then
                  if (lconstant_kernel_test) then
                    kernel=kernel_cst
                  elseif (llinear_kernel_test) then
                    kernel=kernel_lin* &
                       four_pi_rhopmat_over_three*(fp(j,iap)**3+fp(k,iap)**3)
                  elseif (lproduct_kernel_test) then
                    kernel=kernel_pro* &
                       four_pi_rhopmat_over_three2*fp(j,iap)**3*fp(k,iap)**3
                  endif
                  if (j==k .and. lnoselfcollision) kernel=0.0
                  dt1_coag_par=dt1_coag_par+kernel*min(npswarmj,npswarmk)
                else
!
!  Time-step for physical kernel.
!
                  if (sum((vpk-vpj)*(xpk-xpj))<0.0) then
                    deltavjk=sqrt(sum((vpk-vpj)**2))
                    dt1_coag_par=dt1_coag_par+ &
                        pi*(fp(k,iap)+fp(k,iap))**2*deltavjk* &
                        min(npswarmj,npswarmk)
                  endif
!
                endif
                j=kneighbour(j)
                if (j==0) exit
              enddo
!
!  Put particle's time-step into inverse time-step array.
!
              dt1_max(l-nghost)=max(dt1_max(l-nghost),dt1_coag_par*cdtpcoag1)
!
!  Move to next particle in the grid cell.
!
              k=kneighbour(k)
!
            enddo
          endif
        enddo
      endif
!
    endsubroutine particles_coagulation_timestep
!***********************************************************************
    subroutine particles_coagulation_pencils(fp,ineargrid)
!
!  Calculate outcome of superparticle collisions by comparing the collision
!  time-scale to the time-step. A random number is used to determine
!  whether two superparticles collide in this time-step.
!
!  Collisions lead to coagulation, shattering, erosion, or bouncing.
!
!  24-nov-10/anders: coded
!  15-nov-12/KWJ: modified
!
      use Diagnostics
      use General, only: random_number_wrapper
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: xpj, xpk, vpj, vpk
      real :: lambda_mfp1, deltavjk, tau_coll1, prob, r, kernel
      real :: npswarmj, npswarmk
      integer :: l, j, k, ncoll, ncoll_par, npart_par
!
      intent (in) :: ineargrid
      intent (inout) :: fp
!
!  If using the Zsom & Dullemond Monte Carlo method (KWJ)
!
      if(lzsomdullemond) then
        call particles_coag_outcome_zd(fp)
        return
      end if
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
        call shepherd_neighbour_pencil(fp,ineargrid,kshepherd,kneighbour)
!
        do l=l1,l2
!
!  Start with shepherd particle.
!
          k=kshepherd(l-nghost)
          if (k>0) then
            do while (k/=0)
              if (lcoag_simultaneous) then
                j=k
              else
                j=kshepherd(l-nghost)
              endif
              npart_par=0
              ncoll_par=0
!
              if (lparticles_number) then
                npswarmk=fp(k,inpswarm)
              else
                npswarmk=rhop_swarm/(four_pi_rhopmat_over_three*fp(k,iap)**3)
              endif
!
!  Move through neighbours, excluding particles tagged with negative radius
!  (they have coagulated with a sink particle and will be removed later).
!
              do while (fp(k,iap)>=0.0)
                if (lcoag_simultaneous) then
                  j=kneighbour(j)
                  if (j==0) exit
!
!  Do not attempt to collide with particles tagged for removal.
!
                  if (lparticles_number) then
                    do while (fp(j,iap)<0.0)
                      j=kneighbour(j)
                      if (j==0) exit
                    enddo
                    if (j==0) exit
                  endif
                endif
!
!  Particle number density either from f array or from rhop_swarm.
!
                if (lparticles_number) then
                  npswarmj=fp(j,inpswarm)
                else
                  npswarmj=rhop_swarm/(four_pi_rhopmat_over_three*fp(j,iap)**3)
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
!  Only consider collisions between particles approaching each other.
!
                if ((sum((vpk-vpj)*(xpk-xpj))<0.0).or.lkernel_test) then
!
!  Relative particle speed.
!
                  deltavjk=sqrt(sum((vpk-vpj)**2))
!
!  The time-scale for collisions between a representative particle from
!  superparticle k and the particle swarm in superparticle j is
!
!    tau_coll = 1/(n_j*sigma_jk*deltav_jk)
!
!  where n_j is the number density of swarm j and sigma_jk is the collisional
!  cross section [=pi*(a_j+a_k)^2].
!
!  The collision time-scale above is actually the *mass collision time-scale*,
!  i.e. the time for a particle to collide with its own mass. Since super-
!  particles always represent the same mass, that means that
!
!    - The interaction time-scale for a particle k to collide with a swarm
!      of larger particles j is simply the collision time-scale (when one
!      small particle has collid with a large, then they all have collided).
!
!    - The interaction time-scale for a particle k to collide with a swarm
!      of smaller particles is the collision time-scale times the mass ratio
!      m_k/m_j, so that tau_coll = 1/(n_k*sigma_jk*deltav_jk).
!
!  We use lcoag_simultaneous to specify that two colliding superparticles
!  change internal properties at the same time. The default is that the
!  swarms evolve separately.
!
                  if (lkernel_test) then
                    if (lconstant_kernel_test) then
                      kernel=kernel_cst
                    elseif (llinear_kernel_test) then
                      kernel=kernel_lin* &
                        four_pi_rhopmat_over_three*(fp(j,iap)**3+fp(k,iap)**3)
                    elseif (lproduct_kernel_test) then
                      kernel=kernel_pro* &
                        four_pi_rhopmat_over_three2*fp(j,iap)**3*fp(k,iap)**3
                    endif
                    if (j==k .and. lnoselfcollision) kernel=0.0
                    if (fp(k,iap)<fp(j,iap)) then
                      tau_coll1=kernel*npswarmj
                    else
                      tau_coll1=kernel*npswarmk
                    endif
                  else
                    tau_coll1=deltavjk*pi*(fp(k,iap)+fp(j,iap))**2* &
                        min(npswarmj,npswarmk)
                    if (lgravitational_cross_section) then
                      tau_coll1=tau_coll1*(1.0+2*GNewton* &
                          four_pi_rhopmat_over_three* &
                          (fp(j,iap)**3+fp(k,iap)**3)/((fp(k,iap)+fp(k,iap))* &
                          (deltavjk+deltav_grav_floor)**2))
                    endif
                  endif
!
                  if (tau_coll1/=0.0) then
!
!  The probability for a collision in this time-step is dt/tau_coll.
!
                    prob=dt*tau_coll1
                    call random_number_wrapper(r)
                    if (r<=prob) then
!
                      call coagulation_fragmentation(fp,j,k,deltavjk, &
                          xpj,xpk,vpj,vpk)
!
                      if (lparticles_number) then
                        npswarmk=fp(k,inpswarm)
                      else
                        npswarmk=rhop_swarm/(four_pi_rhopmat_over_three* &
                            fp(k,iap)**3)
                      endif
!
                      ncoll=ncoll+1
                      ncoll_par=ncoll_par+1
!
                    endif
                  endif
                endif
!
!  Move to next particle neighbour.
!
                if (.not.lcoag_simultaneous) then
                  j=kneighbour(j)
                  if (j==0) exit
                endif
!
              enddo
!
!  Collision diagnostics. Since this subroutine is called in the last sub-
!  time-step, we can not use ldiagnos. Therefore we calculate collision
!  diagnostics in the preceding time-step. This has the side effect that
!  collision diagnostics are not particle normalized in it==1 or if it1==1.
!
              if (it==1 .or. mod(it,it1)==0) then
                if (idiag_ncoagpm/=0) &
                    call sum_par_name((/float(ncoll_par)/),idiag_ncoagpm)
                if (idiag_ncoagpartpm/=0) &
                    call sum_par_name((/float(npart_par)/),idiag_ncoagpartpm)
              endif
!
!  Move to next particle in the grid cell.
!
              k=kneighbour(k)
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
          if (idiag_ncoagpm/=0) &
              call sum_par_name(fp(1:npar_loc,ixp),idiag_ncoagpm)
          if (idiag_ncoagpartpm/=0) &
              call sum_par_name(fp(1:npar_loc,ixp),idiag_ncoagpartpm)
        endif
      endif
!
!  Remove the particles that have been tagged for removal. We have to remove
!  them after the coagulation step, as we will otherwise confuse the
!  shepherd-neighbour algorithm.
!
      k=1
      do while (.true.)
        if (fp(k,iap)<0.0) then
          fp(k,iap)=-fp(k,iap)
          call remove_particle(fp,ipar,k)
        else
          k=k+1
          if (k>npar_loc) exit
        endif
      enddo
!
    endsubroutine particles_coagulation_pencils
!***********************************************************************
    subroutine particles_coagulation_blocks(fp,ineargrid)
!
!  Calculate outcome of superparticle collisions by comparing the collision
!  time-scale to the time-step. A random number is used to determine
!  whether two superparticles collide in this time-step.
!
!  Collisions lead to coagulation, shattering, erosion, or bouncing.
!
!  24-nov-10/anders: coded
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call fatal_error('particles_coagulation_blocks','not implemented yet')
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine particles_coagulation_blocks
!***********************************************************************
    subroutine coagulation_fragmentation(fp,j,k,deltavjk,xpj,xpk,vpj,vpk)
!
!  Change the particle size to the new size, but keep the total mass in the
!  particle swarm the same.
!
!  A representative particle k colliding with a swarm of larger particles j
!  simply obtains the new mass m_k -> m_k + m_j.
!
!  A representative particle k colliding with a swarm of smaller particles j
!  obtains the new mass m_k -> 2*m_k.
!
!  24-nov-10/anders: coded
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer :: j, k
      real :: deltavjk
      real, dimension (3) :: xpj, xpk, vpj, vpk
!
      real, dimension (3) :: vvcm, vvkcm, vvkcmnew, vvkcm_normal, vvkcm_parall
      real, dimension (3) :: nvec
      real :: mpsma, mpbig, npsma, npbig, npnew, mpnew, apnew
      real :: rhopsma, rhopbig, apsma, apbig
      real :: coeff_restitution, deltav_recoil, escape_speed
!
      if (lparticles_number) then
!
!  Physical particles in the two swarms collide simultaneously.
!
        if (lcoag_simultaneous) then
!
!  Define mpsma, npsma, rhopsma for the superparticle with smaller particles.
!  Define mpbig, npbig, rhopbig for the superparticle with bigger particles.
!
          if (fp(j,iap)<fp(k,iap)) then
            apsma = fp(j,iap)
            mpsma = four_pi_rhopmat_over_three*fp(j,iap)**3
            npsma = fp(j,inpswarm)
            apbig = fp(k,iap)
            mpbig = four_pi_rhopmat_over_three*fp(k,iap)**3
            npbig = fp(k,inpswarm)
          else
            apsma = fp(k,iap)
            mpsma = four_pi_rhopmat_over_three*fp(k,iap)**3
            npsma = fp(k,inpswarm)
            apbig = fp(j,iap)
            mpbig = four_pi_rhopmat_over_three*fp(j,iap)**3
            npbig = fp(j,inpswarm)
          endif
          rhopsma=mpsma*npsma
          rhopbig=mpbig*npbig
!
!  Calculate the coefficient of restitution. We base in here on the
!  experimental fit of Higa et al. (1986) to ice at 100 K.
!
          coeff_restitution=(deltavjk/1.8)**(-alog10(deltavjk/1.8))
!
!  Particles stick:
!
!    a) when the mass ratio between large and small particles is
!       sufficiently high.
!
!    b) when the recoil velocity is less than the escape speed
!
          nvec=(xpj-xpk)
          nvec=nvec/sqrt(sum(nvec**2))
          vvcm=0.5*(vpj+vpk)
          vvkcm=vpk-vvcm
          vvkcm_normal=nvec*(sum(vvkcm*nvec))
          vvkcm_parall=vvkcm-vvkcm_normal
          deltav_recoil=sqrt(sum((vvkcm_parall - &
              coeff_restitution*vvkcm_normal)**2))
          escape_speed =sqrt(2*GNewton*mpbig/apbig)
!
          if (mpbig/mpsma>critical_mass_ratio_sticking .or. &
              deltav_recoil<escape_speed .or. &
              fp(j,inpswarm)==1/(dx*dy*dz) .or. &
              fp(k,inpswarm)==1/(dx*dy*dz)) then
            mpnew=mpbig+rhopsma/npbig
          else
            if (mpsma>2*minimum_particle_mass) then
              mpnew=0.5*mpsma
            else
              mpnew=minimum_particle_mass
            endif
          endif
!
!  The new radius is defined from the new particle mass, while the new
!  particle number in each superparticle comes from total mass conservation.
!
          apnew=(mpnew*three_over_four_pi_rhopmat)**(1.0/3.0)
          npnew=0.5*(rhopsma+rhopbig)/mpnew
!
!  Turn into sink particle if number of physical particles is less than one.
!
          if (npnew*dx*dy*dz<1.0) then
            if (fp(j,iap)<fp(k,iap)) then
              fp(k,ivpx:ivpz)=(rhopsma*fp(j,ivpx:ivpz) + &
                  rhopbig*fp(k,ivpx:ivpz))/(rhopsma+rhopbig)
            else
              fp(k,ivpx:ivpz)=(rhopbig*fp(j,ivpx:ivpz) + &
                  rhopsma*fp(k,ivpx:ivpz))/(rhopsma+rhopbig)
            endif
!
!  Tag particle for removal by making the radius negative.
!
            fp(j,iap)=-fp(j,iap)
            fp(k,iap)=((rhopsma+rhopbig)*dx*dy*dz/ &
                four_pi_rhopmat_over_three)**(1.0/3.0)
            fp(k,inpswarm)=1/(dx*dy*dz)
          else
            fp(j,iap)=apnew
            fp(k,iap)=apnew
            fp(j,inpswarm)=npnew
            fp(k,inpswarm)=npnew
!            fp(j,ivpx:ivpz)=0.5*(fp(j,ivpx:ivpz)+fp(k,ivpx:ivpz))
!            fp(k,ivpx:ivpz)=fp(j,ivpx:ivpz)
          endif
        else
!
!  Physical particles in the two swarms collide asymmetrically.
!
          if (fp(k,iap)<fp(j,iap)) then
            fp(k,inpswarm)=fp(k,inpswarm)* &
                 (1/(1.0+(fp(j,iap)/fp(k,iap))**3))
            fp(k,iap)=(fp(k,iap)**3+fp(j,iap)**3)**(1.0/3.0)
          else
            fp(k,inpswarm)=0.5*fp(k,inpswarm)
            fp(k,iap)=2**(1.0/3.0)*fp(k,iap)
          endif
        endif
      else
!
!  In case we only evolve the particle radius, number their number.
!  We assume that the mass represented by each superparticle is the same.
!
        if (lcoag_simultaneous) then
          fp(j,iap)=2**(1.0/3.0)*max(fp(j,iap),fp(k,iap))
          fp(k,iap)=fp(j,iap)
        else
          if (fp(k,iap)<fp(j,iap)) then
            fp(k,iap)=(fp(k,iap)**3+fp(j,iap)**3)**(1.0/3.0)
          else
            fp(k,iap)=2**(1.0/3.0)*fp(k,iap)
          endif
        endif
      endif
!
    endsubroutine coagulation_fragmentation
!***********************************************************************
    subroutine particles_coag_timestep_zd(fp)
!
!  Calculate the next timestep with the same scheme as in Zsom & Dullemond
!  (2008). Only use for kernel tests. Self-collision included.
!
!  15-nov-12/KWJ: coded
!
      use General, only: random_number_wrapper
      ! Get all variables
      real, dimension (npar,mpvar) :: fp
!
      real :: r, kernel, dt1_coag_par
      integer :: j, l
!
!  If first timestep calculate the r_ik-matrix and the total rates. 
!  Loop over all particle pairs.
!
!  We also need to get the cumulative distributions to pick the particles.
!  N.B. Non-normalized cumulative distributions for easier handle of coagulation
!  outcome.
!
!
      if (it==1) then
        j=1
        l=1
!  Initialize all rates and cumulative functions
        r_i_tot = 0.
        r_total = 0.
        cum_func_first_i = 0.
        cum_func_sec_ik = 0.
!  Loop over all particles
        do while (.true.)
          do while (.true.)
!  Get kernel for jl-pair. Rep. particle j and phys. particle l
            if (lconstant_kernel_test) then
              kernel=kernel_cst
            elseif (llinear_kernel_test) then
              kernel=kernel_lin* &
                four_pi_rhopmat_over_three*(fp(j,iap)**3+fp(l,iap)**3)
            elseif (lproduct_kernel_test) then
              kernel=kernel_pro* &
                four_pi_rhopmat_over_three2*fp(j,iap)**3*fp(l,iap)**3
            endif
!  Calculate rates and cum. functions
            r_ik_mat(j,l) = fp(l,inpswarm)*kernel    ! rate matrix element
            r_i_tot(j) = r_i_tot(j) + r_ik_mat(j,l)  ! total rate for part. j
            cum_func_sec_ik(j,l) = r_i_tot(j)        ! cum. func. for phys. part.
!
            l = l + 1                      ! go to next particle
            if(l == npar + 1) exit             ! exit if last particle
          enddo
          r_total = r_total + r_i_tot(j)   ! total rate for all particles
          cum_func_first_i(j) = r_total    ! cum. func. for rep. part.
!
          l = 1
          j = j + 1                        ! go to next particle
          if(j == npar + 1) exit
        enddo
      endif
!
!  Now get the timestep (time until next collision) and put in 
!  inverse time-step array. Use eq. 5 i Zsom & Dullemond 2008
!
      call random_number_wrapper(r)
      dt1_coag_par = -r_total/log(r)
      dt1_max = dt1_coag_par
!
    endsubroutine particles_coag_timestep_zd
!***********************************************************************
    subroutine particles_coag_outcome_zd(fp)
!
!  Calculate which two particles collide and the outcome with the same 
!  scheme as in Zsom & Dullemond (2008). Only use for kernel tests.
!  Assume perfect sticking for now.
!
!  15-nov-12/KWJ: coded
!
      use General, only: random_number_wrapper
      ! Get all variables
      real, dimension (npar,mpvar) :: fp
!
      real :: r, r_i_old, kernel
      integer :: i, j, k, l
      real :: delta_r       ! change in r_i rate
!
!  Get the particles that is involved in the collision with 
!  random numbers and cumulative distribution functions. 
!  N.B. Using un-normalized cumulative functions.
!  Use the bisection method to find the right particle.
!
!  Find the representative particle, j.
      call random_number_wrapper(r)
      r = r*r_total
      call particles_coagulation_bisection(cum_func_first_i, r, j)
!  Find the physical particle, l.
      call random_number_wrapper(r)
      r = r*r_i_tot(j)
      call particles_coagulation_bisection(cum_func_sec_ik(j,:), r, l)
!
!  Change properties of the representative particle (particle j). Change
!  mass (radius) and particle density in superparticle j.
!
      fp(j,inpswarm)=fp(j,inpswarm)* &
           (1.0/(1.0+(fp(l,iap)/fp(j,iap))**3))
      fp(j,iap)=(fp(j,iap)**3+fp(l,iap)**3)**(1.0/3.0)
!
!  Update r_ik matrix. Only need to update elements with representative 
!  particle i = j (column) and physical particle k = j (row).
!
!  Loop over i and update row with k = j. Skip k = i = j 
!  (value in i = j column).
!
      i = 1
      do while (.true.)
        if (i /= j) then
          if (lconstant_kernel_test) then
            kernel=kernel_cst
          elseif (llinear_kernel_test) then
            kernel=kernel_lin* &
                   four_pi_rhopmat_over_three*(fp(i,iap)**3+fp(j,iap)**3)
          elseif (lproduct_kernel_test) then
            kernel=kernel_pro* &
                   four_pi_rhopmat_over_three2*fp(i,iap)**3*fp(j,iap)**3
          endif
          delta_r = fp(j,inpswarm)*kernel - &
                    r_ik_mat(i,j)               ! change in matrix element
          r_ik_mat(i,j) = fp(j,inpswarm)*kernel ! rate matrix element
!  Update rate and cumulative function
          r_i_tot(i) = r_i_tot(i) + delta_r  ! total rate for rep. particle i
!  Only need to update cumulative function from element j and onwards
          cum_func_sec_ik(i, j:) = cum_func_sec_ik(i, j:) + &
                                   delta_r
        end if
! 
        i = i + 1            ! go to next superparticle
        if(i == npar + 1) exit
      enddo
!
!  Update i = j column.
!
      k = 1
      do while (.true.)
        if (lconstant_kernel_test) then
          kernel=kernel_cst
        elseif (llinear_kernel_test) then
          kernel=kernel_lin* &
                 four_pi_rhopmat_over_three*(fp(j,iap)**3+fp(k,iap)**3)
        elseif (lproduct_kernel_test) then
          kernel=kernel_pro* &
                 four_pi_rhopmat_over_three2*fp(j,iap)**3*fp(k,iap)**3
        endif
        r_ik_mat(j,k) = fp(k,inpswarm)*kernel   ! rate matrix element
! Update rate for representative particle j
        if (k == 1) then
          r_i_tot(j) = r_ik_mat(j,k)
        else
          r_i_tot(j) = r_i_tot(j) + r_ik_mat(j,k)
        endif
        cum_func_sec_ik(j,k) = r_i_tot(j)
!
        k = k + 1
        if (k == npar + 1) exit
      enddo
!
!  Update total rate and cumulative function for rep. particles
!   
    i = 1
    do while (.true.)
      if (i == 1) then
        r_total = r_i_tot(i)
      else
        r_total = r_total + r_i_tot(i)
      end if
      cum_func_first_i(i) = r_total
!
      i = i + 1
      if (i == npar + 1) exit
    end do
!
    endsubroutine particles_coag_outcome_zd
!***********************************************************************
    subroutine particles_coagulation_bisection(qArr, qVal, iPart)
!
!  Given random value qVal [0,max_rate[ find particle iPart through the bisection 
!  method given cumulative function qArr.
!
!  15-nov-12/KWJ: coded
!
      real, dimension (npar) :: qArr
      real :: qVal
      integer :: iPart, jl, ju, jm
!
      intent (in) :: qArr,qVal
      intent (out) :: iPart
!
      jl = 1                         ! lower index limit
      ju = size(qArr)                ! upper index limit
      if (qVal <= qArr(1)) then      ! qVal is in first bin
        iPart = 1
      else 
        do while ((ju-jl) > 1)       ! not yet done
          jm = (ju+jl)/2             ! midpoint
          if (qVal >= qArr(jm)) then
            jl = jm                  ! new lower limit
          else
            ju = jm                  ! or new upper limit
          endif
        enddo
        iPart = ju                   ! qVal is in the ju-bin
      endif
!
    endsubroutine particles_coagulation_bisection
!***********************************************************************
    subroutine read_particles_coag_run_pars(unit,iostat)
!
!  Read run parameters from run.in.
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,nml=particles_coag_run_pars,err=99,iostat=iostat)
      else
        read(unit,nml=particles_coag_run_pars,err=99)
      endif
!
99    return
!
    endsubroutine read_particles_coag_run_pars
!***********************************************************************
    subroutine write_particles_coag_run_pars(unit)
!
!  Write run parameters to param.nml.
!
      integer, intent(in) :: unit
!
      write(unit,NML=particles_coag_run_pars)
!
    endsubroutine write_particles_coag_run_pars
!*******************************************************************
    subroutine rprint_particles_coagulation(lreset,lwrite)
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
        idiag_ncoagpm=0; idiag_ncoagpartpm=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ncoagpm',idiag_ncoagpm)
        call parse_name(iname,cname(iname),cform(iname), &
            'ncoagpartpm',idiag_ncoagpartpm)
      enddo
!
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_particles_coagulation
!***********************************************************************
endmodule Particles_coagulation
