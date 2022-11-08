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
  real :: tstart_droplet_coagulation=impossible
  real :: reference_radius=5e-5
  logical :: ldroplet_coagulation_runtime=.false.
  logical :: lcoag_simultaneous=.false., lnoselfcollision=.true.
  logical :: lshear_in_vp=.true.
  logical :: lkernel_test=.false., lconstant_kernel_test=.false.
  logical :: llinear_kernel_test=.false., lproduct_kernel_test=.false.
  logical :: lgravitational_cross_section=.false.
  logical :: lzsomdullemond=.false.     ! use Zsom and Dullemond method
  logical :: lconstant_deltav=.false.   ! use constant relative velocity
  logical :: lmaxwell_deltav=.false.    ! use maxwellian relative velocity
  logical :: ldroplet_coagulation=.false., normal_coagulation=.false.
  logical :: lcollision_output=.false., lcollision_output_swapped=.false.,luser_random_number_wrapper=.false.
  logical :: lrelabelling=.false.
  logical :: kernel_output=.false., radius_output=.false.
  logical :: sphericalKernel=.false.
  logical :: lcheck_reference_radius=.false.
  logical :: lremove_particle_phys=.true., lremove_particle=.false., lremove_particle2=.false.
  character (len=labellen) :: droplet_coagulation_model='standard'
!
  real, dimension(:,:), allocatable :: r_ik_mat, cum_func_sec_ik
  real, dimension(:), allocatable :: r_i_tot, cum_func_first_i
  real :: r_total = 0.0
  real :: delta_r = 0.0
  real :: deltav = 1.0          ! relative velocity
  real :: maxwell_param = 1.0   ! alpha parameter for maxwell distribution
  real :: rdifference = 1.0
!
  integer :: idiag_ncoagpm=0, idiag_ncoagpartpm=0, idiag_dt1_coag_par=0
!
  real :: deltad = 1., a0 = 1.
  real :: r1, r2, r3, r4, r5, r6, r7, r8, r_diff
  integer :: idiag_k100_100=0, idiag_k100_80=0, idiag_k100_60=0, idiag_k100_50=0, &
             idiag_k100_40=0, idiag_k100_30=0, idiag_k100_20=0, idiag_k100_10=0, &
             idiag_k80_80=0, idiag_k80_60=0, idiag_k80_50=0, idiag_k80_40=0, &
             idiag_k80_30=0, idiag_k80_20=0, idiag_k80_10=0, idiag_k60_60=0, &
             idiag_k60_50=0, idiag_k60_40=0, idiag_k60_30=0, idiag_k60_20=0, &
             idiag_k60_10=0, idiag_k50_50=0, idiag_k50_40=0, idiag_k50_30=0, &
             idiag_k50_20=0, idiag_k50_10=0, idiag_k40_40=0, idiag_k40_30=0, &
             idiag_k40_20=0, idiag_k40_10=0, idiag_k30_30=0, idiag_k30_20=0, &
             idiag_k30_10=0, idiag_k20_20=0, idiag_k20_10=0, idiag_k10_10=0
  integer :: rbin = 0
!
  namelist /particles_coag_run_pars/ &
      cdtpcoag, lcoag_simultaneous, lshear_in_vp, lconstant_kernel_test, &
      kernel_cst, llinear_kernel_test, kernel_lin, lproduct_kernel_test, &
      kernel_pro, lnoselfcollision, lgravitational_cross_section, &
      GNewton, deltav_grav_floor, critical_mass_ratio_sticking, &
      minimum_particle_mass, minimum_particle_radius, lzsomdullemond, &
      lconstant_deltav, lmaxwell_deltav, deltav, maxwell_param, &
      ldroplet_coagulation, droplet_coagulation_model, lcollision_output, &
      luser_random_number_wrapper, lrelabelling, rdifference, &
      kernel_output, deltad, a0, rbin, &
      radius_output, r1, r2, r3, r4, r5, r6, r7, r8, r_diff, &
      sphericalKernel, normal_coagulation, tstart_droplet_coagulation, &
      lcheck_reference_radius, reference_radius, &
      lremove_particle_phys, lremove_particle, lremove_particle2, &
      lcollision_output_swapped
!
  contains
!***********************************************************************
    subroutine initialize_particles_coag(f)
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
      real, dimension (mpar_loc,mparray) :: fp
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
!  Special treatment for kernel tests.
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
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: xpj, xpk, vpj, vpk
      real :: lambda_mfp1, deltavjk, tau_coll1, prob, rran, kernel
      real :: npswarmj, npswarmk
      integer :: l, j, k, ncoll, ncoll_par, npart_par
!
      intent (in) :: ineargrid
      intent (inout) :: fp
!
      real, dimension (10,10) :: kernel_array
      real, dimension (10) :: radius_all
      real, dimension (10) :: radius_ratio
      real :: rmin, rmax
      integer :: ibin, ik, ij, ikernel, row
      integer, parameter :: max_rows = 10
      real, parameter :: radius_diff=5.e-6
      real, dimension (rbin) :: adplus, adminus, ad
      character(len=50) :: itn,filename
      real :: r100p, r100m, r80p, r80m, r60p, r60m, r50p, r50m, r40p, r40m, r30p, r30m, &
              r20p, r20m, r10p, r10m 
      integer :: k100_100, k100_80, k100_60, k100_50, k100_40, k100_30, k100_20, k100_10, &
                 k80_80, k80_60, k80_50, k80_40, k80_30, k80_20, k80_10, &
                 k60_60, k60_50, k60_40, k60_30, k60_20, k60_10, &
                 k50_50, k50_40, k50_30, k50_20, k50_10, &
                 k40_40, k40_30, k40_20, k40_10, &
                 k30_30, k30_20, k30_10, &
                 k20_20, k20_10, &
                 k10_10
!
! a flag to start collision on the fly
      if (t >= tstart_droplet_coagulation) ldroplet_coagulation=.true. 
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
              k100_100=0; k100_80=0; k100_60=0; k100_50=0; k100_40=0; k100_30=0; k100_20=0; k100_10=0
              k80_80=0; k80_60=0; k80_50=0; k80_40=0; k80_30=0; k80_20=0; k80_10=0
              k60_60=0; k60_50=0; k60_40=0; k60_30=0; k60_20=0; k60_10=0
              k50_50=0; k50_40=0; k50_30=0; k50_20=0; k50_10=0
              k40_40=0; k40_30=0; k40_20=0; k40_10=0
              k30_30=0; k30_20=0; k30_10=0
              k20_20=0; k20_10=0
              k10_10=0
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
                if (ldroplet_coagulation .or. &
                    (sum((vpk-vpj)*(xpk-xpj))<0.0) .or.&
                    lkernel_test) then
!
!  Relative particle speed.
!
                  if (sphericalKernel) then
                    deltavjk=2*sqrt(sum(((vpk-vpj)*(xpk-xpj))**2))/sqrt(sum((xpk-xpj)**2))
                  else
                    deltavjk=sqrt(sum((vpk-vpj)**2))
                  endif
!
! Check for consistensy
!
                  if (droplet_coagulation_model=='standard') then
                    if (lcoag_simultaneous) then
                      call fatal_error('particles_coagulation_pencils',&
                          'One can not set droplet_coagulation_model=standard for lcoag_simultaneous=T')
                    endif
                  else if (droplet_coagulation_model=='shima') then
                    if (.not. lcoag_simultaneous) then
                      call fatal_error('particles_coagulation_pencils',&
                          'One can not set droplet_coagulation_model=shima for lcoag_simultaneous=F')
                    endif
                  endif
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
!      of lighter particles is the collision time-scale times the mass ratio
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
                    if (ldroplet_coagulation .and. .not. lcoag_simultaneous) then
                      tau_coll1=deltavjk*pi*(fp(k,iap)+fp(j,iap))**2*npswarmj
                    elseif (ldroplet_coagulation .and. &
                            lcoag_simultaneous   .and. &
                            droplet_coagulation_model=='shima') then
                      tau_coll1=deltavjk*pi*(fp(k,iap)+fp(j,iap))**2* &
                          max(npswarmj,npswarmk)
                    elseif (normal_coagulation) then 
                      tau_coll1=deltavjk*pi*(fp(k,iap)+fp(j,iap))**2* &
                          min(npswarmj,npswarmk)
                    else
                      tau_coll1=0     
                      if (lgravitational_cross_section) then
                        tau_coll1=tau_coll1*(1.0+2*GNewton* &
                            four_pi_rhopmat_over_three* &
                            (fp(j,iap)**3+fp(k,iap)**3)/((fp(k,iap)+fp(k,iap))* &
                            (deltavjk+deltav_grav_floor)**2))
                      endif
                    endif
                  endif
!
                  if (tau_coll1/=0.0) then
!
!  The probability for a collision in this time-step is dt/tau_coll.
!
                    prob=dt*tau_coll1
!NILS: I found that when the particles_coagulation module was turned on, 
!NILS: the fluid turbulence was affected. By inspecting the energy spectrum
!NILS: it became apparent that it was primarily the small scales that were
!NILS: affected. This happened even though all back reactions from the 
!NILS: particles to the fluid were turned off. After some intense debugging 
!NILS: I finally found that the problem was due to a call to
!NILS: the random_number_wrapper routine. When I replaced this call with a call
!NILS: to the intinsic random_number subroutine, the energy spectrum was fine.
!NILS: 
!NILS: I am not able to understand why the call to the random_number_wrapper
!NILS: routine should cause any problems, but I have reproduce the problem
!NILS: both on our linux cluster (with pfg90) and on my laptop (with gfortran).
                    if (luser_random_number_wrapper) then
                      call random_number_wrapper(rran)
                    else
                      call random_number(rran)
                    endif
                    if (rran<=prob) then
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
!17-06-21: Xiang-Yu coded: kernel of ri rj, diagnostics as time series
                      if (radius_output) then
                        !logorithmic binning
                        do ibin = 1,rbin
                          ad(ibin) = a0*deltad**ibin
                        enddo
                        do ibin = 1,rbin
                          if (abs(ad(ibin)-r1) == minval(abs(ad-r1))) then
                            r100p = ad(ibin+1)
                            r100m = ad(ibin-1)
                          endif
                          if (abs(ad(ibin)-r2) == minval(abs(ad-r2))) then
                            r80p = ad(ibin+1)
                            r80m = ad(ibin-1)
                          endif
                          if (abs(ad(ibin)-r3) == minval(abs(ad-r3))) then
                            r60p = ad(ibin+1)
                            r60m = ad(ibin-1)
                          endif
                          if (abs(ad(ibin)-r4) == minval(abs(ad-r4))) then
                            r50p = ad(ibin+1)
                            r50m = ad(ibin-1)
                          endif
                          if (abs(ad(ibin)-r5) == minval(abs(ad-r5))) then
                            r40p = ad(ibin+1)
                            r40m = ad(ibin-1)
                          endif
                          if (abs(ad(ibin)-r6) == minval(abs(ad-r6))) then
                            r30p = ad(ibin+1)
                            r30m = ad(ibin-1)
                          endif
                          if (abs(ad(ibin)-r7) == minval(abs(ad-r7))) then
                            r20p = ad(ibin+1)
                            r20m = ad(ibin-1)
                          endif
                          if (abs(ad(ibin)-r8) == minval(abs(ad-r8))) then
                            r10p = ad(ibin+1)
                            r10m = ad(ibin-1)
                          endif
                        enddo
                        if (fp(j,iap)>=fp(k,iap)) then
                          !r1=100 mum
                          if (fp(j,iap)>=r100m .and. fp(j,iap)<=r100p .and. fp(k,iap)>=r100m &
                             .and. fp(k,iap)<=r100p) k100_100 = k100_100 + 1
                          if (fp(j,iap)>=r100m .and. fp(j,iap)<=r100p .and. fp(k,iap)>=r80m &
                             .and. fp(k,iap)<=r80p) k100_80 = k100_80 + 1
                          if (fp(j,iap)>=r100m .and. fp(j,iap)<=r100p .and. fp(k,iap)>=r60m &
                             .and. fp(k,iap)<=r60p) k100_60 = k100_60 + 1
                          if (fp(j,iap)>=r100m .and. fp(j,iap)<=r100p .and. fp(k,iap)>=r50m &
                             .and. fp(k,iap)<=r50p) k100_50 = k100_50 + 1
                          if (fp(j,iap)>=r100m .and. fp(j,iap)<=r100p .and. fp(k,iap)>=r40m &
                             .and. fp(k,iap)<=r40p) k100_40 = k100_40 + 1
                          if (fp(j,iap)>=r100m .and. fp(j,iap)<=r100p .and. fp(k,iap)>=r30m &
                             .and. fp(k,iap)<=r30p) k100_30 = k100_30 + 1
                          if (fp(j,iap)>=r100m .and. fp(j,iap)<=r100p .and. fp(k,iap)>=r20m &
                             .and. fp(k,iap)<=r20p) k100_20 = k100_20 + 1
                          if (fp(j,iap)>=r100m .and. fp(j,iap)<=r100p .and. fp(k,iap)>=r10m &
                             .and. fp(k,iap)<=r10p) k100_10 = k100_10 + 1
                          !r2=80 mum
                          if (fp(j,iap)>=r80m .and. fp(j,iap)<=r80p .and. fp(k,iap)>=r80m .and. fp(k,iap)<=r80p) k80_80 = k80_80 + 1
                          if (fp(j,iap)>=r80m .and. fp(j,iap)<=r80p .and. fp(k,iap)>=r60m .and. fp(k,iap)<=r60p) k80_60 = k80_60 + 1
                          if (fp(j,iap)>=r80m .and. fp(j,iap)<=r80p .and. fp(k,iap)>=r50m .and. fp(k,iap)<=r50p) k80_50 = k80_50 + 1
                          if (fp(j,iap)>=r80m .and. fp(j,iap)<=r80p .and. fp(k,iap)>=r40m .and. fp(k,iap)<=r40p) k80_40 = k80_40 + 1
                          if (fp(j,iap)>=r80m .and. fp(j,iap)<=r80p .and. fp(k,iap)>=r30m .and. fp(k,iap)<=r30p) k80_30 = k80_30 + 1
                          if (fp(j,iap)>=r80m .and. fp(j,iap)<=r80p .and. fp(k,iap)>=r20m .and. fp(k,iap)<=r20p) k80_20 = k80_20 + 1
                          if (fp(j,iap)>=r80m .and. fp(j,iap)<=r80p .and. fp(k,iap)>=r10m .and. fp(k,iap)<=r10p) k80_10 = k80_10 + 1
                          !r3=60 mum
                          if (fp(j,iap)>=r60m .and. fp(j,iap)<=r60p .and. fp(k,iap)>=r60m .and. fp(k,iap)<=r60p) k60_60 = k60_60 + 1
                          if (fp(j,iap)>=r60m .and. fp(j,iap)<=r60p .and. fp(k,iap)>=r50m .and. fp(k,iap)<=r50p) k60_50 = k60_50 + 1
                          if (fp(j,iap)>=r60m .and. fp(j,iap)<=r60p .and. fp(k,iap)>=r40m .and. fp(k,iap)<=r40p) k60_40 = k60_40 + 1
                          if (fp(j,iap)>=r60m .and. fp(j,iap)<=r60p .and. fp(k,iap)>=r30m .and. fp(k,iap)<=r30p) k60_30 = k60_30 + 1
                          if (fp(j,iap)>=r60m .and. fp(j,iap)<=r60p .and. fp(k,iap)>=r20m .and. fp(k,iap)<=r20p) k60_20 = k60_20 + 1
                          if (fp(j,iap)>=r60m .and. fp(j,iap)<=r60p .and. fp(k,iap)>=r10m .and. fp(k,iap)<=r10p) k60_10 = k60_10 + 1
                          !r4=50 mum                                                                           
                          if (fp(j,iap)>=r50m .and. fp(j,iap)<=r50p .and. fp(k,iap)>=r50m .and. fp(k,iap)<=r50p) k50_50 = k50_50 + 1
                          if (fp(j,iap)>=r50m .and. fp(j,iap)<=r50p .and. fp(k,iap)>=r40m .and. fp(k,iap)<=r40p) k50_40 = k50_40 + 1
                          if (fp(j,iap)>=r50m .and. fp(j,iap)<=r50p .and. fp(k,iap)>=r30m .and. fp(k,iap)<=r30p) k50_30 = k50_30 + 1
                          if (fp(j,iap)>=r50m .and. fp(j,iap)<=r50p .and. fp(k,iap)>=r20m .and. fp(k,iap)<=r20p) k50_20 = k50_20 + 1
                          if (fp(j,iap)>=r50m .and. fp(j,iap)<=r50p .and. fp(k,iap)>=r10m .and. fp(k,iap)<=r10p) k50_10 = k50_10 + 1
                          !r5=40 mum                                                                           
                          if (fp(j,iap)>=r40m .and. fp(j,iap)<=r40p .and. fp(k,iap)>=r40m .and. fp(k,iap)<=r40p) k40_40 = k40_40 + 1
                          if (fp(j,iap)>=r40m .and. fp(j,iap)<=r40p .and. fp(k,iap)>=r30m .and. fp(k,iap)<=r30p) k40_30 = k40_30 + 1
                          if (fp(j,iap)>=r40m .and. fp(j,iap)<=r40p .and. fp(k,iap)>=r20m .and. fp(k,iap)<=r20p) k40_20 = k40_20 + 1
                          if (fp(j,iap)>=r40m .and. fp(j,iap)<=r40p .and. fp(k,iap)>=r10m .and. fp(k,iap)<=r10p) k40_10 = k40_10 + 1
                          !r6=30 mum                                                                           
                          if (fp(j,iap)>=r30m .and. fp(j,iap)<=r30p .and. fp(k,iap)>=r30m .and. fp(k,iap)<=r30p) k30_30 = k30_30 + 1
                          if (fp(j,iap)>=r30m .and. fp(j,iap)<=r30p .and. fp(k,iap)>=r20m .and. fp(k,iap)<=r20p) k30_20 = k30_20 + 1
                          if (fp(j,iap)>=r30m .and. fp(j,iap)<=r30p .and. fp(k,iap)>=r10m .and. fp(k,iap)<=r10p) k30_10 = k30_10 + 1
                          !r7=20 mum
                          if (fp(j,iap)>=r20m .and. fp(j,iap)<=r20p .and. fp(k,iap)>=r20m .and. fp(k,iap)<=r20p) k20_20 = k20_20 + 1
                          if (fp(j,iap)>=r20m .and. fp(j,iap)<=r20p .and. fp(k,iap)>=r10m .and. fp(k,iap)<=r10p) k20_10 = k20_10 + 1
                          !r7=10 mum
                          if (fp(j,iap)>=r10m .and. fp(j,iap)<=r10p .and. fp(k,iap)>=r10m .and. fp(k,iap)<=r10p) k10_10 = k10_10 + 1
                        else
                          !r1=100 mum
                          if (fp(k,iap)>=r100m .and. fp(k,iap)<=r100p .and. fp(j,iap)>=r100m &
                             .and. fp(j,iap)<=r100p) k100_100 = k100_100 + 1
                          if (fp(k,iap)>=r100m .and. fp(k,iap)<=r100p .and. fp(j,iap)>=r80m &
                             .and. fp(j,iap)<=r80p) k100_80 = k100_80 + 1
                          if (fp(k,iap)>=r100m .and. fp(k,iap)<=r100p .and. fp(j,iap)>=r60m &
                             .and. fp(j,iap)<=r60p) k100_60 = k100_60 + 1
                          if (fp(k,iap)>=r100m .and. fp(k,iap)<=r100p .and. fp(j,iap)>=r50m &
                             .and. fp(j,iap)<=r50p) k100_50 = k100_50 + 1
                          if (fp(k,iap)>=r100m .and. fp(k,iap)<=r100p .and. fp(j,iap)>=r40m &
                             .and. fp(j,iap)<=r40p) k100_40 = k100_40 + 1
                          if (fp(k,iap)>=r100m .and. fp(k,iap)<=r100p .and. fp(j,iap)>=r30m &
                             .and. fp(j,iap)<=r30p) k100_30 = k100_30 + 1
                          if (fp(k,iap)>=r100m .and. fp(k,iap)<=r100p .and. fp(j,iap)>=r20m &
                             .and. fp(j,iap)<=r20p) k100_20 = k100_20 + 1
                          if (fp(k,iap)>=r100m .and. fp(k,iap)<=r100p .and. fp(j,iap)>=r10m &
                             .and. fp(j,iap)<=r10p) k100_10 = k100_10 + 1
                          !r2=80 mum
                          if (fp(k,iap)>=r80m .and. fp(k,iap)<=r80p .and. fp(j,iap)>=r80m .and. fp(j,iap)<=r80p) k80_80 = k80_80 + 1
                          if (fp(k,iap)>=r80m .and. fp(k,iap)<=r80p .and. fp(j,iap)>=r60m .and. fp(j,iap)<=r60p) k80_60 = k80_60 + 1
                          if (fp(k,iap)>=r80m .and. fp(k,iap)<=r80p .and. fp(j,iap)>=r50m .and. fp(j,iap)<=r50p) k80_50 = k80_50 + 1
                          if (fp(k,iap)>=r80m .and. fp(k,iap)<=r80p .and. fp(j,iap)>=r40m .and. fp(j,iap)<=r40p) k80_40 = k80_40 + 1
                          if (fp(k,iap)>=r80m .and. fp(k,iap)<=r80p .and. fp(j,iap)>=r30m .and. fp(j,iap)<=r30p) k80_30 = k80_30 + 1
                          if (fp(k,iap)>=r80m .and. fp(k,iap)<=r80p .and. fp(j,iap)>=r20m .and. fp(j,iap)<=r20p) k80_20 = k80_20 + 1
                          if (fp(k,iap)>=r80m .and. fp(k,iap)<=r80p .and. fp(j,iap)>=r10m .and. fp(j,iap)<=r10p) k80_10 = k80_10 + 1
                          !r3=60 mum
                          if (fp(k,iap)>=r60m .and. fp(k,iap)<=r60p .and. fp(j,iap)>=r60m .and. fp(j,iap)<=r60p) k60_60 = k60_60 + 1
                          if (fp(k,iap)>=r60m .and. fp(k,iap)<=r60p .and. fp(j,iap)>=r50m .and. fp(j,iap)<=r50p) k60_50 = k60_50 + 1
                          if (fp(k,iap)>=r60m .and. fp(k,iap)<=r60p .and. fp(j,iap)>=r40m .and. fp(j,iap)<=r40p) k60_40 = k60_40 + 1
                          if (fp(k,iap)>=r60m .and. fp(k,iap)<=r60p .and. fp(j,iap)>=r30m .and. fp(j,iap)<=r30p) k60_30 = k60_30 + 1
                          if (fp(k,iap)>=r60m .and. fp(k,iap)<=r60p .and. fp(j,iap)>=r20m .and. fp(j,iap)<=r20p) k60_20 = k60_20 + 1
                          if (fp(k,iap)>=r60m .and. fp(k,iap)<=r60p .and. fp(j,iap)>=r10m .and. fp(j,iap)<=r10p) k60_10 = k60_10 + 1
                          !r4=50 mum                                                                           
                          if (fp(k,iap)>=r50m .and. fp(k,iap)<=r50p .and. fp(j,iap)>=r50m .and. fp(j,iap)<=r50p) k50_50 = k50_50 + 1
                          if (fp(k,iap)>=r50m .and. fp(k,iap)<=r50p .and. fp(j,iap)>=r40m .and. fp(j,iap)<=r40p) k50_40 = k50_40 + 1
                          if (fp(k,iap)>=r50m .and. fp(k,iap)<=r50p .and. fp(j,iap)>=r30m .and. fp(j,iap)<=r30p) k50_30 = k50_30 + 1
                          if (fp(k,iap)>=r50m .and. fp(k,iap)<=r50p .and. fp(j,iap)>=r20m .and. fp(j,iap)<=r20p) k50_20 = k50_20 + 1
                          if (fp(k,iap)>=r50m .and. fp(k,iap)<=r50p .and. fp(j,iap)>=r10m .and. fp(j,iap)<=r10p) k50_10 = k50_10 + 1
                          !r5=40 mum                                                                           
                          if (fp(k,iap)>=r40m .and. fp(k,iap)<=r40p .and. fp(j,iap)>=r40m .and. fp(j,iap)<=r40p) k40_40 = k40_40 + 1
                          if (fp(k,iap)>=r40m .and. fp(k,iap)<=r40p .and. fp(j,iap)>=r30m .and. fp(j,iap)<=r30p) k40_30 = k40_30 + 1
                          if (fp(k,iap)>=r40m .and. fp(k,iap)<=r40p .and. fp(j,iap)>=r20m .and. fp(j,iap)<=r20p) k40_20 = k40_20 + 1
                          if (fp(k,iap)>=r40m .and. fp(k,iap)<=r40p .and. fp(j,iap)>=r10m .and. fp(j,iap)<=r10p) k40_10 = k40_10 + 1
                          !r6=30 mum                                                                           
                          if (fp(k,iap)>=r30m .and. fp(k,iap)<=r30p .and. fp(j,iap)>=r30m .and. fp(j,iap)<=r30p) k30_30 = k30_30 + 1
                          if (fp(k,iap)>=r30m .and. fp(k,iap)<=r30p .and. fp(j,iap)>=r20m .and. fp(j,iap)<=r20p) k30_20 = k30_20 + 1
                          if (fp(k,iap)>=r30m .and. fp(k,iap)<=r30p .and. fp(j,iap)>=r10m .and. fp(j,iap)<=r10p) k30_10 = k30_10 + 1
                          !r7=20 mum
                          if (fp(k,iap)>=r20m .and. fp(k,iap)<=r20p .and. fp(j,iap)>=r20m .and. fp(j,iap)<=r20p) k20_20 = k20_20 + 1
                          if (fp(k,iap)>=r20m .and. fp(k,iap)<=r20p .and. fp(j,iap)>=r10m .and. fp(j,iap)<=r10p) k20_10 = k20_10 + 1
                          !r7=10 mum
                          if (fp(k,iap)>=r10m .and. fp(k,iap)<=r10p .and. fp(j,iap)>=r10m .and. fp(j,iap)<=r10p) k10_10 = k10_10 + 1
                        endif 
                      endif
!17-06-21: Xiang-Yu

!17-06-18: Xiang-Yu coded: kernel of ri rj, diagnostics as 2-D matrix
                      if (kernel_output) then
! read radius and radius ratio
                        open(unit=11,file="radius.txt")
                        do row = 1, max_rows
                        read(11,*) radius_all(row)
                        enddo
                        close(unit=11)
                        rmin = minval(radius_all)
                        rmax = maxval(radius_all)
                        !open(unit=12,file="ratio.txt")
                        !do row = 1, max_rows
                        !        read(12,*) radius_ratio(row)
                        !enddo
                        !close(unit=12)
                        !logorithmic binning
                        do ibin = 1,rbin
                          adminus(ibin) = a0*deltad**(ibin-1)
                          adplus(ibin) = a0*deltad**ibin
                          ad(ibin) = (adminus(ibin)+adplus(ibin))*.5
                        enddo
                        ! search for collector and collected particles and bin them in the kernel
                        if (fp(k,iap) >= rmin .and. fp(k,iap) <= rmax) then
                        !       ikernel=0
                          do ibin=1,rbin 
                            if (abs(fp(k,iap)-radius_all(ibin)) == minval(abs(fp(k,iap)-radius_all))) then
                            ik=ibin
                            endif
                          enddo
                          do ibin=1,rbin
                            if (abs(fp(j,iap)-radius_all(ibin)) == minval(abs(fp(j,iap)-radius_all))) then
                            ij=ibin
                            endif
                          enddo
                          ! ikernel=ikernel+1
                          !kernel_array(ik,ij)=(kernel_array(ik,ij)+tau_coll1)/ikernel
                          kernel_array(ik,ij) = tau_coll1
                          ! output the kernel             
                          write(itn,'(I5)') it  ! convert  integer to char
                          filename = adjustl(trim(adjustr(itn)))
                          open(13, file=trim(directory_dist)//'/'//trim(filename)//'.dat')
                          write(13,"(10e11.3)") kernel_array
                          close(13)
                          !       print*,'kernel_array=',kernel_array
                          !       print*,'tau_coll1=',tau_coll1
                        endif
                      endif
!17-06-18:XY
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
                    call sum_par_name((/real(ncoll_par)/),idiag_ncoagpm)
                if (idiag_ncoagpartpm/=0) &
                    call sum_par_name((/real(npart_par)/),idiag_ncoagpartpm)
                if (idiag_dt1_coag_par/=0) &
                    call sum_par_name((/real(ncoll_par)/),idiag_dt1_coag_par)
                if (idiag_k100_100/=0) &        
                          call sum_par_name((/real(k100_100)/),idiag_k100_100) 
                if (idiag_k100_80/=0) & 
                          call sum_par_name((/real(k100_80)/),idiag_k100_80) 
                if (idiag_k100_60/=0) & 
                          call sum_par_name((/real(k100_60)/),idiag_k100_60) 
                if (idiag_k100_50/=0) & 
                          call sum_par_name((/real(k100_50)/),idiag_k100_50) 
                if (idiag_k100_40/=0) & 
                          call sum_par_name((/real(k100_40)/),idiag_k100_40) 
                if (idiag_k100_30/=0) & 
                          call sum_par_name((/real(k100_30)/),idiag_k100_30) 
                if (idiag_k100_20/=0) & 
                          call sum_par_name((/real(k100_20)/),idiag_k100_20) 
                if (idiag_k100_10/=0) & 
                          call sum_par_name((/real(k100_10)/),idiag_k100_10) 
                if (idiag_k80_80/=0) call sum_par_name((/real(k80_80)/),idiag_k80_80) 
                if (idiag_k80_60/=0) call sum_par_name((/real(k80_60)/),idiag_k80_60) 
                if (idiag_k80_50/=0) call sum_par_name((/real(k80_50)/),idiag_k80_50) 
                if (idiag_k80_40/=0) call sum_par_name((/real(k80_40)/),idiag_k80_40) 
                if (idiag_k80_30/=0) call sum_par_name((/real(k80_30)/),idiag_k80_30) 
                if (idiag_k80_20/=0) call sum_par_name((/real(k80_20)/),idiag_k80_20) 
                if (idiag_k80_10/=0) call sum_par_name((/real(k80_10)/),idiag_k80_10) 
                if (idiag_k60_60/=0) call sum_par_name((/real(k60_60)/),idiag_k60_60) 
                if (idiag_k60_50/=0) call sum_par_name((/real(k60_50)/),idiag_k60_50) 
                if (idiag_k60_40/=0) call sum_par_name((/real(k60_40)/),idiag_k60_40) 
                if (idiag_k60_30/=0) call sum_par_name((/real(k60_30)/),idiag_k60_30) 
                if (idiag_k60_20/=0) call sum_par_name((/real(k60_20)/),idiag_k60_20) 
                if (idiag_k60_10/=0) call sum_par_name((/real(k60_10)/),idiag_k60_10) 
                if (idiag_k50_50/=0) call sum_par_name((/real(k50_50)/),idiag_k50_50) 
                if (idiag_k50_40/=0) call sum_par_name((/real(k50_40)/),idiag_k50_40) 
                if (idiag_k50_30/=0) call sum_par_name((/real(k50_30)/),idiag_k50_30) 
                if (idiag_k50_20/=0) call sum_par_name((/real(k50_20)/),idiag_k50_20) 
                if (idiag_k50_10/=0) call sum_par_name((/real(k50_10)/),idiag_k50_10) 
                if (idiag_k40_40/=0) call sum_par_name((/real(k40_40)/),idiag_k40_40) 
                if (idiag_k40_30/=0) call sum_par_name((/real(k40_30)/),idiag_k40_30) 
                if (idiag_k40_20/=0) call sum_par_name((/real(k40_20)/),idiag_k40_20) 
                if (idiag_k40_10/=0) call sum_par_name((/real(k40_10)/),idiag_k40_10) 
                if (idiag_k30_30/=0) call sum_par_name((/real(k30_30)/),idiag_k30_30) 
                if (idiag_k30_20/=0) call sum_par_name((/real(k30_20)/),idiag_k30_20) 
                if (idiag_k30_10/=0) call sum_par_name((/real(k30_10)/),idiag_k30_10) 
                if (idiag_k20_20/=0) call sum_par_name((/real(k20_20)/),idiag_k20_20) 
                if (idiag_k20_10/=0) call sum_par_name((/real(k20_10)/),idiag_k20_10) 
                if (idiag_k10_10/=0) call sum_par_name((/real(k10_10)/),idiag_k10_10) 
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
          if (idiag_dt1_coag_par/=0) &
              call sum_par_name(fp(1:npar_loc,ixp),idiag_dt1_coag_par)
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
      real, dimension (mpar_loc,mparray) :: fp
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
!  Change the particle size to the new size. 
!  If droplet_coagulation_model does not equal 'shima' the total mass in the
!  particle swarm is kept constant.
!
!  A representative particle k colliding with a swarm of larger particles j
!  simply obtains the new mass m_k -> m_k + m_j.
!
!  A representative particle k colliding with a swarm of lighter particles j
!  obtains the new mass m_k -> 2*m_k (this is true only when 
!  ldroplet_coagulation=F.)
!
!  24-nov-10/anders: coded
!  25-jan-16/Xiang-yu and nils: added the model of Shima et al. (2009)
!
      use General
!
      real, dimension (mpar_loc,mparray) :: fp
      integer :: j, k, swarm_index1, swarm_index2, iswap, ipar_j_
      real :: deltavjk
      real, dimension (3) :: xpj, xpk, vpj, vpk, vpnew
!
      real, dimension (3) :: vvcm, vvkcm, vvkcmnew, vvkcm_normal, vvkcm_parall
      real, dimension (3) :: nvec
      real :: mpsma, mpbig, npsma, npbig, npnew, mpnew, apnew
      real :: rhopsma, rhopbig, apsma, apbig
      real :: coeff_restitution, deltav_recoil, escape_speed
      real :: apj,mpj,npj,apk,mpk,npk
      logical :: lswap, exists
      character (len=fnlen) :: directory=''
!
      if (lparticles_number) then
!
!  Physical particles in the two swarms collide simultaneously.
!
        if (lcoag_simultaneous) then
!
!  Define mpsma, npsma, rhopsma for the superparticle with lighter particles.
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
!  If we are working on droploet coagulation all colliding droplets are
!  assumed to coalesce. I.e. the mass of the new representative droplet equals
!  the sum of the two initial droplet masses.  
!
          if (ldroplet_coagulation) then
            if (droplet_coagulation_model=='standard') then
              mpnew=mpsma+mpbig
              apnew=(mpnew*three_over_four_pi_rhopmat)**(1.0/3.0)
              npnew=(rhopsma+rhopbig)/(2.*mpnew)
              vpnew=(&
                  fp(j,ivpx:ivpz)*four_pi_rhopmat_over_three*fp(j,iap)**3+&
                  fp(k,ivpx:ivpz)*four_pi_rhopmat_over_three*fp(k,iap)**3)/mpnew
              fp(j,iap)=apnew
              fp(k,iap)=apnew
              fp(j,inpswarm)=npnew
              fp(k,inpswarm)=npnew
              fp(j,ivpx:ivpz)=vpnew
              fp(k,ivpx:ivpz)=vpnew
            else if (droplet_coagulation_model=='shima') then
              mpj = four_pi_rhopmat_over_three*fp(j,iap)**3
              npj = fp(j,inpswarm)
              mpk = four_pi_rhopmat_over_three*fp(k,iap)**3
              npk = fp(k,inpswarm)
!              
!
!  Identify the swarm with the highest particle number density
!
              if (npk<npj) then
                swarm_index1=k
                swarm_index2=j
                lswap=mpk<mpj
              else
                swarm_index1=j
                swarm_index2=k
                lswap=mpk>mpj
              endif
              if (lcollision_output) then
                open(99,POSITION='append', &
               !   FILE=trim(directory_dist)//'/collisions.dat')
                  FILE=trim(directory_snap)//'/collisions.dat') !21-08-14/XYLI
                write(99,"(f18.6,2i10,2f12.8,1p,2e11.3)") t,ipar(j),ipar(k),fp(j,iap),fp(k,iap),fp(j,inpswarm),fp(k,inpswarm)
                close(99)
              endif
!
!  Check if we have the special case where the particle number 
!  densities are equal. This will typically be the case in the beginning of 
!  a simulation.
!
              if (npk == npj) then
!
! 13-June-18/Xiang-Yu: coded
!                if (npk*dx*dy*dz<1.0 .and. lremove_particle) then
                if (lremove_particle.and.lremove_particle_phys) then
                  if (fp(swarm_index1,iap)>=fp(swarm_index2,iap)) then
                    fp(swarm_index1,ivpx:ivpz)=(rhopbig*fp(swarm_index1,ivpx:ivpz) + &
                        rhopsma*fp(swarm_index2,ivpx:ivpz))/(rhopsma+rhopbig)
!                    fp(swarm_index1,inpswarm)=1.0/(dx*dy*dz)
                    fp(swarm_index1,iap)=((rhopsma+rhopbig)/fp(swarm_index1,inpswarm)/ &
                        four_pi_rhopmat_over_three)**(1.0/3.0)
                    fp(swarm_index2,iap)=-fp(swarm_index2,iap)
                  else
                    fp(swarm_index2,ivpx:ivpz)=(rhopbig*fp(swarm_index2,ivpx:ivpz) + &
                        rhopsma*fp(swarm_index1,ivpx:ivpz))/(rhopsma+rhopbig)
!                    fp(swarm_index2,inpswarm)=1.0/(dx*dy*dz)
                    fp(swarm_index2,iap)=((rhopsma+rhopbig)/fp(swarm_index2,inpswarm)/ &
                        four_pi_rhopmat_over_three)**(1.0/3.0)
                    fp(swarm_index1,iap)=-fp(swarm_index1,iap)
                  endif
                else
!
! 13-jun-18/Xiang-Yu.
! 10-aug-21/Axel+Nils: added lremove_particle_phys option
!
                  if (lremove_particle_phys) then
                    fp(swarm_index1,inpswarm)=fp(swarm_index1,inpswarm)/2.
                    fp(swarm_index2,inpswarm)=fp(swarm_index2,inpswarm)/2.
                  endif
                  mpnew=mpj+mpk
                  fp(swarm_index1,iap)&
                      =(mpnew*three_over_four_pi_rhopmat)**(1.0/3.0)
                  fp(swarm_index2,iap)&
                      =(mpnew*three_over_four_pi_rhopmat)**(1.0/3.0)
                  vpnew=(fp(j,ivpx:ivpz)*mpj+fp(k,ivpx:ivpz)*mpk)/mpnew
                  fp(swarm_index1,ivpx:ivpz)=vpnew
                  fp(swarm_index2,ivpx:ivpz)=vpnew
                endif
              else
!
!  Set particle radius. The radius of the swarm with the highes particle 
!  number density does not change since its mass doesn't change.
!
                mpnew=mpj+mpk
                fp(swarm_index1,iap)&
                    =(mpnew*three_over_four_pi_rhopmat)**(1.0/3.0)
!
!  Set particle number densities. Only the swarm with the highest 
!  particle number density, changes its particle number density.
!  Added possibility to not do this (which is then no longer mass-conserving)
!
                if (lremove_particle_phys) then
                  fp(swarm_index2,inpswarm)=abs(npj-npk)
                endif
!
!  Set particle velocities. Only the swarm with the lowest 
!  particle number density, changes its particle velcity.
!
                vpnew=(fp(j,ivpx:ivpz)*mpj+fp(k,ivpx:ivpz)*mpk)/mpnew
                fp(swarm_index1,ivpx:ivpz)=vpnew
              endif
!
!11-Feb-20/Xiang-Yu:for the purpose of particle removing scheme, for all cases, not only for "nx=ny".              
              npnew=(rhopsma+rhopbig)/(2.*mpnew)
!11-Feb-20/Xiang-Yu.
!
! 10-July-18/Xiang-Yu coded: remove superparticle containing less than one droplet.
              if (npnew*dx*dy*dz<1.0 .and. lremove_particle2) then
                if (fp(j,iap)<fp(k,iap)) then
                  fp(k,ivpx:ivpz)=(rhopsma*fp(j,ivpx:ivpz) + &
                      rhopbig*fp(k,ivpx:ivpz))/(rhopsma+rhopbig)
                else
                  fp(k,ivpx:ivpz)=(rhopbig*fp(j,ivpx:ivpz) + &
                      rhopsma*fp(k,ivpx:ivpz))/(rhopsma+rhopbig)
                endif
                fp(j,iap)=-fp(j,iap)
                fp(k,iap)=((rhopsma+rhopbig)*dx*dy*dz/ &
                    four_pi_rhopmat_over_three)**(1.0/3.0)
                fp(k,inpswarm)=1/(dx*dy*dz)
              endif
! 10-July-18/Xiang-Yu.
!
!  swap names ipar to follow lucky droplets
!
              if (lrelabelling) then
                !if (lswap.and.(fp(j,iap)/=fp(k,iap))) then
                if (lswap.and.(abs(fp(j,iap)-fp(k,iap)))>rdifference) then
                  ipar_j_=ipar(k)
                  ipar(k)=ipar(j)
                  ipar(j)=ipar_j_
                  iswap=1
                else
                  iswap=0
                endif
              endif

!
              if (lcollision_output_swapped) then
                open(99,POSITION='append', &
                  FILE=trim(directory_dist)//'/collisions_swapped.dat')
                  write(99,"(f18.6,2i10,2f12.8,1p,2e11.3)") &
                    t,ipar(j),ipar(k),fp(j,iap),fp(k,iap), &
                    fp(j,inpswarm),fp(k,inpswarm),iswap
                close(99)
              endif
!
            else
              call fatal_error('','No such droplet_coagulation_model')
            endif
          else
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
            endif
          endif ! if (ldroplet_coagulation)
        else
!
!  If we are working on droploet coagulation all colliding droplets are
!  assumed to coalesce. I.e. the mass of the new representative droplet equals
!  the sum of the two initial droplet masses.  
!
          if (ldroplet_coagulation) then
            apj = fp(j,iap)
            mpj = four_pi_rhopmat_over_three*fp(j,iap)**3
            npj = fp(j,inpswarm)
            apk = fp(k,iap)
            mpk = four_pi_rhopmat_over_three*fp(k,iap)**3
            npk = fp(k,inpswarm)
            if (droplet_coagulation_model=='standard') then
              if (lcollision_output) then
                open(99,POSITION='append', &
                  FILE=trim(directory_dist)//'/collisions.dat')
                write(99,"(f18.6,2i10,2f12.8,1p,2e11.3)") t,ipar(j),ipar(k),fp(j,iap),fp(k,iap),fp(j,inpswarm),fp(k,inpswarm)
                close(99)
              endif
              mpnew=mpj+mpk
              apnew=(mpnew*three_over_four_pi_rhopmat)**(1.0/3.0)
!            npnew=(mpj*npk+mpk*npk)/(2.*mpnew)
              npnew=mpk*npk/mpnew
              vpnew=(fp(j,ivpx:ivpz)*mpj+fp(k,ivpx:ivpz)*mpk)/mpnew
              fp(k,iap)=apnew
              fp(k,inpswarm)=npnew
              fp(k,ivpx:ivpz)=vpnew
            else if (droplet_coagulation_model=='shima') then
!
!  Identify the swarm with the highest particle number density
!
              if (npk<npj) then
                swarm_index1=k
                swarm_index2=j
              else
                swarm_index1=j
                swarm_index2=k
              endif
!
!  Check if we have the special case where the particle number 
!  densities are equal. This will typically be the case in the beginning of 
!  a simulation.
!
              if (npk == npj) then
                fp(swarm_index1,inpswarm)=fp(swarm_index1,inpswarm)/2.
                fp(swarm_index2,inpswarm)=fp(swarm_index2,inpswarm)/2.
                mpnew=mpj+mpk
                fp(swarm_index1,iap)&
                    =(mpnew*three_over_four_pi_rhopmat)**(1.0/3.0)
                fp(swarm_index2,iap)&
                    =(mpnew*three_over_four_pi_rhopmat)**(1.0/3.0)
                vpnew=(fp(j,ivpx:ivpz)*mpj+fp(k,ivpx:ivpz)*mpk)/mpnew
                fp(swarm_index1,ivpx:ivpz)=vpnew
                fp(swarm_index2,ivpx:ivpz)=vpnew
              else
!
!  Set particle radius. The radius of the swarm with the highes particle 
!  number density does not change since its mass doesn't change.
!
                mpnew=mpj+mpk
                fp(swarm_index1,iap)&
                    =(mpnew*three_over_four_pi_rhopmat)**(1.0/3.0)
!
!  Set particle number densities. Only the the swarm with the highest 
!  particle number density, changes its particle number density.
!
                fp(swarm_index2,inpswarm)=abs(npj-npk)
!
!  Set particle velocities. Only the swarm with the lowest 
!  particle number density, changes its particle velcity.
!
                vpnew=(fp(j,ivpx:ivpz)*mpj+fp(k,ivpx:ivpz)*mpk)/mpnew
                fp(swarm_index1,ivpx:ivpz)=vpnew
              endif
            else
              call fatal_error('','No such droplet_coagulation_model')
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
!  Check whether the largest article has reached a certain reference radius
!
      if (lcheck_reference_radius) then
        if (fp(j,iap)>reference_radius.or.fp(k,iap)>reference_radius) then
          call safe_character_assign(directory,trim(datadir)//'/proc'//itoa(iproc))
          inquire(FILE=trim(directory)//'/reached_reference_radius.dat',EXIST=exists)
          if (.not.exists) then
            open(unit=11,file=trim(directory)//'/reached_reference_radius.dat',STATUS='unknown')
            write(11,*) t,fp(j,iap),fp(k,iap)
            close(11)
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
!  08-jan-13/KWJ: modified
!
      use General, only: random_number_wrapper
      ! Get all variables
!NILS: Shouldn't npar be mpar_loc?
!      real, dimension (npar,mparray) :: fp
      real, dimension (mpar_loc,mparray) :: fp
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
!  Theoretical kernels with analytic solutions
            if (lconstant_kernel_test) then
              kernel=kernel_cst
            elseif (llinear_kernel_test) then
              kernel=kernel_lin* &
                four_pi_rhopmat_over_three*(fp(j,iap)**3+fp(l,iap)**3)
            elseif (lproduct_kernel_test) then
              kernel=kernel_pro* &
                four_pi_rhopmat_over_three2*fp(j,iap)**3*fp(l,iap)**3
!  Physical kernels
            elseif (lconstant_deltav) then    ! constant relative velocity
              kernel=deltav*pi*(fp(j,iap)+fp(l,iap))**2
            elseif (lmaxwell_deltav) then  ! maxwellian relative velocity
              call particles_coag_maxwell(deltav, maxwell_param)
              kernel=deltav*pi*(fp(j,iap)+fp(l,iap))**2
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
      dt1_max = max(dt1_max,dt1_coag_par)
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
!NILS: Shouldn't npar be mpar_loc?
!      real, dimension (npar,mparray) :: fp
      real, dimension (mpar_loc,mparray) :: fp
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
!  Theoretical kernels with analytic solutions
          if (lconstant_kernel_test) then
            kernel=kernel_cst
          elseif (llinear_kernel_test) then
            kernel=kernel_lin* &
                   four_pi_rhopmat_over_three*(fp(i,iap)**3+fp(j,iap)**3)
          elseif (lproduct_kernel_test) then
            kernel=kernel_pro* &
                   four_pi_rhopmat_over_three2*fp(i,iap)**3*fp(j,iap)**3
!  Physical kernels
          elseif (lconstant_deltav) then    ! constant relative velocity
            kernel=deltav*pi*(fp(j,iap)+fp(l,iap))**2
          elseif (lmaxwell_deltav) then  ! maxwellian relative velocity
              call particles_coag_maxwell(deltav, maxwell_param)
            kernel=deltav*pi*(fp(j,iap)+fp(l,iap))**2
          endif
!
          delta_r = fp(j,inpswarm)*kernel - &
                    r_ik_mat(i,j)               ! change in matrix element
          r_ik_mat(i,j) = fp(j,inpswarm)*kernel ! rate matrix element
!  Update rate and cumulative function
          r_i_tot(i) = r_i_tot(i) + delta_r  ! total rate for rep. particle i
!  Only need to update cumulative function from element j and onwards
          cum_func_sec_ik(i, j:) = cum_func_sec_ik(i, j:) + &
                                   delta_r
        endif
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
!  Physical kernels
        elseif (lconstant_deltav) then    ! constant relative velocity
          kernel=deltav*pi*(fp(j,iap)+fp(l,iap))**2
        elseif (lmaxwell_deltav) then  ! maxwellian relative velocity
          call particles_coag_maxwell(deltav, maxwell_param)
          kernel=deltav*pi*(fp(j,iap)+fp(l,iap))**2
        endif
!
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
    subroutine particles_coag_maxwell_johnk(dv, alpha)
!
!  Generate a Maxwell-Boltzmann distributed relative velocity dv with
!  parameter alpha through the Johnk's algorithm.
!
!  For particles following the gas: alpha = sqrt(k_B*T/m)
!  Mean velocity = 2*sqrt(2/pi)*alpha
!
!  14-jan-13/KWJ: coded
!
      use General, only: random_number_wrapper
!
      real :: dv, alpha
      real :: r, w, w1, w2
!
      intent (in) :: alpha
      intent (out) :: dv
!
      call random_number_wrapper(r)
      dv = -log(r)
!
      do while (.true.)
        call random_number_wrapper(r)
        w1 = r*r
        call random_number_wrapper(r)
        w2 = r*r
        w = w1 + w2
        if(w <= 1.) exit
      enddo
!
      call random_number_wrapper(r)
      dv = dv - w1/w*log(r)
!
      dv = alpha*sqrt(2.*dv)
!
    endsubroutine particles_coag_maxwell_johnk
!***********************************************************************
    subroutine particles_coag_maxwell(dv, alpha)
!
!  Generate a Maxwell-Boltzmann distributed relative velocity dv with
!  parameter alpha through the algorithm described in Mohamed (2011).
!
!  For particles following the gas: alpha = sqrt(k_B*T/m)
!  Mean velocity = 2*sqrt(2/pi)*alpha
!
!  14-jan-13/KWJ: coded
!
      use General, only: random_number_wrapper
!
      real :: dv, alpha
      real :: y, r, tmp1, tmp2
!     For findig random nr.
!     C = 0.5, y0 = 1, k = 2/(C*sqrt(pi))*y0^(1/2)*e^(-C*y0)
!     k = 4/sqrt(e*pi) ~ 1.37
!     g = 2/(k*C*sqrt(pi)) = sqrt(e) ~ 1.649
!     real :: g = 1.648721271
      real :: e = 2.718281828
!     Only really interested in g^2 = e
!
      intent (in) :: alpha
      intent (out) :: dv
!
!
      do while (.true.)
        call random_number_wrapper(r)
        y = -2.*log(r)
        tmp1 = e*y*r*r
        call random_number_wrapper(r)
        tmp2 = r*r
        if(tmp1.ge.tmp2) exit
      enddo
!
      dv = alpha*sqrt(2.*y)
!
    endsubroutine particles_coag_maxwell
!***********************************************************************
    subroutine read_particles_coag_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_coag_run_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_coag_run_pars
!***********************************************************************
    subroutine write_particles_coag_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_coag_run_pars)
!
    endsubroutine write_particles_coag_run_pars
!***********************************************************************
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
        idiag_ncoagpm=0; idiag_ncoagpartpm=0; idiag_dt1_coag_par=0
        idiag_k100_100=0; idiag_k100_80=0; idiag_k100_60=0; idiag_k100_50=0
        idiag_k100_40=0; idiag_k100_30=0; idiag_k100_20=0; idiag_k100_10=0
        idiag_k80_80=0; idiag_k80_60=0; idiag_k80_50=0; idiag_k80_40=0
        idiag_k80_30=0; idiag_k80_20=0; idiag_k80_10=0
        idiag_k60_60=0; idiag_k60_50=0; idiag_k60_40=0
        idiag_k60_30=0; idiag_k60_20=0; idiag_k60_10=0
        idiag_k50_50=0; idiag_k50_40=0; idiag_k50_30=0
        idiag_k50_20=0; idiag_k50_10=0; idiag_k40_40=0
        idiag_k40_30=0; idiag_k40_20=0; idiag_k40_10=0
        idiag_k30_30=0; idiag_k30_20=0; idiag_k30_10=0
        idiag_k20_20=0; idiag_k20_10=0; idiag_k10_10=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ncoagpm',idiag_ncoagpm)
        call parse_name(iname,cname(iname),cform(iname), &
            'ncoagpartpm',idiag_ncoagpartpm)
        call parse_name(iname,cname(iname),cform(iname),'dt1_coag_par',idiag_dt1_coag_par)
        call parse_name(iname,cname(iname),cform(iname),'k100_100',idiag_k100_100)
        call parse_name(iname,cname(iname),cform(iname),'k100_80',idiag_k100_80)
        call parse_name(iname,cname(iname),cform(iname),'k100_60',idiag_k100_60)
        call parse_name(iname,cname(iname),cform(iname),'k100_50',idiag_k100_50)
        call parse_name(iname,cname(iname),cform(iname),'k100_40',idiag_k100_40)
        call parse_name(iname,cname(iname),cform(iname),'k100_30',idiag_k100_30)
        call parse_name(iname,cname(iname),cform(iname),'k100_20',idiag_k100_20)
        call parse_name(iname,cname(iname),cform(iname),'k100_10',idiag_k100_10)
        call parse_name(iname,cname(iname),cform(iname),'k80_80',idiag_k80_80)
        call parse_name(iname,cname(iname),cform(iname),'k80_60',idiag_k80_60)
        call parse_name(iname,cname(iname),cform(iname),'k80_50',idiag_k80_50)
        call parse_name(iname,cname(iname),cform(iname),'k80_40',idiag_k80_40)
        call parse_name(iname,cname(iname),cform(iname),'k80_30',idiag_k80_30)
        call parse_name(iname,cname(iname),cform(iname),'k80_20',idiag_k80_20)
        call parse_name(iname,cname(iname),cform(iname),'k80_10',idiag_k80_10)
        call parse_name(iname,cname(iname),cform(iname),'k60_60',idiag_k60_60)
        call parse_name(iname,cname(iname),cform(iname),'k60_50',idiag_k60_50)
        call parse_name(iname,cname(iname),cform(iname),'k60_40',idiag_k60_40)
        call parse_name(iname,cname(iname),cform(iname),'k60_30',idiag_k60_30)
        call parse_name(iname,cname(iname),cform(iname),'k60_20',idiag_k60_20)
        call parse_name(iname,cname(iname),cform(iname),'k60_10',idiag_k60_10)
        call parse_name(iname,cname(iname),cform(iname),'k60_10',idiag_k60_10)
        call parse_name(iname,cname(iname),cform(iname),'k50_50',idiag_k50_50)
        call parse_name(iname,cname(iname),cform(iname),'k50_40',idiag_k50_40)
        call parse_name(iname,cname(iname),cform(iname),'k50_30',idiag_k50_30)
        call parse_name(iname,cname(iname),cform(iname),'k50_20',idiag_k50_20)
        call parse_name(iname,cname(iname),cform(iname),'k50_10',idiag_k50_10)
        call parse_name(iname,cname(iname),cform(iname),'k40_40',idiag_k40_40)
        call parse_name(iname,cname(iname),cform(iname),'k40_30',idiag_k40_30)
        call parse_name(iname,cname(iname),cform(iname),'k40_20',idiag_k40_20)
        call parse_name(iname,cname(iname),cform(iname),'k40_10',idiag_k40_10)
        call parse_name(iname,cname(iname),cform(iname),'k30_30',idiag_k30_30)
        call parse_name(iname,cname(iname),cform(iname),'k30_20',idiag_k30_20)
        call parse_name(iname,cname(iname),cform(iname),'k30_10',idiag_k30_10)
        call parse_name(iname,cname(iname),cform(iname),'k20_20',idiag_k20_20)
        call parse_name(iname,cname(iname),cform(iname),'k20_10',idiag_k20_10)
        call parse_name(iname,cname(iname),cform(iname),'k10_10',idiag_k10_10)
      enddo
!
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_particles_coagulation
!***********************************************************************
endmodule Particles_coagulation
