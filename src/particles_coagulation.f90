! $Id: particles_collisions.f90 13730 2010-04-23 15:05:00Z sven.bingert $
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
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'particles_coagulation.h'
!
  real :: cdtpcoag=0.2, cdtpcoag1=5.0
  real :: kernel_cst=1.0, kernel_lin=1.0, kernel_pro=1.0
  real :: four_pi_rhopmat_over_three2=0.0
  logical :: lcoag_simultaneous=.false., lnoselfcollision=.true.
  logical :: lshear_in_vp=.true.
  logical :: lkernel_test=.false., lconstant_kernel_test=.false.
  logical :: llinear_kernel_test=.false., lproduct_kernel_test=.false.
!
  integer :: idiag_ncoagpm=0, idiag_ncoagpartpm=0
!
  namelist /particles_coag_run_pars/ &
      cdtpcoag, lcoag_simultaneous, lshear_in_vp, lconstant_kernel_test, &
      kernel_cst, llinear_kernel_test, kernel_lin, lproduct_kernel_test, &
      kernel_pro, lnoselfcollision
!
  contains
!***********************************************************************
    subroutine initialize_particles_coag(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-10/anders: coded
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
!  Precalculate inverse of coagulation time-step parameter.
!
      cdtpcoag1=1/cdtpcoag
!
!  Short hand for any kernel test.
!
      lkernel_test=lconstant_kernel_test.or.llinear_kernel_test.or. &
          lproduct_kernel_test
!
!  Squared volume factor needed for product kernel test.
!
      four_pi_rhopmat_over_three2=four_pi_rhopmat_over_three**2
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
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: xpj, xpk, vpj, vpk
      real :: deltavjk, dt1_coag_par, kernel
      real :: npswarmj, npswarmk
      integer :: j, k, l
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
                  deltavjk=sqrt(sum((vpk-vpj)**2))
                  if (sum((vpk-vpj)*(xpk-xpj))<0.0) then
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
      real :: mpsma, mpbig, npsma, npbig, npnew, mpnew, apnew
      real :: rhopsma, rhopbig
      integer :: l, j, k, ncoll, ncoll_par, npart_par
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
!  Change the particle size to the new size, but keep the total mass in the
!  particle swarm the same.
!
!  A representative particle k colliding with a swarm of larger particles j
!  simply obtains the new mass m_k -> m_k + m_j.
!
!  A representative particle k colliding with a swarm of smaller particles j
!  obtains the new mass m_k -> 2*m_k.
!
                      if (lparticles_number) then
                        if (lcoag_simultaneous) then
                          if (fp(j,iap) < fp(k,iap)) then
                            mpsma = four_pi_rhopmat_over_three*fp(j,iap)**3
                            npsma = fp(j,inpswarm)
                            mpbig = four_pi_rhopmat_over_three*fp(k,iap)**3
                            npbig = fp(k,inpswarm)
                          else
                            mpsma = four_pi_rhopmat_over_three*fp(k,iap)**3
                            npsma = fp(k,inpswarm)
                            mpbig = four_pi_rhopmat_over_three*fp(j,iap)**3
                            npbig = fp(j,inpswarm)
                          endif
                          rhopsma=mpsma*npsma
                          rhopbig=mpbig*npbig
                          mpnew=mpbig+rhopsma/npbig
                          apnew=(mpnew/four_pi_rhopmat_over_three)**(1.0/3.0)
                          npnew=0.5*(rhopsma+rhopbig)/mpnew
                          if (npnew*dx*dy*dz<1.0) then
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
                        else
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
                      if (lparticles_number) then
                        npswarmk=fp(k,inpswarm)
                      else
                        npswarmk=rhop_swarm/(four_pi_rhopmat_over_three*fp(k,iap)**3)
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
