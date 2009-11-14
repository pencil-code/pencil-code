! $Id: particles_viscosity.f90 10506 2009-03-19 12:42:20Z ajohan@strw.leidenuniv.nl $
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
  use Particles_sub
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'particles_collisions.h'
!
  real :: lambda_mfp_single=1.0, coeff_restitution=1.0
  integer :: ncoll_max_par=-1, npart_max_par=-1
  logical :: lcollision_random_angle=.false., lcollision_big_ball=.false.
  logical :: lshear_in_vp=.true.
  character (len=labellen) :: icoll='random-angle'
!
  integer :: idiag_ncollpm=0, idiag_npartpm=0
!
  namelist /particles_coll_run_pars/ &
      lambda_mfp_single, coeff_restitution, icoll, lshear_in_vp, &
      ncoll_max_par, npart_max_par
!
  contains
!***********************************************************************
    subroutine initialize_particles_collisions(f,lstarting)
!
!  07-oct-08/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, intent(in) :: lstarting
!
      select case(icoll)
      case ('random-angle'); lcollision_random_angle=.true.
      case ('big-ball');     lcollision_big_ball=.true.
      case default
        if (lroot) print*, 'No such value for icoll: ', trim(icoll)
        call fatal_error('initialize_particles_collisions','')
      endselect
!
      allocate(kneighbour(mpar_loc))
!
      if (npart_max_par/=-1 .and. (.not.lrandom_particle_pencils)) then
        if (lroot) then
          print*, 'initialize_particles_collisions: '// &
              'with npart_max_par/=-1 one should set'
          print*, '    lrandom_particle_pencils=T in particles_run_pars'
          print*, '    to avoid artificial sub domains in the collisions'
        endif
        call fatal_error('initialize_particles_collisions','')
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_collisions
!***********************************************************************
    subroutine calc_particles_collisions(fp,ineargrid)
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
      use General
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: xpj, xpk, vpj, vpk
      real :: deltavjk, tau_coll1, prob, r
      integer, dimension (nx) :: np_pencil
      integer :: l, j, k, np_point, ncoll, ncoll_par, npart_par
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
        call shepherd_neighbour(fp,ineargrid,kshepherd,kneighbour)
!
!  Calculate number of particles per grid point. This is only needed in order
!  to limit the number of collision partners. Note that f(:,:,:,inp) is not
!  up-to-date at this time.
!
        if (npart_max_par/=-1) then
          np_pencil=0
          do k=k1_imn(imn),k2_imn(imn)
            np_pencil(ineargrid(k,1)-nghost)=np_pencil(ineargrid(k,1)-nghost)+1
          enddo
        endif
!
        do l=l1,l2
          k=kshepherd(l-nghost)
          if (k>0) then
            np_point=np_pencil(l-nghost)
            do while (k/=0)
              j=k
              npart_par=0
              ncoll_par=0
              do while (.true.)
                j=kneighbour(j)
                if (j==0) then
                  if (npart_max_par/=-1 .and. npart_max_par<np_point) then
                    j=kshepherd(l-nghost)
                  else
                    exit
                  endif
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
                deltavjk=sqrt(sum((vpk-vpj)**2))
!
!  The time-scale for collisions between a representative particle from
!  superparticle k and the particle cluster in superparticle j is
!
!    tau_coll = 1/(n*sigma*dv) = lambda/dv
!
!  where lambda is the mean free path of a particle relative to a single
!  superparticle (this is a constant).
!
                tau_coll1=deltavjk/lambda_mfp_single
!
!  Increase collision rate artificially for fewer collisions.
!
                if (npart_max_par/=-1 .and. npart_max_par<np_point) then
                  tau_coll1=tau_coll1*np_point/npart_max_par
                endif
!
                if (tau_coll1/=0.0) then
!
!  The probability for a collision in this time-step is dt/tau_coll.
!
                  prob=dt*tau_coll1
                  call random_number_wrapper(r)
                  if (r<=prob) then
                    if (lshear .and. lshear_in_vp) then
                      vpk(2)=vpk(2)+qshear*Omega*xpk(1)
                      vpj(2)=vpj(2)+qshear*Omega*xpj(1)
                    endif
                    call particle_collision(xpj,xpk,vpj,vpk,j,k)
                    fp(k,ivpx:ivpz)=vpk
                    fp(j,ivpx:ivpz)=vpj
                    ncoll=ncoll+1
                    ncoll_par=ncoll_par+1
                  endif
                endif
                npart_par=npart_par+1
                if (ncoll_max_par/=-1 .and. ncoll_par==ncoll_max_par) exit
                if (npart_max_par/=-1 .and. npart_par==npart_max_par) exit
              enddo
              k=kneighbour(k)
!
!  Collision diagnostics. Since this subroutine is called in the last sub-
!  time-step, we can not use ldiagnos. Therefore we calculate collision
!  diagnostics in the preceding time-step. This has the side effect that
!
!    a) Collision diagnostics for time-step zero are all zero
!    b) Collision diagnostics are not particle normalized if it1==1
!
              if (mod(it,it1)==0) then
                if (idiag_ncollpm/=0) &
                    call sum_par_name((/float(ncoll_par)/),idiag_ncollpm)
                if (idiag_npartpm/=0) &
                    call sum_par_name((/float(npart_par)/),idiag_npartpm)
              endif
!
            enddo
          endif
        enddo
      enddo
!
    endsubroutine calc_particles_collisions
!***********************************************************************
    subroutine particle_collision(xpj,xpk,vpj,vpk,j,k)
!
!  Calculate collision between two particles.
!
!  13-nov-09/anders: coded
!
      use General
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
!  Dissipate some of the relative energy.
!
      if (coeff_restitution/=1.0) &
          vvkcmnew=vvkcmnew*coeff_restitution
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
        if (coeff_restitution/=1.0) &
            vvkcm_normal=vvkcm_normal*coeff_restitution
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
    endsubroutine particle_collision
!***********************************************************************
    subroutine read_particles_coll_run_pars(unit,iostat)
!
!  Read run parameters from run.in.
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
    endsubroutine read_particles_coll_run_pars
!***********************************************************************
    subroutine write_particles_coll_run_pars(unit)
!
!  Write run parameters to param.nml.
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
        idiag_ncollpm=0; idiag_npartpm=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ncollpm',idiag_ncollpm)
        call parse_name(iname,cname(iname),cform(iname),'npartpm',idiag_npartpm)
      enddo
!
    endsubroutine rprint_particles_collisions
!***********************************************************************
endmodule Particles_collisions
