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
!
  namelist /particles_coll_run_pars/ &
      lambda_mfp_single, coeff_restitution
!
  contains
!***********************************************************************
    subroutine initialize_particles_collisions(lstarting)
!
!  07-oct-08/anders: coded
!
      logical, intent(in) :: lstarting
!
      allocate(kneighbour(mpar_loc))
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
      use General
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: vvpbar, vvkcmt, vvkcmtnew
      real :: deltavjk, tau_coll1, prob, r, theta_rot, phi_rot
      integer :: l, j, k
!
!  Pencil loop.
!
      do imn=1,ny*nz
!
!  Create list of shepherd and neighbour particles for each grid cell in the
!  current pencil.
!
        call shepherd_neighbour(fp,ineargrid,kshepherd,kneighbour)
        do l=l1,l2
          k=kshepherd(l-nghost)
          if (k>0) then
            do while (k/=0)
              j=k
              do while (kneighbour(j)/=0)
                j=kneighbour(j)
!
!  Calculate the relative speed of particles j and k.
!
                deltavjk=sqrt( &
                    (fp(k,ivpx)-fp(j,ivpx))**2 + &
                    (fp(k,ivpy)-fp(j,ivpy))**2 + &
                    (fp(k,ivpz)-fp(j,ivpz))**2 )
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
                if (tau_coll1/=0.0) then
!
!  The probability for a collision in this time-step is dt/tau_coll.
!
                  prob=dt*tau_coll1
                  call random_number_wrapper(r)
                  if (r<=prob) then
                    if (ip<=6) then
                      print*, 'calc_particles_collisions: collision between'// &
                          ' superparticles ', ipar(k), ' and ', ipar(j)
                      print*, 'calc_particles_collisions: fp(j,ivpx:ivpz)='
                      print*, '  ', fp(j,ivpx:ivpz)
                      print*, 'calc_particles_collisions: fp(k,ivpx:ivpz)='
                      print*, '  ', fp(k,ivpx:ivpz)
                      print*, 'calc_particles_collisions: Ej, Ek, Ej+Ek='
                      print*, '  ', 0.5*(sum(fp(j,ivpx:ivpz)**2)), &
                          0.5*(sum(fp(k,ivpx:ivpz)**2)), &
                          0.5*(sum(fp(j,ivpx:ivpz)**2))+ &
                          0.5*(sum(fp(k,ivpx:ivpz)**2))
                    endif
!
!  Choose random unit vector direction in COM frame. Here theta_rot is the
!  pole angle and phi_rot is the azimuthal angle. While we can choose phi_rot
!  randomly, theta_rot must be chosen with care so that equal areas on the
!  unit sphere are equally probable.
!    (see http://mathworld.wolfram.com/SpherePointPicking.html)
!
                    call random_number_wrapper(theta_rot)
                    theta_rot=acos(2*theta_rot-1)
                    call random_number_wrapper(phi_rot)
                    phi_rot  =phi_rot  *2*pi
                    vvpbar=0.5*(fp(j,ivpx:ivpz)+fp(k,ivpx:ivpz))
                    vvkcmt=fp(k,ivpx:ivpz)-vvpbar
                    vvkcmtnew(1)=+sin(theta_rot)*cos(phi_rot)
                    vvkcmtnew(2)=+sin(theta_rot)*sin(phi_rot)
                    vvkcmtnew(3)=+cos(theta_rot)
!
!  Multiply unit vector by length of velocity vector.
!
                    vvkcmtnew=vvkcmtnew* &
                        sqrt(vvkcmt(1)**2+vvkcmt(2)**2+vvkcmt(3)**2)
!
!  Dissipate some of the relative energy.
!
                    if (coeff_restitution/=1.0) &
                        vvkcmtnew=vvkcmtnew*coeff_restitution
!
!  Change velocity vectors in normal frame.
!
                    fp(k,ivpx:ivpz)=+vvkcmtnew+vvpbar
                    fp(j,ivpx:ivpz)=-vvkcmtnew+vvpbar
                    if (ip<=6) then
                      print*, 'calc_particles_collisions: Ej2, Ek2, Ej2+Ek2='
                      print*, '  ', 0.5*(sum(fp(j,ivpx:ivpz)**2)), &
                          0.5*(sum(fp(k,ivpx:ivpz)**2)), &
                          0.5*(sum(fp(j,ivpx:ivpz)**2))+ &
                          0.5*(sum(fp(k,ivpx:ivpz)**2))
                    endif
                  endif
                endif
              enddo
              k=kneighbour(k)
            enddo
!  "if (k>0) then"
          endif
!
        enddo
      enddo
!
    endsubroutine calc_particles_collisions
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
!***********************************************************************
endmodule Particles_collisions
