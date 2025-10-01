! $Id$
!
!  This module takes care of stirring of particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_stirring=.true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Particles_stirring
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_sub
!
  implicit none
!
  include 'particles_stirring.h'
!
  real :: deltat_stir=1.0, deltav_stir=0.0
!
  namelist /particles_stirring_run_pars/ &
      deltat_stir, deltav_stir
!
  contains
!***********************************************************************
    subroutine register_particles_stirring()
!
    endsubroutine register_particles_stirring
!***********************************************************************
    subroutine particle_stirring(fp,ineargrid)
!
!  Particle stirring by random, uncorrelated kicks.
!
!  30-dec-10/anders+michiel: coded
!
      use General
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real :: theta, phi, prob, r
      integer :: k
!
      do k=1,npar_loc
        prob=dt/deltat_stir
        call random_number_wrapper(r)
        if (r<=prob) then
          call random_number_wrapper(theta)
          call random_number_wrapper(phi)
          theta=acos(2*theta-1)
          phi  =phi*2*pi
          fp(k,ivpx)=fp(k,ivpx)+sin(theta)*cos(phi)*deltav_stir
          fp(k,ivpy)=fp(k,ivpy)+sin(theta)*sin(phi)*deltav_stir
          fp(k,ivpz)=fp(k,ivpz)+cos(theta)         *deltav_stir
        endif
      enddo
!
    endsubroutine particle_stirring
!***********************************************************************
    subroutine read_particles_stir_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_stirring_run_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_stir_run_pars
!***********************************************************************
    subroutine write_particles_stir_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_stirring_run_pars)
!
    endsubroutine write_particles_stir_run_pars
!***********************************************************************
endmodule Particles_stirring
