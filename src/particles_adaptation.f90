! $Id: particles_coagulation.f90 19828 2012-11-27 09:58:06Z kalle.jansson.89 $
!
!  This modules takes care of adapting the number of particles in a grid cell
!  to a desired value. This module is based on an original idea by Jacob Trier
!  Frederiksen and was developed by Anders Johansen and Chao-Chin Yang.
!
!  EXPERIMENTAL MODULE - PLEASE DO NOT USE WITHOUT CONTACTING THE AUTHORS
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_adaptation= .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Particles_adaptation
!
  use Cdata
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_sub
!
  implicit none
!
  include 'particles_adaptation.h'
!
  integer :: npar_target=100, npar_deviation=50
  character (len=labellen) :: adaptation_method='random'
!
  namelist /particles_adapt_run_pars/ &
      npar_target, npar_deviation, adaptation_method
!
  contains
!***********************************************************************
    subroutine initialize_particles_adaptation(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  03-apr-13/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, intent(in) :: lstarting
!
!  Report fatal error if Particle_mass module not used.
!
      if (.not.lparticles_mass) &
          call fatal_error('initialize_particles_adaptation', &
          'must use Particles_mass module for particle adaptation')
!
!  We must be flexible about the particle number.
!
      if (mpar_loc<2*npar_loc) &
          call fatal_error_local('initialize_particles_adaptation', &
          'must have mpar_loc > 2*npar_loc for particle adaptation')
      call fatal_error_local_collect
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_adaptation
!***********************************************************************
    subroutine particles_adaptation_pencils(f,fp,dfp,ipar,ineargrid)
!
!  Adapt the number of particles in each grid cell to a desired value, under
!  conservation of mass and momentum.
!
!  03-apr-13/anders: coded
!
      use General, only: random_number_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc) :: ipar
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real :: new_weight, vpxm, vpym, vpzm, vpxrms, vpyrms, vpzrms
      real :: r, p
      integer, dimension(nx) :: np, k1_l, k2_l
      integer :: k, l, ix, iy, iz, ix1
!
!  Do particle adaptation pencil by pencil.
!
      do imn=ny*nz,1,-1
        if (npar_imn(imn)/=0) then
          iy=mm(imn)
          iz=nn(imn)
!
!  Sort particles by the grid cell within each pencil.
!
          np=0
          do k=k1_imn(imn),k2_imn(imn)
            ix=ineargrid(k,1)
            ix1=ix-nghost
            np(ix1)=np(ix1)+1
          enddo
!
!  Create the beginning index of particles in each cell.
!
          k1_l(1)=k1_imn(imn)
          do ix1=2,nx
            k1_l(ix1)=k1_l(ix1-1)+np(ix1-1)
          enddo
!
!  Place sorted particles temporarily in the dfp array.
!
          k2_l=k1_l
          do k=k1_imn(imn),k2_imn(imn)
            ix=ineargrid(k,1)
            ix1=ix-nghost
            dfp(k2_l(ix1),:)=fp(k,:)
            k2_l(ix1)=k2_l(ix1)+1
          enddo
          k2_l=k1_l+np-1
!
!  Put sorted particles back into the fp array.
!
          fp(k1_imn(imn):k2_imn(imn),:)=dfp(k1_imn(imn):k2_imn(imn),:)
!
!  Particle adaptation can be done using several methods. The simplest one is
!  'random' which removes particles randomly and creates particles with random
!  positions and velocities.
!
          select case (adaptation_method)
!
          case ('random')
!
!  We need to count backwards in the array of sorted particles in order to be
!  able to remove particles without damaging the particle index arrays k1_l
!  and k2_l.
!
            do ix=l2,l1,-1
              ix1=ix-nghost
!
!  Destroy particles if there are too many.
!
              if (np(ix1)>npar_target+npar_deviation) then 
                new_weight=sum(fp(k1_l(ix1):k2_l(ix1),irhopswarm))/npar_target
                do k=k2_l(ix1),k1_l(ix1)+npar_target,-1
                  call remove_particle(fp,ipar,k)
                enddo
                fp(k1_l(ix1):k1_l(ix1)+npar_target-1,irhopswarm)=new_weight
              endif
!
!  Create particles if there are too few.
!
              if (np(ix1)<npar_target-npar_deviation) then 
                new_weight=sum(fp(k1_l(ix1):k2_l(ix1),irhopswarm))/npar_target
                fp(k1_l(ix1):k2_l(ix1),irhopswarm)=new_weight
                vpxm=sum(fp(k1_l(ix1):k2_l(ix1),ivpx))/np(ix1)
                vpym=sum(fp(k1_l(ix1):k2_l(ix1),ivpy))/np(ix1)
                vpzm=sum(fp(k1_l(ix1):k2_l(ix1),ivpz))/np(ix1)
                vpxrms=sqrt(sum((fp(k1_l(ix1):k2_l(ix1),ivpx)-vpxm)**2)/np(ix1))
                vpyrms=sqrt(sum((fp(k1_l(ix1):k2_l(ix1),ivpy)-vpym)**2)/np(ix1))
                vpzrms=sqrt(sum((fp(k1_l(ix1):k2_l(ix1),ivpz)-vpzm)**2)/np(ix1))
                do k=npar_loc+1,npar_loc+npar_target-np(ix1)
                  call random_number_wrapper(fp(k,ixp))
                  call random_number_wrapper(fp(k,iyp))
                  call random_number_wrapper(fp(k,izp))
                  fp(k,ixp)=x(ix)+(fp(k,ixp)-0.5)*dx
                  fp(k,iyp)=y(iy)+(fp(k,iyp)-0.5)*dy
                  fp(k,izp)=z(iz)+(fp(k,izp)-0.5)*dz
                  call random_number_wrapper(r)
                  call random_number_wrapper(p)
                  fp(k,ivpx)=vpxm+vpxrms*sqrt(-2*log(r))*sin(2*pi*p)
                  call random_number_wrapper(r)
                  call random_number_wrapper(p)
                  fp(k,ivpy)=vpym+vpyrms*sqrt(-2*log(r))*sin(2*pi*p)
                  call random_number_wrapper(r)
                  call random_number_wrapper(p)
                  fp(k,ivpz)=vpzm+vpzrms*sqrt(-2*log(r))*sin(2*pi*p)
                  fp(k,irhopswarm)=new_weight
                enddo
                npar_loc=npar_loc+npar_target-np(ix1)
              endif
            enddo
!
          case ('LBG') ! to be implemented
!
          endselect
!
       endif
!
      enddo
!
    endsubroutine particles_adaptation_pencils
!***********************************************************************
    subroutine read_particles_adapt_run_pars(unit,iostat)
!
!  Read run parameters from run.in.
!
!  03-apr-13/anders: adapted
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,nml=particles_adapt_run_pars,err=99,iostat=iostat)
      else
        read(unit,nml=particles_adapt_run_pars,err=99)
      endif
!
99    return
!
    endsubroutine read_particles_adapt_run_pars
!***********************************************************************
    subroutine write_particles_adapt_run_pars(unit)
!
!  Write run parameters to param.nml.
!
!  03-apr-13/anders: adapted
!
      integer, intent(in) :: unit
!
      write(unit,NML=particles_adapt_run_pars)
!
    endsubroutine write_particles_adapt_run_pars
!*******************************************************************
    subroutine rprint_particles_adaptation(lreset,lwrite)
!
!  Read and register diagnostic parameters.
!
!  03-apr-13/anders: adapted
!
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_particles_adaptation
!***********************************************************************
endmodule Particles_adaptation
