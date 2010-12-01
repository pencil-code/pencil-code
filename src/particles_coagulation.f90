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
  real :: kernel_cst=1.0, cdtpcoag=0.2, cdtpcoag1=5.0
  logical :: lshear_in_vp=.true., lconstant_kernel_test=.false.
!
  integer :: idiag_ncoagpm=0, idiag_ncoagpartpm=0
!
  namelist /particles_coag_run_pars/ &
      lshear_in_vp, lconstant_kernel_test, kernel_cst
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
!  Allocate neighbour array necessary for identifying collisions.
!
      if (.not.allocated(kneighbour)) allocate(kneighbour(mpar_loc))
!
!  Precalculate inverse of coagulation time-step parameter.
!
      cdtpcoag1=1/cdtpcoag
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
      real :: deltavjk, dt1_coag_par
      integer :: j, k, l
!
!  Create list of shepherd and neighbour particles for each grid cell in the
!  current pencil.
!
      call shepherd_neighbour_pencil(fp,ineargrid,kshepherd,kneighbour)
!
      do l=l1,l2
        k=kshepherd(l-nghost)
        if (k>0) then
          do while (k/=0)
            dt1_coag_par=0.0
            j=kshepherd(l-nghost)
            do while (.true.)
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
              if (lconstant_kernel_test) then
                dt1_coag_par=dt1_coag_par+kernel_cst* &
                    min(fp(j,inpswarm),fp(k,inpswarm))
              else
!
!  Only consider collisions between particles approaching each other.
!
                deltavjk=sqrt(sum((vpk-vpj)**2))
                if (sum((vpk-vpj)*(xpk-xpj))<0.0) then
                  dt1_coag_par=dt1_coag_par+ &
                      pi*(fp(k,iap)+fp(k,iap))**2*deltavjk* &
                      min(fp(j,inpswarm),fp(k,inpswarm))
                endif
!
              endif
              j=kneighbour(j)
              if (j==0) exit
            enddo
!
            dt1_max(l-nghost)=max(dt1_max(l-nghost),dt1_coag_par*cdtpcoag1)
!
            k=kneighbour(k)
!
          enddo
        endif
      enddo
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
      real :: lambda_mfp1, deltavjk, tau_coll1, prob, r
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
          k=kshepherd(l-nghost)
          if (k>0) then
            do while (k/=0)
              j=kshepherd(l-nghost)
              npart_par=0
              ncoll_par=0
              do while (.true.)
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
                if ((sum((vpk-vpj)*(xpk-xpj))<0.0).or. &
                     lconstant_kernel_test) then
!
!  Relative particle speed.
!
                  deltavjk=sqrt(sum((vpk-vpj)**2))
!
!  The time-scale for collisions between a representative particle from
!  superparticle k and the particle swarm in superparticle j is
!
!    tau_coll = 1/(n*sigma*dv) = lambda/dv
!
!  where lambda is the mean free path of a particle relative to a single
!  superparticle and sigma is the collisional cross section.
!
                  if (lconstant_kernel_test) then
                    tau_coll1=kernel_cst*min(fp(j,inpswarm),fp(k,inpswarm))
                  else
                    if (lparticles_number) then
                      tau_coll1=deltavjk*pi*(fp(k,iap)+fp(j,iap))**2* &
                          min(fp(j,inpswarm),fp(k,inpswarm))
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
!  Change the particle size to the new size, but keep the total mass in the
!  particle swarm the same.
!
                      if (lparticles_number) then
                        fp(j,iap)=2**(1.0/3.0)*max(fp(k,iap),fp(k,iap))
                        fp(k,iap)=fp(j,iap)
                        fp(j,inpswarm)=0.5*min(fp(k,inpswarm),fp(k,inpswarm))
                        fp(k,inpswarm)=fp(j,inpswarm)
                      endif
!
                      ncoll=ncoll+1
                      ncoll_par=ncoll_par+1
!
                    endif
                  endif
                endif
                j=kneighbour(j)
                if (j==0) exit
              enddo
              k=kneighbour(k)
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
