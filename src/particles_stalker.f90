! $Id$
!
!  This module writes information about the local state of the gas at
!  the positions of a selected number of particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_stalker=.true.
!
! MSCRATCH CONTRIBUTION 1
!
!***************************************************************
module Particles_stalker
!
  use Cdata
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
!
  implicit none
!
  include 'particles_stalker.h'
!
  real :: dstalk=0.1, tstalk=0.0
  integer :: iscratch=0, nout=0, nvar_stalk=0
  logical :: linterpolate_cic=.false., linterpolate_tsc=.true.
  logical :: lstalk_xx=.true., lstalk_vv=.true.
  logical :: lstalk_uu=.true., lstalk_guu=.false.
  logical :: lstalk_rho=.true., lstalk_grho=.false.
  logical :: lstalk_bb=.true., lstalk_ap=.true.
  logical :: lstalk_npswarm=.false., lstalk_rhopswarm=.true., lstalk_potself=.true.
  logical :: lstalk_aps=.true.
  logical :: lstalk_sink_particles=.false.
  logical :: lstalk_relvel=.false.
  logical :: lstalk_gTT=.false.
!
  real, dimension (:), allocatable :: xp, yp, zp, vpx, vpy, vpz, ux, uy, uz
  real, dimension (:), allocatable :: rho, drhodx, drhody, drhodz
  real, dimension (:), allocatable :: duxdx, duxdy, duxdz, duydx, duydy, duydz, duzdx, duzdy, duzdz
  real, dimension (:), allocatable :: gTTx, gTTy, gTTz
  real, dimension (:), allocatable :: bx, by, bz, ap, npswarm, rhopswarm
  real, dimension (:), allocatable :: potself, aps, relvel
!
  namelist /particles_stalker_init_pars/ &
      dstalk, linterpolate_cic, linterpolate_tsc, &
      lstalk_xx, lstalk_vv, lstalk_uu, lstalk_guu, lstalk_rho, lstalk_grho, &
      lstalk_bb, lstalk_ap, lstalk_npswarm, lstalk_rhopswarm, lstalk_potself, &
      lstalk_aps, lstalk_sink_particles, lstalk_relvel,lstalk_gTT
!
  namelist /particles_stalker_run_pars/ &
      dstalk, linterpolate_cic, linterpolate_tsc, lstalk_sink_particles
!
  contains
!***********************************************************************
    subroutine initialize_particles_stalker(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  13-nov-07/anders: coded
!
      use General, only: keep_compiler_quiet
      use FArrayManager
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Stop if no particles are stalked.
!
      if (npar_stalk==0) then
        if (lroot) print*, 'initialize_particles_stalker: npar_stalk is zero - set it in cparam.local'
        call fatal_error('initialize_particles_stalker','')
      endif
!
!  Turn off stalking if physics not selected.
!
      if (ivpx==0)       lstalk_vv=.false.
      if (iap==0)        lstalk_ap=.false.
      if (inpswarm==0)   lstalk_npswarm=.false.
      if (irhopswarm==0) lstalk_rhopswarm=.false.
      if (iaps==0)       lstalk_aps=.false.
      if (iuu==0)        lstalk_uu=.false.
      if (iuu==0)        lstalk_guu=.false.
      if (ilnrho==0)     lstalk_rho=.false.
      if (ilnrho==0)     lstalk_grho=.false.
      if (iaa==0)        lstalk_bb=.false.
      if (ipotself==0)   lstalk_potself=.false.
      if (iuu==0 .or. ivpx==0) lstalk_relvel=.false.
      if (iogTTx==0)     lstalk_gTT=.false.
!
!  Need scratch slot in f array to interpolate derived variables.
!
      if (lstalk_guu .or. lstalk_grho .or. lstalk_bb) &
          call farray_acquire_scratch_area('scratch',iscratch)
!
!  Count the number of variables to be stalked.
!
      nvar_stalk=0
      if (lstalk_xx)        nvar_stalk=nvar_stalk+3
      if (lstalk_vv)        nvar_stalk=nvar_stalk+3
      if (lstalk_ap)        nvar_stalk=nvar_stalk+1
      if (lstalk_npswarm)   nvar_stalk=nvar_stalk+1
      if (lstalk_rhopswarm) nvar_stalk=nvar_stalk+1
      if (lstalk_aps)       nvar_stalk=nvar_stalk+1
      if (lstalk_uu)        nvar_stalk=nvar_stalk+3
      if (lstalk_guu)       nvar_stalk=nvar_stalk+9
      if (lstalk_rho)       nvar_stalk=nvar_stalk+1
      if (lstalk_grho)      nvar_stalk=nvar_stalk+3
      if (lstalk_bb)        nvar_stalk=nvar_stalk+3
      if (lstalk_potself)   nvar_stalk=nvar_stalk+1
      if (lstalk_relvel)    nvar_stalk=nvar_stalk+1
      if (lstalk_gTT)       nvar_stalk=nvar_stalk+3
!
!  Reserve space for particle arrays.
!
      if (lstalk_xx) allocate (xp(npar_stalk), yp(npar_stalk), zp(npar_stalk))
      if (lstalk_vv) allocate (vpx(npar_stalk), vpy(npar_stalk), vpz(npar_stalk))
      if (lstalk_ap) allocate (ap(npar_stalk))
      if (lstalk_npswarm) allocate (npswarm(npar_stalk))
      if (lstalk_rhopswarm) allocate (rhopswarm(npar_stalk))
      if (lstalk_aps) allocate (aps(npar_stalk))
      if (lstalk_uu) allocate (ux(npar_stalk), uy(npar_stalk), uz(npar_stalk))
      if (lstalk_guu) allocate (duxdx(npar_stalk), duxdy(npar_stalk), duxdz(npar_stalk), duydx(npar_stalk), duydy(npar_stalk), &
          duydz(npar_stalk), duzdx(npar_stalk), duzdy(npar_stalk), duzdz(npar_stalk))
      if (lstalk_rho) allocate (rho(npar_stalk))
      if (lstalk_grho) allocate (drhodx(npar_stalk), drhody(npar_stalk), drhodz(npar_stalk))
      if (lstalk_bb) allocate (bx(npar_stalk), by(npar_stalk), bz(npar_stalk))
      if (lstalk_potself) allocate (potself(npar_stalk))
      if (lstalk_relvel) allocate (relvel(npar_stalk))
      if (lstalk_gTT) allocate (gTTx(npar_stalk), gTTy(npar_stalk), gTTz(npar_stalk))
!
!  Write information on which variables are stalked to file.
!
      if (lroot) then
        open(1,file=trim(datadir)//'/particles_stalker_header.dat', &
            status='unknown')
          if (lstalk_xx)        write(1,'(A)',advance='no') 'xp,yp,zp,'
          if (lstalk_vv)        write(1,'(A)',advance='no') 'vpx,vpy,vpz,'
          if (lstalk_ap)        write(1,'(A)',advance='no') 'ap,'
          if (lstalk_npswarm)   write(1,'(A)',advance='no') 'npswarm,'
          if (lstalk_rhopswarm) write(1,'(A)',advance='no') 'rhopswarm,'
          if (lstalk_aps)       write(1,'(A)',advance='no') 'aps,'
          if (lstalk_uu)        write(1,'(A)',advance='no') 'ux,uy,uz,'
          if (lstalk_guu)       write(1,'(A)',advance='no') 'duxdx,duxdy,duxdz,'
          if (lstalk_guu)       write(1,'(A)',advance='no') 'duydx,duydy,duydz,'
          if (lstalk_guu)       write(1,'(A)',advance='no') 'duzdx,duzdy,duzdz,'
          if (lstalk_rho)       write(1,'(A)',advance='no') 'rho,'
          if (lstalk_grho)      write(1,'(A)',advance='no') 'drhodx,drhody,drhodz,'
          if (lstalk_bb)        write(1,'(A)',advance='no') 'bx,by,bz,'
          if (lstalk_potself)   write(1,'(A)',advance='no') 'potself,'
          if (lstalk_relvel)    write(1,'(A)',advance='no') 'relvel,'
          if (lstalk_gTT)       write(1,'(A)',advance='no') 'dTTdx,dTTdy,dTTdz,'
         close (1)
      endif
!
!  Read time of next stalking from file.
!
      if (lrun) then
        open(1,file=trim(datadir)//'/tstalk.dat',form='formatted', &
            status='unknown')
          read(1,*) tstalk, nout
        close(1)
      else
        tstalk=0.0
        if (lroot) then
          open(1,file=trim(datadir)//'/tstalk.dat',form='formatted')
            write(1,'(e14.7,i8)') tstalk, nout
          close(1)
        endif
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_stalker
!***********************************************************************
    subroutine finalize_particles_stalker
!
!  Finalize module.
!
!  01-May-2019/PABourdin: coded
!
      if (lstalk_xx) deallocate (xp, yp, zp)
      if (lstalk_vv) deallocate (vpx, vpy, vpz)
      if (lstalk_ap) deallocate (ap)
      if (lstalk_npswarm) deallocate (npswarm)
      if (lstalk_rhopswarm) deallocate (rhopswarm)
      if (lstalk_aps) deallocate (aps)
      if (lstalk_uu) deallocate (ux, uy, uz)
      if (lstalk_guu) deallocate (duxdx, duxdy, duxdz, duydx, duydy, duydz, duzdx, duzdy, duzdz)
      if (lstalk_rho) deallocate (rho)
      if (lstalk_grho) deallocate (drhodx, drhody, drhodz)
      if (lstalk_bb) deallocate (bx, by, bz)
      if (lstalk_potself) deallocate (potself)
      if (lstalk_relvel) deallocate (relvel)
      if (lstalk_gTT) deallocate (gTTx, gTTy, gTTz)
!
    endsubroutine finalize_particles_stalker
!***********************************************************************
    subroutine particles_stalker_sub(f,fp,ineargrid)
!
!  Find local state of the gas at the position of a limited number of
!  particles. This information is saved to a file and can be read in
!  with IDL.
!
!  13-nov-07/anders: coded
!
      use IO, only: output_stalker_init, output_stalker, output_part_finalize
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      real :: t_sp
!
      integer, dimension (npar_stalk) :: k_stalk
      integer :: i, k, npar_stalk_loc, ivalue
!
!  Only stalk particles every dstalk.
!
      if (t<tstalk) then
!  Do nothing.
      else
        t_sp=t
!
!  Find out how many stalked particles are at the local processor.
!
        npar_stalk_loc=0
        if (lstalk_sink_particles) then
          do k=1,npar_loc
            if (fp(k,iaps)>0.0) then
              npar_stalk_loc=npar_stalk_loc+1
              if (npar_stalk_loc>npar_stalk) then
                print*, 'particles_stalker_sub: too many sink particles '// &
                    'are stalked'
                print*, 'iproc, it, itsub, npar_stalk, npar_stalk_loc=', &
                    it, itsub, npar_stalk, npar_stalk_loc
                call fatal_error_local('particles_stalker_sub','')
              else
                k_stalk(npar_stalk_loc)=k
              endif
            endif
          enddo
          call fatal_error_local_collect()
        else
          do k=1,npar_loc
            if (ipar(k)<=npar_stalk) then
              npar_stalk_loc=npar_stalk_loc+1
              k_stalk(npar_stalk_loc)=k
            endif
          enddo
        endif
!
!  Gather environment information for each stalked particle, starting with the
!  position.
!
        if (lstalk_xx) then
          do i=1,npar_stalk_loc
            xp(i)=fp(k_stalk(i),ixp)
            yp(i)=fp(k_stalk(i),iyp)
            zp(i)=fp(k_stalk(i),izp)
          enddo
        endif
!
!  Velocity vector (only relevant for inertial particles, not for tracers).
!
        if (lstalk_vv) then
          do i=1,npar_stalk_loc
            vpx(i)=fp(k_stalk(i),ivpx)
            vpy(i)=fp(k_stalk(i),ivpy)
            vpz(i)=fp(k_stalk(i),ivpz)
          enddo
        endif
!
!  Particle radius.
!
        if (lstalk_ap) then
          do i=1,npar_stalk_loc
            ap(i)=fp(k_stalk(i),iap)
           enddo
        endif
!
!  Particle number density.
!
        if (lstalk_npswarm) then
          do i=1,npar_stalk_loc
            npswarm(i)=fp(k_stalk(i),inpswarm)
           enddo
        endif
!
!  Particle mass density.
!
        if (lstalk_rhopswarm) then
          do i=1,npar_stalk_loc
            rhopswarm(i)=fp(k_stalk(i),irhopswarm)
           enddo
        endif
!
!  Sink particle radius.
!
        if (lstalk_aps) then
          do i=1,npar_stalk_loc
            aps(i)=fp(k_stalk(i),iaps)
           enddo
        endif
!
!  Local gas velocity.
!
        if (lstalk_uu) then
          call stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,iux,ux)
          call stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,iuy,uy)
          call stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,iuz,uz)
        endif
!
!  Temperature gradient
!  Note that when iogTT is already the gradient of the temperature,
!  calculated and stored for use in another module, so it is stalked
!  like the local velocity
!        
        if (lstalk_gTT) then
          call stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,iogTTx,gTTx)
          call stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,iogTTy,gTTy)
          call stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,iogTTz,gTTz)
        endif        
!
!  Local gradient of gas velocity.
!
        if (lstalk_guu) then
          call stalk_gradient(f,fp,k_stalk,npar_stalk_loc,ineargrid,iux,1,duxdx)
          call stalk_gradient(f,fp,k_stalk,npar_stalk_loc,ineargrid,iux,2,duxdy)
          call stalk_gradient(f,fp,k_stalk,npar_stalk_loc,ineargrid,iux,3,duxdz)
          call stalk_gradient(f,fp,k_stalk,npar_stalk_loc,ineargrid,iuy,1,duydx)
          call stalk_gradient(f,fp,k_stalk,npar_stalk_loc,ineargrid,iuy,2,duydy)
          call stalk_gradient(f,fp,k_stalk,npar_stalk_loc,ineargrid,iuy,3,duydz)
          call stalk_gradient(f,fp,k_stalk,npar_stalk_loc,ineargrid,iuz,1,duzdx)
          call stalk_gradient(f,fp,k_stalk,npar_stalk_loc,ineargrid,iuz,2,duzdy)
          call stalk_gradient(f,fp,k_stalk,npar_stalk_loc,ineargrid,iuz,3,duzdz)
        endif
!
!  Local gas density (logarithmic or not logarithmic).
!
        if (lstalk_rho) then
          call stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,ilnrho,rho)
        endif
!
!  Local density gradient.
!
        if (lstalk_grho) then
          call stalk_gradient(f,fp,k_stalk,npar_stalk_loc,ineargrid,ilnrho,1, &
              drhodx)
          call stalk_gradient(f,fp,k_stalk,npar_stalk_loc,ineargrid,ilnrho,2, &
              drhody)
          call stalk_gradient(f,fp,k_stalk,npar_stalk_loc,ineargrid,ilnrho,3, &
              drhodz)
        endif
!
!  Local magnetic field.
!
        if (lstalk_bb) then
          call stalk_magnetic(f,fp,k_stalk,npar_stalk_loc,ineargrid,bx,by,bz)
        endif
!
!  Local gravitational potential.
!
        if (lstalk_potself) then
          call stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,ipotself, &
              potself)
        endif
!
!  Relative velocity of gas and particle
!
        if (lstalk_relvel) then
          allocate (relvel(npar_stalk))
          do i=1,npar_stalk_loc
            vpx(i)=fp(k_stalk(i),ivpx)
            vpy(i)=fp(k_stalk(i),ivpy)
            vpz(i)=fp(k_stalk(i),ivpz)
          enddo
          call stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,iux,ux)
          call stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,iuy,uy)
          call stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,iuz,uz)
          do i=1,npar_stalk_loc
             relvel(i)=sqrt(((vpx(i)-ux(i))**2)+((vpy(i)-uy(i))**2)+((vpz(i)-uz(i))**2))
          enddo
        endif
!
!  Write information to a file
!
        call output_stalker_init (nvar_stalk, npar_stalk_loc, nout, ipar(k_stalk(1:npar_stalk_loc)))
        if (lstalk_xx) then
          call output_stalker ('xp', npar_stalk, npar_stalk_loc, xp)
          call output_stalker ('yp', npar_stalk, npar_stalk_loc, yp)
          call output_stalker ('zp', npar_stalk, npar_stalk_loc, zp)
        endif
        if (lstalk_vv) then
          call output_stalker ('vpx', npar_stalk, npar_stalk_loc, vpx)
          call output_stalker ('vpy', npar_stalk, npar_stalk_loc, vpy)
          call output_stalker ('vpz', npar_stalk, npar_stalk_loc, vpz)
        endif
        if (lstalk_ap) call output_stalker ('ap', npar_stalk, npar_stalk_loc, ap)
        if (lstalk_npswarm) call output_stalker ('npswarm', npar_stalk, npar_stalk_loc, npswarm)
        if (lstalk_rhopswarm) call output_stalker ('rhopswarm', npar_stalk, npar_stalk_loc, rhopswarm)
        if (lstalk_aps) call output_stalker ('aps', npar_stalk, npar_stalk_loc, aps)
        if (lstalk_uu) then
          call output_stalker ('ux', npar_stalk, npar_stalk_loc, ux)
          call output_stalker ('uy', npar_stalk, npar_stalk_loc, uy)
          call output_stalker ('uz', npar_stalk, npar_stalk_loc, uz)
        endif
        if (lstalk_guu) then
          call output_stalker ('duxdx', npar_stalk, npar_stalk_loc, duxdx)
          call output_stalker ('duxdy', npar_stalk, npar_stalk_loc, duxdy)
          call output_stalker ('duxdz', npar_stalk, npar_stalk_loc, duxdz)
          call output_stalker ('duydx', npar_stalk, npar_stalk_loc, duydx)
          call output_stalker ('duydy', npar_stalk, npar_stalk_loc, duydy)
          call output_stalker ('duydz', npar_stalk, npar_stalk_loc, duydz)
          call output_stalker ('duzdx', npar_stalk, npar_stalk_loc, duzdx)
          call output_stalker ('duzdy', npar_stalk, npar_stalk_loc, duzdy)
          call output_stalker ('duzdz', npar_stalk, npar_stalk_loc, duzdz)
        endif
        if (lstalk_rho) call output_stalker ('rho', npar_stalk, npar_stalk_loc, rho)
        if (lstalk_grho) then
          call output_stalker ('drhodx', npar_stalk, npar_stalk_loc, drhodx)
          call output_stalker ('drhody', npar_stalk, npar_stalk_loc, drhody)
          call output_stalker ('drhodz', npar_stalk, npar_stalk_loc, drhodz)
        endif
        if (lstalk_bb) then
          call output_stalker ('bx', npar_stalk, npar_stalk_loc, bx)
          call output_stalker ('by', npar_stalk, npar_stalk_loc, by)
          call output_stalker ('bz', npar_stalk, npar_stalk_loc, bz)
        endif
        if (lstalk_potself) call output_stalker ('potself', npar_stalk, npar_stalk_loc, potself)
        if (lstalk_relvel) call output_stalker ('relvel', npar_stalk, npar_stalk_loc, relvel)
        if (lstalk_gTT) then
          call output_stalker ('gTTx', npar_stalk, npar_stalk_loc, gTTx)
          call output_stalker ('gTTy', npar_stalk, npar_stalk_loc, gTTy)
          call output_stalker ('gTTz', npar_stalk, npar_stalk_loc, gTTz)
        endif
        call output_part_finalize
!
!  Next stalking time is dstalk later.
!
        do while (tstalk<=t_sp)
          tstalk=tstalk+dstalk
        enddo
        nout=nout+1
        if (lroot) then
          open(1,file=trim(datadir)//'/tstalk.dat',form='formatted')
            write(1,'(e14.7,i8)') tstalk, nout
          close(1)
        endif
      endif
!
    endsubroutine particles_stalker_sub
!***********************************************************************
    subroutine stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,ivar,value)
!
!  Find local state of the gas at the position of a limited number of
!  particles. This subroutines take care of interpolation all simple
!  physical variables, i.e. the ones that are directly present in the f
!  array.
!
!  13-nov-07/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (npar_stalk) :: k_stalk
      integer :: npar_stalk_loc
      integer, dimension (mpar_loc,3) :: ineargrid
      integer :: ivar
      real, dimension (npar_stalk) :: value
!
      real, dimension (1) :: value_loc
      integer :: i, k, ix0, iy0, iz0, iblock
!
      do i=1,npar_stalk_loc
        k=k_stalk(i)
        ix0=ineargrid(k_stalk(i),1)
        iy0=ineargrid(k_stalk(i),2)
        iz0=ineargrid(k_stalk(i),3)
!
        if (lparticles_blocks) then
          iblock=inearblock(k)
        else
          iblock=0
        endif
!
!  Interpolation is either zeroth, first or second order spline interpolation.
!
        if (lparticlemesh_cic) then
          call interpolate_linear( &
              f,ivar,ivar,fp(k,ixp:izp),value_loc,ineargrid(k,:),iblock,ipar(k))
        elseif (lparticlemesh_tsc) then
          call interpolate_quadratic_spline( &
              f,ivar,ivar,fp(k,ixp:izp),value_loc,ineargrid(k,:),iblock,ipar(k))
        else
          value_loc=(/f(ix0,iy0,iz0,ivar)/)
        endif
        value(i)=value_loc(1)
      enddo
!
    endsubroutine stalk_variable
!***********************************************************************
    subroutine stalk_gradient(f,fp,k_stalk,npar_stalk_loc,ineargrid,ivar,ider,value)
!
!  Take simple derivative of physical variable in f and store the result
!  in scratch slot in f. Then find the local state of the gas using the usual
!  subroutine 'stalk_variable'.
!
!  13-nov-07/anders: coded
!
      use Boundcond, only: bc_per_x, bc_per_y, bc_per_z
      use Deriv, only: der
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (npar_stalk) :: k_stalk
      integer :: npar_stalk_loc
      integer, dimension (mpar_loc,3) :: ineargrid
      integer :: ivar, ider
      real, dimension (npar_stalk) :: value
!
      real, dimension (nx) :: der_pencil
!
!  Calculate derivative of variable number ivar with respect to
!  direction number ider.
!
      do n=n1,n2; do m=m1,m2
        call der(f(:,:,:,ivar),der_pencil,ider)
        f(l1:l2,m,n,iscratch)=der_pencil
      enddo; enddo
      call bc_per_x(f,'top',iscratch); call bc_per_x(f,'bot',iscratch)
      call bc_per_y(f,'top',iscratch); call bc_per_y(f,'bot',iscratch)
      call bc_per_z(f,'top',iscratch); call bc_per_z(f,'bot',iscratch)
!
!  Now that the derivative is stored in the f array, we can use the usual
!  subroutine to find the local state of the gas at the positions of the
!  particles.
!
      call stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,iscratch,value)
!
    endsubroutine stalk_gradient
!***********************************************************************
    subroutine stalk_magnetic(f,fp,k_stalk,npar_stalk_loc,ineargrid,bx,by,bz)
!
!  Calculate magnetic field B from vector potential A and interpolate
!  to position of stalked particles.
!
!  14-nov-07/anders: coded
!
      use Boundcond, only: bc_per_x, bc_per_y, bc_per_z
      use Deriv, only: der
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (npar_stalk) :: k_stalk
      integer :: npar_stalk_loc
      integer, dimension (mpar_loc,3) :: ineargrid
      real, dimension (npar_stalk) :: bx, by, bz
!
      real, dimension (nx) :: der1_pencil, der2_pencil
!
      do n=n1,n2; do m=m1,m2
        call der(f,iaa-1+3,der1_pencil,2)
        call der(f,iaa-1+2,der2_pencil,3)
        f(l1:l2,m,n,iscratch)=der1_pencil-der2_pencil
      enddo; enddo
      call bc_per_x(f,'top',iscratch); call bc_per_x(f,'bot',iscratch)
      call bc_per_y(f,'top',iscratch); call bc_per_y(f,'bot',iscratch)
      call bc_per_z(f,'top',iscratch); call bc_per_z(f,'bot',iscratch)
      call stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,iscratch,bx)
!
      do n=n1,n2; do m=m1,m2
        call der(f,iaa-1+1,der1_pencil,3)
        call der(f,iaa-1+3,der2_pencil,1)
        f(l1:l2,m,n,iscratch)=der1_pencil-der2_pencil
      enddo; enddo
      call bc_per_x(f,'top',iscratch); call bc_per_x(f,'bot',iscratch)
      call bc_per_y(f,'top',iscratch); call bc_per_y(f,'bot',iscratch)
      call bc_per_z(f,'top',iscratch); call bc_per_z(f,'bot',iscratch)
      call stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,iscratch,by)
!
      do n=n1,n2; do m=m1,m2
        call der(f,iaa-1+2,der1_pencil,1)
        call der(f,iaa-1+1,der2_pencil,2)
        f(l1:l2,m,n,iscratch)=der1_pencil-der2_pencil
      enddo; enddo
      call bc_per_x(f,'top',iscratch); call bc_per_x(f,'bot',iscratch)
      call bc_per_y(f,'top',iscratch); call bc_per_y(f,'bot',iscratch)
      call bc_per_z(f,'top',iscratch); call bc_per_z(f,'bot',iscratch)
      call stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,iscratch,bz)
!
    endsubroutine stalk_magnetic
!***********************************************************************
!    subroutine stalk_rep(fp,npar_stalk_loc,rep)
!      real, dimension (mpar_loc,mparray) :: fp
!      real, dimension (npar_stalk) :: rep,nu
!      integer :: npar_stalk_loc
!      integer :: k
!
!      do k=1,npar_stalk
!        rep(k)=2.0*fp(k,iap)*sqrt(sum((interp_uu(k,:)-fp(k,ivpx:ivpz))**2))/nu(k)
!      enddo
!    endsubroutine stalk_rep
!***********************************************************************
    subroutine read_pstalker_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_stalker_init_pars, IOSTAT=iostat)
!
    endsubroutine read_pstalker_init_pars
!***********************************************************************
    subroutine write_pstalker_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_stalker_init_pars)
!
    endsubroutine write_pstalker_init_pars
!***********************************************************************
    subroutine read_pstalker_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_stalker_run_pars, IOSTAT=iostat)
!
    endsubroutine read_pstalker_run_pars
!***********************************************************************
    subroutine write_pstalker_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_stalker_run_pars)
!
    endsubroutine write_pstalker_run_pars
!***********************************************************************
endmodule Particles_stalker
