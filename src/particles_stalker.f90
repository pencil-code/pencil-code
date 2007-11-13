! $Id: particles_stalker.f90,v 1.1 2007-11-13 13:44:43 ajohan Exp $
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

  use Cdata
  use Particles_cdata
  use Particles_sub
  use Messages

  implicit none

  include 'particles_stalker.h'

  real :: dstalk=0.1, tstalk=0.0
  integer :: iscratch=0, nout=0
  logical :: linterpolate_cic=.false., linterpolate_tsc=.true.
  logical :: lstalk_xx=.true., lstalk_vv=.true.
  logical :: lstalk_uu=.true., lstalk_guu=.false.
  logical :: lstalk_rho=.true., lstalk_grho=.false.

  namelist /particles_stalker_init_pars/ &
      dstalk, linterpolate_cic, linterpolate_tsc, &
      lstalk_xx, lstalk_vv, lstalk_uu, lstalk_guu, lstalk_rho, lstalk_grho

  namelist /particles_stalker_run_pars/ &
      dstalk, linterpolate_cic, linterpolate_tsc

  contains

!***********************************************************************
    subroutine initialize_particles_stalker(lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  13-nov-07/anders: coded
!
      use FArrayManager
!
      logical :: lstarting
!
!  Need scratch slot in f array to interpolate derived variables.
!
      call farray_acquire_scratch_area('scratch',iscratch)
!
!  Turn off stalking if physics not selected.
!
      if (ivpx==0)        lstalk_vv=.false.
      if (.not. lhydro)   lstalk_uu=.false.
      if (.not. ldensity) lstalk_rho=.false.
      if (.not. ldensity) lstalk_grho=.false.
!
!  Write information on which variables are stalked to file.
!
      if (lroot) then
        open(1,file=trim(datadir)//'/particles_stalker_header.dat', &
            status='unknown')
          if (lstalk_xx)   write(1,'(A)',advance='no') 'xp,yp,zp,'
          if (lstalk_vv)   write(1,'(A)',advance='no') 'vpx,vpy,vpz,'
          if (lstalk_uu)   write(1,'(A)',advance='no') 'ux,uy,uz,'
          if (lstalk_guu)  write(1,'(A)',advance='no') 'duxdx,duxdy,duxdz,'
          if (lstalk_guu)  write(1,'(A)',advance='no') 'duydx,duydy,duydz,'
          if (lstalk_guu)  write(1,'(A)',advance='no') 'duzdx,duzdy,duzdz,'
          if (lstalk_rho)  write(1,'(A)',advance='no') 'rho,'
          if (lstalk_grho) write(1,'(A)',advance='no') 'drhodx,drhody,drhodz,'
        close (1)
      endif
!
!  Read time of next stalking from file.
!
      if (.not. lstarting) then
        open(1,file=trim(datadir)//'/tstalk.dat',form='formatted',status='unknown')
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
    endsubroutine initialize_particles_stalker
!***********************************************************************
    subroutine particles_stalker_sub(f,fp,ineargrid)
!
!  Find local state of the gas at the position of a limited number of
!  particles. This information is saved to a file and can be read in
!  with IDL.
!
!  13-nov-07/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (npar_stalk) :: xp, yp, zp, vpx, vpy, vpz, ux, uy, uz
      real, dimension (npar_stalk) :: rho, drhodx, drhody, drhodz
      real, dimension (npar_stalk) :: duxdx, duxdy, duxdz
      real, dimension (npar_stalk) :: duydx, duydy, duydz
      real, dimension (npar_stalk) :: duzdx, duzdy, duzdz
      real, dimension (100) :: values
      real, dimension (3) :: uu_loc
      integer, dimension (npar_stalk) :: k_stalk
      integer :: i, k, ix0, iy0, iz0, npar_stalk_loc, ivalue
!
!  Only stalk particles every dstalk.
!
      if (t<tstalk) then
!  Do nothing.
      else
!
!  Find out how many stalked particles are at the local processor.
!
        npar_stalk_loc=0
        do k=1,npar_loc
          if (ipar(k)<=npar_stalk) then
            npar_stalk_loc=npar_stalk_loc+1 
            k_stalk(npar_stalk_loc)=k
          endif
        enddo
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
!  Local gas velocity.
!
        if (lstalk_uu) then
          call stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,iux,ux)
          call stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,iuy,uy)
          call stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,iuz,uz)
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
          call stalk_gradient(f,fp,k_stalk,npar_stalk_loc,ineargrid,ilnrho,1,drhodx)
          call stalk_gradient(f,fp,k_stalk,npar_stalk_loc,ineargrid,ilnrho,2,drhody)
          call stalk_gradient(f,fp,k_stalk,npar_stalk_loc,ineargrid,ilnrho,3,drhodz)
        endif
!
!  Write information to a file
!
        open(1,file=trim(directory_snap)//'/particles_stalker.dat', &
            form='unformatted',position='append')
!
!  Write the time and the number of stalked particles at this processor.
!
          write(1) t, npar_stalk_loc
!
!  Collect environment information in single array.
!        
          do i=1,npar_stalk_loc
            ivalue=0
            if (lstalk_xx) then
              ivalue=ivalue+1; values(ivalue)=xp(i)
              ivalue=ivalue+1; values(ivalue)=yp(i)
              ivalue=ivalue+1; values(ivalue)=zp(i)
            endif
            if (lstalk_vv) then
              ivalue=ivalue+1; values(ivalue)=vpx(i)
              ivalue=ivalue+1; values(ivalue)=vpy(i)
              ivalue=ivalue+1; values(ivalue)=vpz(i)
            endif
            if (lstalk_uu) then
              ivalue=ivalue+1; values(ivalue)=ux(i)
              ivalue=ivalue+1; values(ivalue)=uy(i)
              ivalue=ivalue+1; values(ivalue)=uz(i)
            endif
            if (lstalk_rho) then
              ivalue=ivalue+1; values(ivalue)=rho(i)
            endif
            if (lstalk_grho) then
              ivalue=ivalue+1; values(ivalue)=drhodx(i)
              ivalue=ivalue+1; values(ivalue)=drhody(i)
              ivalue=ivalue+1; values(ivalue)=drhodz(i)
            endif
!
!  Write to file.
!          
            write(1) ipar(k_stalk(i)), values(1:ivalue)
          enddo
        close (1)
!
!  Next stalking time is dstalk later.
!
        tstalk=tstalk+dstalk
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
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (npar_stalk) :: k_stalk
      integer :: npar_stalk_loc
      integer, dimension (mpar_loc,3) :: ineargrid
      integer :: ivar
      real, dimension (npar_stalk) :: value
!
      real, dimension (1) :: value_loc
      integer :: i, k, ix0, iy0, iz0
!      
      do i=1,npar_stalk_loc
        k=k_stalk(i)
        ix0=ineargrid(k_stalk(i),1)
        iy0=ineargrid(k_stalk(i),2)
        iz0=ineargrid(k_stalk(i),3)
!
!  Interpolation is either zeroth, first or second order spline interpolation.
!        
        if (lparticlemesh_cic) then
          call interpolate_linear( &
              f,ivar,ivar,fp(k,ixp:izp),value_loc,ineargrid(k,:),ipar(k))
        elseif (lparticlemesh_tsc) then
          call interpolate_quadratic_spline( &
              f,ivar,ivar,fp(k,ixp:izp),value_loc,ineargrid(k,:),ipar(k))
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
      use Deriv, only: der
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
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
!
!  Now that the derivative is stored in the f array, we can use the usual
!  subroutine to find the local state of the gas at the positions of the
!  particles.
!
      call stalk_variable(f,fp,k_stalk,npar_stalk_loc,ineargrid,iscratch,value)
!   
    endsubroutine stalk_gradient
!***********************************************************************
    subroutine read_particles_stalker_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_stalker_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_stalker_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_stalker_init_pars
!***********************************************************************
    subroutine write_particles_stalker_init_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_stalker_init_pars)
!
    endsubroutine write_particles_stalker_init_pars
!***********************************************************************
    subroutine read_particles_stalker_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_stalker_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_stalker_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_stalker_run_pars
!***********************************************************************
    subroutine write_particles_stalker_run_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_stalker_run_pars)
!
    endsubroutine write_particles_stalker_run_pars
!***********************************************************************
endmodule Particles_stalker
