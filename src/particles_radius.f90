! $Id: particles_radius.f90,v 1.2 2005-08-22 14:03:12 ajohan Exp $
!
!  This module takes care of everything related to particle radius.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 1
! CPARAM logical, parameter :: lparticles_radius=.true.
!
!***************************************************************
module Particles_radius

  use Cdata
  use Particles_cdata
  use Particles_sub
  use Messages

  implicit none

  include 'particles_radius.h'

  real :: ap0=0.0, rhops=1.0e10
  character (len=labellen), dimension(ninit) :: initap='nothing'

  namelist /particles_radius_init_pars/ &
      initap, ap0, rhops

  namelist /particles_radius_run_pars/ &
      rhops

  integer :: idiag_apm=0, idiag_ap2m=0

  contains

!***********************************************************************
    subroutine register_particles_radius()
!
!  Set up indices for access to the fp and dfp arrays
!
!  22-aug-05/anders: coded
!
      use Mpicomm, only: stop_it
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_particles_radius: called twice')
      first = .false.
!
      if (lroot) call cvs_id( &
           "$Id: particles_radius.f90,v 1.2 2005-08-22 14:03:12 ajohan Exp $")
!
!  Indix for particle radius.
!
      iap=npvar+1
!
!  Increase npvar accordingly.
!
      npvar=npvar+1
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call stop_it('register_particles: npvar > mpvar')
      endif
!
    endsubroutine register_particles_radius
!***********************************************************************
    subroutine initialize_particles_radius(lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  22-aug-05/anders: coded
!
      logical :: lstarting
!
    endsubroutine initialize_particles_radius
!***********************************************************************
    subroutine init_particles_radius(f,fp)
!
!  Initial radius of particles.
!
!  22-aug-05/anders: coded
!
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mpar_loc,mpvar) :: fp
!
      integer :: j
!
      do j=1,ninit
        
        select case(initap(j))

        case('nothing')
          if (lroot.and.j==1) print*, 'init_particles_radius: nothing'

        case('constant')
          if (lroot) print*, 'init_particles_radius: constant radius'
          fp(1:npar_loc,iap)=ap0

        endselect

      enddo
!
    endsubroutine init_particles_radius
!***********************************************************************
    subroutine dap_dt(f,fp,dfp)
!
!  Evolution of particle radius.
!
!  22-aug-05/anders: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mpar_loc,mpvar) :: fp, dfp
!
      real, dimension(1) :: rho_ip
      integer :: k
!
      intent (in) :: f, fp
      intent (out) :: dfp
!
      do k=1,npar_loc
        call interpolate_3d_1st(f,ilnrho,ilnrho,fp(k,ixp:izp),rho_ip,ipar(k))
        if (.not. ldensity_nolog) rho_ip=exp(rho_ip)
        dfp(k,iap) = dfp(k,iap) + &
          0.25*sqrt(fp(k,ivpx)**2+fp(k,ivpy)**2+fp(k,ivpz)**2)*rho_ip(1)/rhops
      enddo
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_apm/=0)  call sum_par_name(fp(1:npar_loc,iap),idiag_apm)
        if (idiag_ap2m/=0)  call sum_par_name(fp(1:npar_loc,iap)**2,idiag_ap2m)
      endif
!
    endsubroutine dap_dt
!***********************************************************************
    subroutine read_particles_radius_init_pars(unit,iostat)
!    
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_radius_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_radius_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_radius_init_pars
!***********************************************************************
    subroutine write_particles_radius_init_pars(unit)
!    
      integer, intent (in) :: unit
!
      write(unit,NML=particles_radius_init_pars)
!
    endsubroutine write_particles_radius_init_pars
!***********************************************************************
    subroutine read_particles_radius_run_pars(unit,iostat)
!    
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_radius_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_radius_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_radius_run_pars
!***********************************************************************
    subroutine write_particles_radius_run_pars(unit)
!    
      integer, intent (in) :: unit
!
      write(unit,NML=particles_radius_run_pars)
!
    endsubroutine write_particles_radius_run_pars
!***********************************************************************
    subroutine rprint_particles_radius(lreset,lwrite)
!   
!  Read and register print parameters relevant for particles radius.
!
!  22-aug-05/anders: coded
!
      use Cdata
      use Sub, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
      logical :: lwr
!
!  Write information to index.pro
! 
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
      if (lwr) write(3,*) 'iap=', iap
!
!  Reset everything in case of reset
!
      if (lreset) then
        idiag_apm=0; idiag_ap2m=0
      endif
!
!  Run through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) &
          print*, 'rprint_particles_radius: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'apm',idiag_apm)
        call parse_name(iname,cname(iname),cform(iname),'ap2m',idiag_ap2m)
      enddo
!
    endsubroutine rprint_particles_radius
!***********************************************************************

endmodule Particles_radius
