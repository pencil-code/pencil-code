! $Id: particles_radius.f90,v 1.12 2006-03-29 13:41:46 ajohan Exp $
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

  implicit none

  include 'particles_radius.h'

  real :: ap0=0.0, vthresh_sweepup=0.0
  character (len=labellen), dimension(ninit) :: initap='nothing'

  namelist /particles_radius_init_pars/ &
      initap, ap0, rhops, vthresh_sweepup

  namelist /particles_radius_run_pars/ &
      rhops, vthresh_sweepup

  integer :: idiag_apm=0, idiag_ap2m=0, idiag_apmin=0, idiag_apmax=0
  integer :: idiag_dvp12m=0

  contains

!***********************************************************************
    subroutine register_particles_radius()
!
!  Set up indices for access to the fp and dfp arrays
!
!  22-aug-05/anders: coded
!
      use Messages, only: fatal_error, cvs_id
!
      logical, save :: first=.true.
!
      if (.not. first) call fatal_error('register_particles_radius: called twice','')
      first = .false.
!
      if (lroot) call cvs_id( &
           "$Id: particles_radius.f90,v 1.12 2006-03-29 13:41:46 ajohan Exp $")
!
!  Index for particle radius.
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
        call fatal_error('register_particles: npvar > mpvar','')
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
!  Calculate the number density of bodies within a superparticle.
!
      mp_tilde=4/3.*pi*rhops*ap0**3
      if (lroot) print*, 'initialize_particles_radius: '// &
          'mass per dust grain mp_tilde=', mp_tilde
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
    subroutine dap_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle radius.
!
!  22-aug-05/anders: coded
!
      use Messages, only: fatal_error
      use Particles_number
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: uu
      real :: rho, deltav, cc, np_tilde
      integer :: k, ix0, iy0, iz0
      logical :: lheader, lfirstcall=.true.
!
      intent (in) :: f, fp
      intent (out) :: dfp
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
!
!  Identify module and boundary conditions.
!
      if (lheader) print*,'dap_dt: Calculate dap_dt'
!
!  Increase in particle radius due to sweep-up of small grains in the gas.
!
      do k=1,npar_loc
        ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
!  No interpolation needed here.
        rho=f(ix0,iy0,iz0,ilnrho)
        if (.not. ldensity_nolog) rho=exp(rho)
        uu=f(ix0,iy0,iz0,iux:iuz)
!  Relative speed.
        deltav=sqrt( (fp(k,ivpx)-uu(1))**2 + (fp(k,ivpy)-uu(2))**2 + (fp(k,ivpz)-uu(3))**2 )
!  Allow boulders to sweep up small grains if relative velocity not too high.
        if (deltav<=vthresh_sweepup) then
          if (.not. lpscalar) then
            call fatal_error('dap_dt','must have passive scalar module for sweep-up')
          else
            cc=f(ix0,iy0,iz0,ilncc)
            if (.not. lpscalar_nolog) cc=exp(cc)
!  Radius increase due to sweep-up          
            dfp(k,iap) = dfp(k,iap) + 0.25*deltav*cc*rho/rhops
!
!  Deplete gas of small grains.
!
            call get_nptilde(fp,k,np_tilde)
            if (lpscalar_nolog) then 
              df(ix0,iy0,iz0,ilncc) = df(ix0,iy0,iz0,ilncc) - &
                  np_tilde*pi*fp(k,iap)**2*deltav*cc
            else
              df(ix0,iy0,iz0,ilncc) = df(ix0,iy0,iz0,ilncc) - &
                  np_tilde*pi*fp(k,iap)**2*deltav
            endif
          endif
        endif
!
        if (ldiagnos) then
          if (idiag_dvp12m/=0) call sum_par_name((/deltav/),idiag_dvp12m)
        endif
      enddo
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_apm/=0) call sum_par_name(fp(1:npar_loc,iap),idiag_apm)
        if (idiag_ap2m/=0) call sum_par_name(fp(1:npar_loc,iap)**2,idiag_ap2m,lsqrt=.true.)
        if (idiag_apmin/=0) &
            call max_par_name(-fp(1:npar_loc,iap),idiag_apmin,lneg=.true.)
        if (idiag_apmax/=0) &
            call max_par_name(fp(1:npar_loc,iap),idiag_apmax)
      endif
!
      lfirstcall=.false.
!
    endsubroutine dap_dt
!***********************************************************************
    subroutine read_particles_rad_init_pars(unit,iostat)
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
    endsubroutine read_particles_rad_init_pars
!***********************************************************************
    subroutine write_particles_rad_init_pars(unit)
!    
      integer, intent (in) :: unit
!
      write(unit,NML=particles_radius_init_pars)
!
    endsubroutine write_particles_rad_init_pars
!***********************************************************************
    subroutine read_particles_rad_run_pars(unit,iostat)
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
    endsubroutine read_particles_rad_run_pars
!***********************************************************************
    subroutine write_particles_rad_run_pars(unit)
!    
      integer, intent (in) :: unit
!
      write(unit,NML=particles_radius_run_pars)
!
    endsubroutine write_particles_rad_run_pars
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
        idiag_apm=0; idiag_ap2m=0; idiag_apmin=0; idiag_apmax=0
        idiag_dvp12m=0
      endif
!
!  Run through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) &
          print*, 'rprint_particles_radius: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'apm',idiag_apm)
        call parse_name(iname,cname(iname),cform(iname),'ap2m',idiag_ap2m)
        call parse_name(iname,cname(iname),cform(iname),'apmin',idiag_apmin)
        call parse_name(iname,cname(iname),cform(iname),'apmax',idiag_apmax)
        call parse_name(iname,cname(iname),cform(iname),'dvp12m',idiag_dvp12m)
      enddo
!
    endsubroutine rprint_particles_radius
!***********************************************************************

endmodule Particles_radius
