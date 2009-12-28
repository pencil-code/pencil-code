! $Id$
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
!
  use Cdata
  use Messages
  use Particles_cdata
  use Particles_sub
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'particles_radius.h'
!
  real :: vthresh_sweepup=-1.0, deltavp12_floor=0.0
  real, dimension (ninit) :: ap0=0.0
  real :: tstart_sweepup_par=0.0, cdtps=0.2
  logical :: lsweepup_par=.true.
  character (len=labellen), dimension(ninit) :: initap='nothing'
!
  namelist /particles_radius_init_pars/ &
      initap, ap0, rhops, vthresh_sweepup, deltavp12_floor, &
      lsweepup_par, tstart_sweepup_par, cdtps
!
  namelist /particles_radius_run_pars/ &
      rhops, vthresh_sweepup, deltavp12_floor, &
      lsweepup_par, tstart_sweepup_par, cdtps
!
  integer :: idiag_apm=0, idiag_ap2m=0, idiag_apmin=0, idiag_apmax=0
  integer :: idiag_dvp12m=0, idiag_dtsweepp=0
!
  contains
!***********************************************************************
    subroutine register_particles_radius()
!
!  Set up indices for access to the fp and dfp arrays
!
!  22-aug-05/anders: coded
!
      if (lroot) call svn_id( &
          "$Id$")
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
    subroutine initialize_particles_radius(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  22-aug-05/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
!  Calculate the number density of bodies within a superparticle.
!
      if (npart_radii > 1 .and. &
          (.not. lcartesian_coords .or. &
          lparticles_nbody .or. &
          lparticles_number .or. &
          lparticles_spin)) then 
        call fatal_error('initialize_particles_radius: npart_radii > 1','')
      else
        mp_tilde=4/3.*pi*rhops*ap0(1)**3
        if (lroot) print*, 'initialize_particles_radius: '// &
            'mass per dust grain mp_tilde=', mp_tilde
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_radius
!***********************************************************************
    subroutine set_particle_radius(f,fp,npar_low,npar_high,init)
!
!  Set radius of new particles.
!
!  18-sep-09/nils: adapted from init_particles_radius
!
      use General, only: random_number_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer :: npar_low,npar_high
      logical, optional :: init
      logical :: initial
      real :: radius_fraction
!
      integer :: j,ind
!
      initial=.false.
      if (present(init)) then
        if (init) initial=.true.
      endif
!
      do j=1,ninit

        select case(initap(j))

        case('nothing')
          if (initial.and.lroot.and.j==1) print*, 'set_particles_radius: nothing'

        case('constant')
          if (initial.and.lroot) print*, 'set_particles_radius: constant radius'
          call random_number_wrapper(radius_fraction)
          ind=ceiling(npart_radii*radius_fraction)
          fp(npar_low:npar_high,iap)=ap0(ind)
        endselect

      enddo
!
      call keep_compiler_quiet(f)
!
    endsubroutine set_particle_radius
!***********************************************************************
    subroutine pencil_criteria_par_radius()
!
!  All pencils that the Particles_radius module depends on are specified here.
!
!  21-nov-06/anders: coded
!
      lpenc_requested(i_uu)=.true.
      lpenc_requested(i_rho)=.true.
      lpenc_requested(i_cc)=.true.
!
    endsubroutine pencil_criteria_par_radius
!***********************************************************************
    subroutine dap_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle radius.
!
!  22-aug-05/anders: coded
!
      use Diagnostics
      use Particles_number
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (nx) :: dt1_sweepup
      real :: deltavp, np_tilde
      integer :: k, ix0
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
      if (lheader) print*,'dap_dt_pencil: Calculate dap_dt'
!
!  Increase in particle radius due to sweep-up of small grains in the gas.
!
      if (lsweepup_par .and. t>=tstart_sweepup_par) then
!
        if (lfirst.and.ldt) dt1_sweepup=0.0
!
        if (npar_imn(imn)/=0) then
          do k=k1_imn(imn),k2_imn(imn)
            ix0=ineargrid(k,1)
!  No interpolation needed here.
!  Relative speed.
            deltavp=sqrt( &
                (fp(k,ivpx)-p%uu(ix0-nghost,1))**2 + &
                (fp(k,ivpy)-p%uu(ix0-nghost,2))**2 + &
                (fp(k,ivpz)-p%uu(ix0-nghost,3))**2 )
            if (deltavp12_floor/=0.0) &
                deltavp=sqrt(deltavp**2+deltavp12_floor**2)
!  Allow boulders to sweep up small grains if relative velocity not too high.
            if (deltavp<=vthresh_sweepup .or. vthresh_sweepup<0.0) then
              if (.not. lpscalar) then
                call fatal_error('dap_dt', &
                    'must have passive scalar module for sweep-up')
              else
!  Radius increase due to sweep-up
                dfp(k,iap) = dfp(k,iap) + &
                    0.25*deltavp*p%cc(ix0-nghost)*p%rho(ix0-nghost)/rhops
!
!  Deplete gas of small grains.
!
                call get_nptilde(fp,k,np_tilde)
                if (lpscalar_nolog) then
                  df(ix0,m,n,icc) = df(ix0,m,n,icc) - &
                      np_tilde*pi*fp(k,iap)**2*deltavp*p%cc(ix0-nghost)
                else
                  df(ix0,m,n,ilncc) = df(ix0,m,n,ilncc) - &
                      np_tilde*pi*fp(k,iap)**2*deltavp
                endif
!  Time-step contribution
                if (lfirst.and.ldt) then
                  dt1_sweepup(ix0-nghost) = dt1_sweepup(ix0-nghost) + &
                      np_tilde*pi*fp(k,iap)**2*deltavp
                endif
!
              endif
            endif
!
            if (ldiagnos) then
              if (idiag_dvp12m/=0) call sum_par_name((/deltavp/),idiag_dvp12m)
            endif
          enddo
        endif
!  Time-step contribution
          if (lfirst.and.ldt) then
            dt1_sweepup=dt1_sweepup/cdtps
            dt1_max=max(dt1_max,dt1_sweepup)
            if (ldiagnos.and.idiag_dtsweepp/=0) &
                call max_mn_name(dt1_sweepup,idiag_dtsweepp,l_dt=.true.)
          endif
!
      endif
!
      lfirstcall=.false.
!
      call keep_compiler_quiet(f)
!
    endsubroutine dap_dt_pencil
!***********************************************************************
    subroutine dap_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle radius.
!
!  21-nov-06/anders: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
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
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dap_dt
!***********************************************************************
    subroutine read_particles_rad_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
      integer :: i
!
      if (present(iostat)) then
        read(unit,NML=particles_radius_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_radius_init_pars,ERR=99)
      endif
!
! Find how many different particle radii we are using
! This must be done because not all parts of the code are adapted to 
! work with more than one particle radius.
!
      do i=1,ninit
        if (ap0(i) .ne. 0) then
          npart_radii=npart_radii+1
        endif
      enddo
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
      use Diagnostics
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
        idiag_dvp12m=0; idiag_dtsweepp=0
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
        call parse_name(iname,cname(iname),cform(iname),'dtsweepp',idiag_dtsweepp)
      enddo
!
    endsubroutine rprint_particles_radius
!***********************************************************************

endmodule Particles_radius
