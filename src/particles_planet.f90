! $Id: particles_planet.f90,v 1.23 2006-03-28 21:52:27 wlyra Exp $
!
!  This module takes care of everything related to planet particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 6
! CPARAM logical, parameter :: lparticles=.true.
! CPARAM logical, parameter :: lparticles_planet=.true.
!
!***************************************************************
module Particles

  use Cdata
  use Particles_cdata
  use Particles_sub
  use Messages

  implicit none

  include 'particles.h'

  real, dimension(npar) :: xp0=0.0, yp0=0.0, zp0=0.0
  real, dimension(npar) :: vpx0=0.0, vpy0=0.0, vpz0=0.0
  real :: delta_vp0=1.0, cdtp=0.2
  real :: gc=0.
  character (len=labellen) :: initxxp='origin', initvvp='nothing'
  logical :: lcalc_orbit=.true.

  namelist /particles_init_pars/ &
      initxxp, initvvp, xp0, yp0, zp0, vpx0, vpy0, vpz0, delta_vp0, &
      bcpx, bcpy, bcpz

  namelist /particles_run_pars/ &
      bcpx, bcpy, bcpz, dsnap_par_minor, cdtp, linterp_reality_check, &
      lcalc_orbit

  integer :: idiag_xpm=0, idiag_ypm=0, idiag_zpm=0
  integer :: idiag_xp2m=0, idiag_yp2m=0, idiag_zp2m=0
  integer :: idiag_vpxm=0, idiag_vpym=0, idiag_vpzm=0
  integer :: idiag_vpx2m=0, idiag_vpy2m=0, idiag_vpz2m=0
  integer :: idiag_smap=0,idiag_eccp=0
  integer :: idiag_xstar=0,idiag_ystar=0
  integer :: idiag_vxstar=0,idiag_vystar=0
  integer :: idiag_vxplanet=0,idiag_vyplanet=0
  integer :: idiag_xplanet=0,idiag_yplanet=0

  contains

!***********************************************************************
    subroutine register_particles()
!
!  Set up indices for access to the fp and dfp arrays
!
!  17-nov-05/anders+wlad: adapted
!
      use Mpicomm, only: stop_it
!
      integer :: k
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_particles: called twice')
      first = .false.
!
      if (lroot) call cvs_id( &
           "$Id: particles_planet.f90,v 1.23 2006-03-28 21:52:27 wlyra Exp $")
!
!  Indices for particle position.
!
      ixp=npvar+1
      iyp=npvar+2
      izp=npvar+3
!
!  Indices for particle velocity.
!
      ivpx=npvar+4
      ivpy=npvar+5
      ivpz=npvar+6
!
!  Increase npvar accordingly.
!
      npvar=npvar+6
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call stop_it('register_particles: npvar > mpvar')
      endif
!
!  Check that we aren't registering too many auxilary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
            call stop_it('register_particles: naux > maux')
      endif
!
!  Set npar_loc=npar for non-parallel implementation of few particles.
!
      if (lroot) then
        npar_loc=npar
      else
        npar_loc=0
      endif
      do k=1,npar
        ipar(k)=k
      enddo
!
    endsubroutine register_particles
!***********************************************************************
    subroutine initialize_particles(lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  17-nov-05/anders+wlad: adapted
!
      logical :: lstarting
!
    endsubroutine initialize_particles
!***********************************************************************
    subroutine init_particles(f,fp,ineargrid)
!
!  Initial positions and velocities of planet particles.
!
!  17-nov-05/anders+wlad: adapted
!
      use Boundcond
      use General, only: random_number_wrapper
      use Mpicomm, only: stop_it
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: uup
      real :: r, p
      integer :: k
!
      intent (in) :: f
      intent (out) :: fp
!
!  Initial particle position.
!
      select case(initxxp)

      case ('origin')
        if (lroot) print*, 'init_particles: All particles at origin'
        fp(1:npar,ixp:izp)=0.

      case ('constant')
        if (lroot) &
            print*, 'init_particles: All particles at x,y,z=', xp0, yp0, zp0
        fp(1:npar_loc,ixp)=xp0
        fp(1:npar_loc,iyp)=yp0
        fp(1:npar_loc,izp)=zp0

      case ('random')
        if (lroot) print*, 'init_particles: Random particle positions'
        do k=1,npar_loc
          if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
          if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
          if (nzgrid/=1) call random_number_wrapper(fp(k,izp))
        enddo
        if (nxgrid/=1) fp(1:npar_loc,ixp)=xyz0(1)+fp(1:npar_loc,ixp)*Lxyz(1)
        if (nygrid/=1) fp(1:npar_loc,iyp)=xyz0(2)+fp(1:npar_loc,iyp)*Lxyz(2)
        if (nzgrid/=1) fp(1:npar_loc,izp)=xyz0(3)+fp(1:npar_loc,izp)*Lxyz(3)

      case default
        if (lroot) print*, 'init_particles: No such such value for initxxp: ', &
            trim(initxxp)
        call stop_it("")

      endselect
!      
!  Particles are not allowed to be present in non-existing dimensions.
!  This would give huge problems with interpolation later.
!
      if (nxgrid==1) fp(1:npar,ixp)=x(nghost+1)
      if (nygrid==1) fp(1:npar,iyp)=y(nghost+1)
      if (nzgrid==1) fp(1:npar,izp)=z(nghost+1)
!
!  Redistribute particles among processors (now that positions are determined).
!
      call boundconds_particles(fp,npar_loc,ipar)
!
!  Initial particle velocity.
!
      select case(initvvp)

      case ('nothing')
        if (lroot) print*, 'init_particles: No particle velocity set'
      case ('zero')
        if (lroot) print*, 'init_particles: Zero particle velocity'
        fp(1:npar_loc,ivpx:ivpz)=0.

      case ('constant')
        if (lroot) print*, 'init_particles: Constant particle velocity'
        if (lroot) print*, 'init_particles: vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
        fp(1:npar_loc,ivpx)=vpx0
        fp(1:npar_loc,ivpy)=vpy0
        fp(1:npar_loc,ivpz)=vpz0

      case ('random')
        if (lroot) print*, 'init_particles: Random particle velocities; '// &
            'delta_vp0=', delta_vp0
        do k=1,npar_loc
          call random_number_wrapper(fp(k,ivpx))
          call random_number_wrapper(fp(k,ivpy))
          call random_number_wrapper(fp(k,ivpz))
        enddo
        fp(1:npar_loc,ivpx) = -delta_vp0 + fp(1:npar_loc,ivpx)*2*delta_vp0
        fp(1:npar_loc,ivpy) = -delta_vp0 + fp(1:npar_loc,ivpy)*2*delta_vp0
        fp(1:npar_loc,ivpz) = -delta_vp0 + fp(1:npar_loc,ivpz)*2*delta_vp0

      case ('follow-gas')
        if (lroot) &
            print*, 'init_particles: Particle velocity equal to gas velocity'
        do k=1,npar_loc
          call interpolate_3d_1st(f,iux,iuz,fp(k,ixp:izp),uup,ineargrid(k,:))
          fp(k,ivpx:ivpz) = uup
        enddo

      case default
        if (lroot) print*, 'init_particles: No such such value for initvvp: ', &
            trim(initvvp)
        call stop_it("")

      endselect
!
    endsubroutine init_particles
!***********************************************************************
    subroutine pencil_interdep_particles(lpencil_in)
!   
!  Interdependency among pencils provided by the Particles module
!  is specified here.
!         
!  16-feb-06/anders: dummy
!
      logical, dimension(npencils) :: lpencil_in
!
      if (NO_WARN) print*, lpencil_in
!
    endsubroutine pencil_interdep_particles
!***********************************************************************
    subroutine calc_pencils_particles(f,p)
!   
!  Calculate particle pencils.
!
!  16-feb-06/anders: dummy
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p
!
      if (NO_WARN) print*, f, p
!
    endsubroutine calc_pencils_particles
!***********************************************************************
    subroutine dxxp_dt(f,fp,dfp,ineargrid)
!
!  Evolution of planet particle position.
!
!  17-nov-05/anders+wlad: adapted
!
      use General, only: random_number_wrapper, random_seed_wrapper
!      
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real :: ran_xp, ran_yp, ran_zp
      integer, dimension (mseed) :: iseed_org
      integer :: k
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
      if (lheader) print*,'dxxp_dt: Calculate dxxp_dt'
      if (lheader) then
        print*, 'dxxp_dt: Particles boundary condition bcpx=', bcpx
        print*, 'dxxp_dt: Particles boundary condition bcpy=', bcpy
        print*, 'dxxp_dt: Particles boundary condition bcpz=', bcpz
      endif
!
      if (lheader) print*, 'dxxp_dt: Set rate of change of particle '// &
          'position equal to particle velocity.'
!
!  The rate of change of a particle's position is the particle's velocity.
!
      if (nxgrid/=1) &
          dfp(1:npar,ixp) = dfp(1:npar,ixp) + fp(1:npar,ivpx)
      if (nygrid/=1) &
          dfp(1:npar,iyp) = dfp(1:npar,iyp) + fp(1:npar,ivpy)
      if (nzgrid/=1) &
          dfp(1:npar,izp) = dfp(1:npar,izp) + fp(1:npar,ivpz)
!
!  With shear there is an extra term due to the background shear flow.
!
      if (lshear.and.nygrid/=0) dfp(1:npar,iyp) = &
          dfp(1:npar,iyp) - qshear*Omega*fp(1:npar,ixp)
!
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine dxxp_dt
!***********************************************************************
    subroutine dvvp_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of planet velocity and star velocity
!
!  17-nov-05/anders+wlad: coded
!  01-feb-06/wlad : disk backreaction moved to planet.f90
!
      use Cdata
      use EquationOfState, only: cs20, gamma
      use Mpicomm, only: stop_it
      use Sub
      use Gravity,only: g0
      use Planet,only: get_ramped_mass
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      real, dimension (mpar_loc) :: sma_planet,ecc_planet
      real, dimension (mpar_loc) :: xstar,ystar,vxstar,vystar
      real, dimension (mpar_loc) :: xplanet,yplanet,vxplanet,vyplanet
      real :: sma,ecc,rv,rrdot,w2,mdot,m2dot,extra_accx,extra_accy
      real :: vx,vy,axr,ayr,vxr,vyr,vxs,vys,mu,tcut,gp_d  
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (nx) :: re, grav_gas
      real :: Omega2,rsep,gs,gp
      integer :: i, k, ix0, iy0, iz0,mg,ng
      logical :: lheader, lfirstcall=.true.
      real :: ax,ay,axs,ays
!
      intent (in) :: fp
      intent (out) :: dfp
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
!
!  Identify module.
!
      if (lheader) print*,'dvvp_dt: Calculate dvvp_dt'
!
!  Add Coriolis force from rotating coordinate frame.
!
      if (Omega/=0.) then
         if (lheader) print*,'dvvp_dt: Add Coriolis force; Omega=', Omega
         Omega2=2*Omega
         dfp(1:npar,ivpx) = dfp(1:npar,ivpx) + Omega2*fp(1:npar,ivpy)
         dfp(1:npar,ivpy) = dfp(1:npar,ivpy) - Omega2*fp(1:npar,ivpx)
!
!  With shear there is an extra term due to the background shear flow.
!          
         if (lshear) dfp(1:npar,ivpy) = &
              dfp(1:npar,ivpy) + qshear*Omega*fp(1:npar,ivpx)
      endif
!
!  More readable variables names
!
      ax = fp(1,ixp)  ; axs = fp(2,ixp)  
      ay = fp(1,iyp)  ; ays = fp(2,iyp) 
!
      vx = fp(1,ivpx) ; vxs = fp(2,ivpx)  
      vy = fp(1,ivpy) ; vys = fp(2,ivpy) 
!
!  Relative positions and velocities
!
      axr = ax - axs ; ayr = ay - ays
      vxr = vx - vxs ; vyr = vy - vys 
!
      if (lroot) then
!
         rsep = sqrt(axr**2 + ayr**2)
!      
!  Planet's gravity on star - must use ramp up as well
!        
         call get_ramped_mass(gp,gs,g0,mdot,m2dot)
!
         extra_accx = m2dot*cos(t) - 2*mdot*sin(t)
         extra_accy = m2dot*sin(t) + 2*mdot*cos(t)
!
         dfp(2,ivpx) = dfp(2,ivpx) - gp/rsep**3 * (axs-ax) + extra_accx 
         dfp(2,ivpy) = dfp(2,ivpy) - gp/rsep**3 * (ays-ay) + extra_accy
!
!  Star's gravity on planet
!
         dfp(1,ivpx) = dfp(1,ivpx) - gs/rsep**3 * (ax-axs) + extra_accx
         dfp(1,ivpy) = dfp(1,ivpy) - gs/rsep**3 * (ay-ays) + extra_accy
!
      endif
!
      if (lcalc_orbit) then

!  Some celestial mechanics can't do any harm
!
!  r2    = x**2  +  y**2 +  z**2
!  rrdot = x*vx  +  y*vy +  z*vz  
!  w2    = vx**2 + vy**2 + vz**2     
!
!  Semimajor axis - star at rest  
! 
! E = 0.5*mu*v**2 - GM*mu/r = -GM*mu/2a
! (E energy; mu reduced mass; a semimajor axis)
!
!    1/a = 2/r - w2/GM      
!  
!  eccentricity
!      
         rv    = sqrt(axr**2 + ayr**2)
         rrdot = axr*vxr + ayr*vyr
         w2    = vxr**2 + vyr**2
!
         sma = (2./rv - w2/g0)**(-1.)
!
         ecc = sqrt(rrdot**2/sma + (1 - rv/sma)**2)
!
         sma_planet(1:npar_loc) = sma
         ecc_planet(1:npar_loc) = ecc
!
         xstar(1:npar_loc) = axs  ; vxstar(1:npar_loc) = vxs
         ystar(1:npar_loc) = ays  ; vystar(1:npar_loc) = vys
!
         xplanet(1:npar_loc) = ax ; vxplanet(1:npar_loc) = vx
         yplanet(1:npar_loc) = ay ; vyplanet(1:npar_loc) = vy
!
      endif
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_xpm/=0)  call sum_par_name(fp(1:npar_loc,ixp),idiag_xpm)
        if (idiag_ypm/=0)  call sum_par_name(fp(1:npar_loc,iyp),idiag_ypm)
        if (idiag_zpm/=0)  call sum_par_name(fp(1:npar_loc,izp),idiag_zpm)
        if (idiag_xp2m/=0) call sum_par_name(fp(1:npar_loc,ixp)**2,idiag_xp2m)
        if (idiag_yp2m/=0) call sum_par_name(fp(1:npar_loc,iyp)**2,idiag_yp2m)
        if (idiag_zp2m/=0) call sum_par_name(fp(1:npar_loc,izp)**2,idiag_zp2m)
        if (idiag_vpxm/=0) call sum_par_name(fp(1:npar_loc,ivpx),idiag_vpxm)
        if (idiag_vpym/=0) call sum_par_name(fp(1:npar_loc,ivpy),idiag_vpym)
        if (idiag_vpzm/=0) call sum_par_name(fp(1:npar_loc,ivpz),idiag_vpzm)
!        
        if (idiag_smap/=0)     call sum_par_name(sma_planet(1:npar_loc),idiag_smap)
        if (idiag_eccp/=0)     call sum_par_name(ecc_planet(1:npar_loc),idiag_eccp)     
        if (idiag_xstar/=0)    call sum_par_name(xstar(1:npar_loc),idiag_xstar)
        if (idiag_ystar/=0)    call sum_par_name(ystar(1:npar_loc),idiag_ystar)
        if (idiag_xplanet/=0)  call sum_par_name(xplanet(1:npar_loc),idiag_xplanet)
        if (idiag_yplanet/=0)  call sum_par_name(yplanet(1:npar_loc),idiag_yplanet)
        if (idiag_vxstar/=0)   call sum_par_name(vxstar(1:npar_loc),idiag_vxstar)
        if (idiag_vystar/=0)   call sum_par_name(vystar(1:npar_loc),idiag_vystar)
        if (idiag_vxplanet/=0) call sum_par_name(vxplanet(1:npar_loc),idiag_vxplanet)
        if (idiag_vyplanet/=0) call sum_par_name(vyplanet(1:npar_loc),idiag_vyplanet)
!
        if (idiag_vpx2m/=0) &
            call sum_par_name(fp(1:npar_loc,ivpx)**2,idiag_vpx2m)
        if (idiag_vpy2m/=0) &
            call sum_par_name(fp(1:npar_loc,ivpy)**2,idiag_vpy2m)
        if (idiag_vpz2m/=0) &
            call sum_par_name(fp(1:npar_loc,ivpz)**2,idiag_vpz2m)
      endif
!
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine dvvp_dt
!***********************************************************************
    subroutine read_particles_init_pars(unit,iostat)
!    
!  17-nov-05/anders+wlad: adapted
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_init_pars
!***********************************************************************
    subroutine write_particles_init_pars(unit)
!    
!  17-nov-05/anders+wlad: adapted
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_init_pars)
!
    endsubroutine write_particles_init_pars
!***********************************************************************
    subroutine read_particles_run_pars(unit,iostat)
!    
!  17-nov-05/anders+wlad: adapted
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_run_pars
!***********************************************************************
    subroutine write_particles_run_pars(unit)
!    
!  17-nov-05/anders+wlad: adapted
!    
      integer, intent (in) :: unit
!
      write(unit,NML=particles_run_pars)
!
    endsubroutine write_particles_run_pars
!***********************************************************************
    subroutine powersnap_particles(f)
!
!  Calculate power spectra of dust particle variables.
!
!  01-feb-06/anders: dummy
!     
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if (NO_WARN) print*, f
!
    endsubroutine powersnap_particles
!***********************************************************************
    subroutine rprint_particles(lreset,lwrite)
!   
!  Read and register print parameters relevant for particles.
!
!  17-nov-05/anders+wlad: adapted
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
      
      if (lwr) then
        write(3,*) 'ixp=', ixp
        write(3,*) 'iyp=', iyp
        write(3,*) 'izp=', izp
        write(3,*) 'ivpx=', ivpx
        write(3,*) 'ivpy=', ivpy
        write(3,*) 'ivpz=', ivpz
      endif
!
!  Reset everything in case of reset
!
      if (lreset) then
        idiag_xpm=0; idiag_ypm=0; idiag_zpm=0
        idiag_xp2m=0; idiag_yp2m=0; idiag_zp2m=0
        idiag_vpxm=0; idiag_vpym=0; idiag_vpzm=0
        idiag_vpx2m=0; idiag_vpy2m=0; idiag_vpz2m=0
        idiag_smap=0; idiag_eccp=0 
        idiag_xstar=0; idiag_ystar=0 
        idiag_vxstar=0; idiag_vystar=0
        idiag_xplanet=0; idiag_yplanet=0 
        idiag_vxplanet=0; idiag_vyplanet=0
      endif
!
!  Run through all possible names that may be listed in print.in
!
      if (lroot .and. ip<14) print*,'rprint_particles: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'xpm',idiag_xpm)
        call parse_name(iname,cname(iname),cform(iname),'ypm',idiag_ypm)
        call parse_name(iname,cname(iname),cform(iname),'zpm',idiag_zpm)
        call parse_name(iname,cname(iname),cform(iname),'xp2m',idiag_xp2m)
        call parse_name(iname,cname(iname),cform(iname),'yp2m',idiag_yp2m)
        call parse_name(iname,cname(iname),cform(iname),'zp2m',idiag_zp2m)
        call parse_name(iname,cname(iname),cform(iname),'vpxm',idiag_vpxm)
        call parse_name(iname,cname(iname),cform(iname),'vpym',idiag_vpym)
        call parse_name(iname,cname(iname),cform(iname),'vpzm',idiag_vpzm)
        call parse_name(iname,cname(iname),cform(iname),'vpx2m',idiag_vpx2m)
        call parse_name(iname,cname(iname),cform(iname),'vpy2m',idiag_vpy2m)
        call parse_name(iname,cname(iname),cform(iname),'vpz2m',idiag_vpz2m)
        call parse_name(iname,cname(iname),cform(iname),'smap',idiag_smap)
        call parse_name(iname,cname(iname),cform(iname),'xstar',idiag_xstar)
        call parse_name(iname,cname(iname),cform(iname),'eccp',idiag_eccp)
        call parse_name(iname,cname(iname),cform(iname),'ystar',idiag_ystar)
        call parse_name(iname,cname(iname),cform(iname),'xplanet',idiag_xplanet)
        call parse_name(iname,cname(iname),cform(iname),'yplanet',idiag_yplanet)
        call parse_name(iname,cname(iname),cform(iname),'vxplanet',idiag_vxplanet)
        call parse_name(iname,cname(iname),cform(iname),'vyplanet',idiag_vyplanet)
        call parse_name(iname,cname(iname),cform(iname),'vxstar',idiag_vxstar)
        call parse_name(iname,cname(iname),cform(iname),'vystar',idiag_vystar)
      enddo
!
    endsubroutine rprint_particles
!***********************************************************************
endmodule Particles
