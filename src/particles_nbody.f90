! $Id: particles_nbody.f90,v 1.10 2006-08-27 20:13:58 wlyra Exp $
!
!  This module takes care of everything related to particle self-gravity.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_nbody=.true.
!
! MPVAR CONTRIBUTION 6
!
!***************************************************************
module Particles_nbody

  use Cdata
  use Messages
  use Particles_cdata
  use Particles_sub

  implicit none

  include 'particles_nbody.h'
  
  real :: dummy
  !real, pointer :: tstart_nbody
  logical :: lnbody_particles=.true.
  real, dimension(nspar) :: xsp0=0.0, ysp0=0.0, zsp0=0.0
  real, dimension(nspar) :: vspx0=0.0, vspy0=0.0, vspz0=0.0
  real :: delta_vsp0=1.0
  real :: gc=0.
  character (len=labellen) :: initxxsp='origin', initvvsp='nothing'
  logical :: lcalc_orbit=.true.
  integer, dimension(nspar) :: ispar

  namelist /particles_nbody_init_pars/ &
       initxxsp, initvvsp, xsp0, ysp0, zsp0, vspx0, vspy0, vspz0, delta_vsp0, &
       bcspx, bcspy, bcspz

  namelist /particles_nbody_run_pars/ &
       bcspx, bcspy, bcspz, dsnap_par_minor, linterp_reality_check, &
       lcalc_orbit

  integer :: idiag_xpm=0, idiag_ypm=0, idiag_zpm=0
  integer :: idiag_xp2m=0, idiag_yp2m=0, idiag_zp2m=0
  integer :: idiag_vpxm=0, idiag_vpym=0, idiag_vpzm=0
  integer :: idiag_vpx2m=0, idiag_vpy2m=0, idiag_vpz2m=0
  integer :: idiag_xstar=0,idiag_ystar=0,idiag_zstar=0
  integer :: idiag_vxstar=0,idiag_vystar=0,idiag_vzstar=0
  integer :: idiag_vxplanet=0,idiag_vyplanet=0,idiag_vzplanet=0
  integer :: idiag_xplanet=0,idiag_yplanet=0,idiag_zplanet=0


  contains

!***********************************************************************
    subroutine register_particles_nbody()
!
!  Set up indices for access to the f, fsp and dfsp arrays.
!
!  27-aug-06/wlad: adapted
!
      use Messages, only: fatal_error, cvs_id
!
      logical, save :: first=.true.
!
      if (.not. first) &
          call fatal_error('register_particles_nbody: called twice','')
      first = .false.
!
      if (lroot) call cvs_id( &
           "$Id: particles_nbody.f90,v 1.10 2006-08-27 20:13:58 wlyra Exp $")
!
!  Check that we aren't registering too many auxiliary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux= ', maux
        call fatal_error('register_shock','naux > maux')
      endif
!     
    endsubroutine register_particles_nbody
!***********************************************************************
    subroutine initialize_particles_nbody(lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
! 
!  27-aug-06/wlad: adapted
!
      !use SharedVariables
!
      integer :: k,ierr
      logical :: lstarting
!
      if (NO_WARN) print*, lstarting
!
      !call get_shared_variable('tstart_nbody',tstart_nbody,ierr)
!
! Tag the sink particles -- The first particles on root are the sink ones
! ipar was initialized on initialize_particles
!
      if (lroot) then
         do k=1,nspar
            ispar(k)=ipar(k)
         enddo
      endif
!
      if (ierr/=0) then
        if (lroot) print*, 'initialize_particles_nbody: '// &
            'there was a problem when getting tstart_nbody!'
        call fatal_error('initialize_particles_nbody','')
      endif
!
    endsubroutine initialize_particles_nbody
!***********************************************************************
    subroutine pencil_criteria_par_nbody()
!   
!  All pencils that the Particles_nbody module depends on are specified here.
!         
!  02-jul-06/anders: adapted
!
    endsubroutine pencil_criteria_par_nbody
!***********************************************************************
    subroutine pencil_interdep_par_nbody(lpencil_in)
!   
!  Interdependency among pencils provided by the Particles_nbody module
!  is specified here.
!         
!  02-jul-06/anders: adapted
!
      logical, dimension(npencils) :: lpencil_in
!
      if (NO_WARN) print*, lpencil_in
!   
    endsubroutine pencil_interdep_par_nbody
!***********************************************************************
    subroutine calc_pencils_par_nbody(f,p)
!   
!  Calculate particle pencils.
!
!  02-jul-06/anders: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      if (NO_WARN) print*, f, p
!   
    endsubroutine calc_pencils_par_nbody
!***********************************************************************
    subroutine init_particles_nbody(f,fp)
!
!  Initial positions and velocities of sink particles.
!
!  17-nov-05/anders+wlad: adapted
!
      use General, only: random_number_wrapper
      use Mpicomm, only: stop_it
      use Sub
      use Gravity, only: g0
      use Planet,  only: get_ramped_mass
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: uup
      real :: gp,gs,gp_init
      real :: r, p
      integer :: k
!
      intent (in) :: f
      intent (out) :: fp
!
!  Initial particle position.
!
      select case(initxxsp)

      case ('origin')
        if (lroot) print*, 'init_particles_nbody: All sink particles at origin'
        fp(1:nspar,ixp:izp)=0.

      case('constant')
         if (lroot) &
              print*, 'init_particles_nbody: All sink particles at x,y,z=', xsp0, ysp0, zsp0
         fp(1:nspar,ixp)=xsp0
         fp(1:nspar,iyp)=ysp0
         fp(1:nspar,izp)=zsp0

      case ('fixed-cm')
!
! center of mass fixed on center of grid
!
         call get_ramped_mass(gp,gs,g0)
         print*,'initparticles_nbody . ispar(1) and ispar(2) equal to positions=',ispar(1),ispar(2)
         fp(ispar(2),ixp)=gp
         fp(ispar(2),iyp)=0.
         fp(ispar(2),izp)=0.
!
         fp(ispar(1),ixp) = -gs
         fp(ispar(1),iyp) = 0.
         fp(ispar(1),izp) = 0.
!
         if (lroot) &
              print*, 'init_particles_nbody: All sink particles at x=', fp(ispar(1),ixp), fp(ispar(2),ixp)

      case default
        if (lroot) print*, 'init_particles: No such such value for initxxsp: ', &
            trim(initxxsp)
        call stop_it("")

      endselect
!
!  Particles are not allowed to be present in non-existing dimensions.
!  This would give huge problems with interpolation later.
!
      if (nxgrid==1) fp(1:nspar,ixp)=x(nghost+1)
      if (nygrid==1) fp(1:nspar,iyp)=y(nghost+1)
      if (nzgrid==1) fp(1:nspar,izp)=z(nghost+1)
!
!  Redistribute particles among processors (now that positions are determined).
!
      call boundconds_sink_particles(fp)
!
!  Initial particle velocity.
!
      select case(initvvsp)

      case ('nothing')
        if (lroot) print*, 'init_particles: No particle velocity set'
      case ('zero')
        if (lroot) print*, 'init_particles: Zero particle velocity'
        fp(1:nspar,ivpx:ivpz)=0.

      case ('constant')
         if (lroot) print*, 'init_particles: Constant particle velocity'
         if (lroot) print*, 'init_particles: vspx0, vspy0, vspz0=', vspx0, vspy0, vspz0
         fp(1:nspar,ivpx)=vspx0
         fp(1:nspar,ivpy)=vspy0
         fp(1:nspar,ivpz)=vspz0

      case ('fixed-cm')
         call get_ramped_mass(gp,gs,g0)
!
         fp(ispar(2),ivpx)=0.
         fp(ispar(2),ivpy)=gp
         fp(ispar(2),ivpz)=0.
!
         fp(ispar(1),ivpx) = 0.
         fp(ispar(1),ivpy) = -gs
         fp(ispar(1),ivpz) = 0.
!
         if (lroot) &
              print*, 'init_particles: All particles at vy=', fp(ispar(1),ivpy), fp(ispar(2),ivpy)

       case default
        if (lroot) print*, 'init_particles: No such such value for initvvsp: ', &
            trim(initvvsp)
        call stop_it("")

      endselect
!
    endsubroutine init_particles_nbody
!***********************************************************************
    subroutine dxxp_dt_nbody_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle position (called from main pencil loop).
!
!  25-apr-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      if (NO_WARN) print*, f, df, fp, dfp, p, ineargrid
!
    endsubroutine dxxp_dt_nbody_pencil
!**********************************************************************
    subroutine dvvp_dt_nbody_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Add self-gravity to particle equation of motion.
!
!  14-jun-06/anders: coded
!
      use Messages, only: fatal_error
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (in) :: f, df, p, fp, dfp, ineargrid
!
      if (NO_WARN) print*, f, df, fp, dfp, p, ineargrid
!
    endsubroutine dvvp_dt_nbody_pencil
!***********************************************************************
    subroutine dxxp_dt_nbody(f,df,fp,dfp,ineargrid)
!
!  Evolution of sink particles position.
!
!  17-nov-05/anders+wlad: adapted
!
      use General, only: random_number_wrapper, random_seed_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real :: ran_xp, ran_yp, ran_zp
      integer, dimension (mseed) :: iseed_org
      integer :: k
      logical :: lheader, lfirstcall=.true.
!
      intent (in) :: f, fp, ineargrid
      intent (inout) :: df, dfp
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
!
!  Identify module and boundary conditions.
!
      if (lheader) print*,'dxxp_dt: Calculate dxxp_dt'
      if (lheader) then
        print*, 'dxxp_dt_nbody: Particles boundary condition bcpx=', bcspx
        print*, 'dxxp_dt_nbody: Particles boundary condition bcpy=', bcspy
        print*, 'dxxp_dt_nbody: Particles boundary condition bcpz=', bcspz
      endif
!
      if (lheader) print*, 'dxxp_dt: Set rate of change of particle '// &
          'position equal to particle velocity.'
!
!  The rate of change of a particle's position is the particle's velocity.
!
      if (nxgrid/=1) &
          dfp(1:nspar,ixp) = dfp(1:nspar,ixp) + fp(1:nspar,ivpx)
      if (nygrid/=1) &
          dfp(1:nspar,iyp) = dfp(1:nspar,iyp) + fp(1:nspar,ivpy)
      if (nzgrid/=1) &
          dfp(1:nspar,izp) = dfp(1:nspar,izp) + fp(1:nspar,ivpz)
!
!  With shear there is an extra term due to the background shear flow.
!
      if (lshear.and.nygrid/=0) dfp(1:nspar,iyp) = &
          dfp(1:nspar,iyp) - qshear*Omega*fp(1:nspar,ixp)
!
!  With masses and disk gravity, the torques of the disk must be
!  counter-balanced to keep the center of mass fixed in space
!
      !if (lmigrate) call reset_center_of_mass(fp,dfp,gp,gs)
!
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine dxxp_dt_nbody
!***********************************************************************
    subroutine dvvp_dt_nbody(f,df,fp,dfp,ineargrid)
!
!  Evolution of sink particles velocities
!
!  27-aug-06/wlad: coded
!
      use Cdata
      use Mpicomm, only: stop_it
      use Sub
      use Gravity, only: g0
      use Planet,only: get_ramped_mass,lmigrate
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
      real, dimension (npar_loc) :: xstar,ystar,zstar,vxstar,vystar,vzstar
      real, dimension (npar_loc) :: xplanet,yplanet,zplanet
      real, dimension (npar_loc) :: vxplanet,vyplanet,vzplanet
      real :: Omega2,rsep,gs,gp,gs_acc,gp_acc,r1sep
      real :: ax,ay,az,axs,ays,azs,vx,vy,vz,vxs,vys,vzs
      real :: axr,ayr,azr,vxr,vyr,vzr
      integer :: i, k
      logical :: lheader, lfirstcall=.true.
!
      intent (in) :: f, fp, ineargrid
      intent (inout) :: df, dfp
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
!
!  Identify module.
!
      if (lheader) print*,'dvvp_dt_nbody: Calculate dvvp_dt_nbody'
!
!  Add Coriolis force from rotating coordinate frame.
!
      if (Omega/=0.) then
         if (lheader) print*,'dvvp_dt: Add Coriolis force; Omega=', Omega
         Omega2=2*Omega
         dfp(1:npar_loc,ivpx) = dfp(1:npar_loc,ivpx) + Omega2*fp(1:npar_loc,ivpy)
         dfp(1:npar_loc,ivpy) = dfp(1:npar_loc,ivpy) - Omega2*fp(1:npar_loc,ivpx)
!
!  With shear there is an extra term due to the background shear flow.
!
         if (lshear) dfp(1:npar_loc,ivpy) = &
              dfp(1:npar_loc,ivpy) + qshear*Omega*fp(1:npar_loc,ivpx)
      endif
!
!  More readable variables names
!
      ax = fp(1,ixp)  ; axs = fp(2,ixp)
      ay = fp(1,iyp)  ; ays = fp(2,iyp)
      az = fp(1,izp)  ; azs = fp(2,izp)
!
!
      vx = fp(1,ivpx) ; vxs = fp(2,ivpx)
      vy = fp(1,ivpy) ; vys = fp(2,ivpy)
      vz = fp(1,ivpz) ; vzs = fp(2,ivpz)
!
!  Relative positions and velocities
!
      axr = ax - axs ; ayr = ay - ays ; azr = az - azs
      vxr = vx - vxs ; vyr = vy - vys ; vzr = vz - vzs
!                                                                               
      if (lroot) then
!                                                                               
         rsep = sqrt(axr**2 + ayr**2 + azr**2)
         r1sep = 1./rsep
!
!  Planet's gravity on star - must use ramp up as well
!
         call get_ramped_mass(gp,gs,g0)
!
         gs_acc = -gs*r1sep**2
         gp_acc = gs_acc * gp/gs
!
         dfp(2,ivpx) = dfp(2,ivpx) + gp_acc*r1sep*(axs-ax)
         dfp(2,ivpy) = dfp(2,ivpy) + gp_acc*r1sep*(ays-ay)
         dfp(2,ivpz) = dfp(2,ivpz) + gp_acc*r1sep*(azs-az)
!
!  Star's gravity on planet
!
         dfp(1,ivpx) = dfp(1,ivpx) + gs_acc*r1sep*(ax-axs)
         dfp(1,ivpy) = dfp(1,ivpy) + gs_acc*r1sep*(ay-ays)
         dfp(1,ivpz) = dfp(1,ivpz) + gs_acc*r1sep*(az-azs)
!
         call reset_center_of_mass(fp,dfp,gp,gs)
!
      endif
!
! At the end of the time-step, check the position of the center of
! mass and update the positions of all particles to keep it at rest
! It shouldn't move much, though. Just that the disk gravity inserts
! some small torques that make it move.
!
      if (lcalc_orbit) then
         xstar(1:npar_loc) = axs  ; vxstar(1:npar_loc) = vxs
         ystar(1:npar_loc) = ays  ; vystar(1:npar_loc) = vys
         zstar(1:npar_loc) = azs  ; vzstar(1:npar_loc) = vzs
!                                                                               
         xplanet(1:npar_loc) = ax ; vxplanet(1:npar_loc) = vx
         yplanet(1:npar_loc) = ay ; vyplanet(1:npar_loc) = vy
         zplanet(1:npar_loc) = az ; vzplanet(1:npar_loc) = vz
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
        if (idiag_xstar/=0)    call sum_par_name(xstar(1:npar_loc),idiag_xstar)
        if (idiag_ystar/=0)    call sum_par_name(ystar(1:npar_loc),idiag_ystar)
        if (idiag_zstar/=0)    call sum_par_name(zstar(1:npar_loc),idiag_zstar)
        if (idiag_xplanet/=0)  call sum_par_name(xplanet(1:npar_loc),idiag_xplanet)
        if (idiag_yplanet/=0)  call sum_par_name(yplanet(1:npar_loc),idiag_yplanet)
        if (idiag_zplanet/=0)  call sum_par_name(zplanet(1:npar_loc),idiag_zplanet)
        if (idiag_vxstar/=0)   call sum_par_name(vxstar(1:npar_loc),idiag_vxstar)
        if (idiag_vystar/=0)   call sum_par_name(vystar(1:npar_loc),idiag_vystar)
        if (idiag_vzstar/=0)   call sum_par_name(vzstar(1:npar_loc),idiag_vzstar)
        if (idiag_vxplanet/=0) call sum_par_name(vxplanet(1:npar_loc),idiag_vxplanet)
        if (idiag_vyplanet/=0) call sum_par_name(vyplanet(1:npar_loc),idiag_vyplanet)
        if (idiag_vzplanet/=0) call sum_par_name(vzplanet(1:npar_loc),idiag_vzplanet)
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
    endsubroutine dvvp_dt_nbody
!***********************************************************************
    subroutine read_particles_nbody_init_pars(unit,iostat)
!    
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_nbody_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_nbody_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_nbody_init_pars
!***********************************************************************
    subroutine write_particles_nbody_init_pars(unit)
!    
      integer, intent (in) :: unit
!
      write(unit,NML=particles_nbody_init_pars)
!
    endsubroutine write_particles_nbody_init_pars
!***********************************************************************
    subroutine read_particles_nbody_run_pars(unit,iostat)
!    
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_nbody_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_nbody_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_nbody_run_pars
!***********************************************************************
    subroutine write_particles_nbody_run_pars(unit)
!    
      integer, intent (in) :: unit
!
      write(unit,NML=particles_nbody_run_pars)
!
    endsubroutine write_particles_nbody_run_pars
!***********************************************************************
    subroutine rprint_particles_nbody(lreset,lwrite)
!
!  Read and register print parameters relevant for sink particles.
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
        idiag_xstar=0; idiag_ystar=0; idiag_zstar=0
        idiag_vxstar=0; idiag_vystar=0; idiag_vzstar=0
        idiag_xplanet=0; idiag_yplanet=0; idiag_zplanet=0
        idiag_vxplanet=0; idiag_vyplanet=0; idiag_vzplanet=0
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
        call parse_name(iname,cname(iname),cform(iname),'xstar',idiag_xstar)
        call parse_name(iname,cname(iname),cform(iname),'ystar',idiag_ystar)
        call parse_name(iname,cname(iname),cform(iname),'zstar',idiag_zstar)
        call parse_name(iname,cname(iname),cform(iname),'xplanet',idiag_xplanet)
        call parse_name(iname,cname(iname),cform(iname),'yplanet',idiag_yplanet)
        call parse_name(iname,cname(iname),cform(iname),'zplanet',idiag_zplanet)
        call parse_name(iname,cname(iname),cform(iname),'vxplanet',idiag_vxplanet)
        call parse_name(iname,cname(iname),cform(iname),'vyplanet',idiag_vyplanet)
        call parse_name(iname,cname(iname),cform(iname),'vzplanet',idiag_vzplanet)
        call parse_name(iname,cname(iname),cform(iname),'vxstar',idiag_vxstar)
        call parse_name(iname,cname(iname),cform(iname),'vystar',idiag_vystar)
        call parse_name(iname,cname(iname),cform(iname),'vzstar',idiag_vzstar)
      enddo
!
    endsubroutine rprint_particles_nbody
!***********************************************************************
    subroutine boundconds_sink_particles(fp,dfp)
!
!  Global boundary conditions for particles.
!
!  27-aug-06/wlad: adapted
!
      use Messages, only: fatal_error_local
      use Mpicomm
!
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension (mpar_loc,mpvar), optional :: dfp
      integer :: k
!
      intent (inout) :: fp,dfp
!
!  Boundary condition in the x-direction.
!
      if (nxgrid/=1) then
         if (bcspx=='p') then
            do k=1,nspar
!  xp < x0
               if (fp(k,ixp)< xyz0(1)) then
                  fp(k,ixp)=fp(k,ixp)+Lxyz(1)
!
!  Particle position must never need more than one addition of Lx to get back
!  in the box. Often a NaN or Inf in the particle position will show up as a
!  problem here.
!
                  if (fp(k,ixp)< xyz0(1)) then
                     print*, 'boundconds_sink_particles: ERROR - sink particle ', &
                          ispar(k), ' was further than Lx outside the simulation box!'
                     print*, 'This must never happen.'
                     print*, 'iproc, ispar, xxp=', iproc, ispar(k), fp(k,ixp:izp)
                     call fatal_error_local('boundconds_sink_particles','')
                  endif
               endif
!  xp > x1
            if (fp(k,ixp)>=xyz1(1)) then
              fp(k,ixp)=fp(k,ixp)-Lxyz(1)
              if (fp(k,ixp)>=xyz1(1)) then
                print*, 'boundconds_sink_particles: ERROR - sink particle ', &
                     ispar(k), ' was further than Lx outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ispar, xxp=', iproc, ispar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
            endif
          enddo
        else
          print*, 'boundconds_particles: No such boundary condition bcpx=', bcspx
          call stop_it('boundconds_particles')
        endif
      endif
!
!  Boundary condition in the y-direction.
!
      if (nygrid/=1) then
        if (bcspy=='p') then
!  yp < y0
          do k=1,nspar
            if (fp(k,iyp)< xyz0(2)) then
              fp(k,iyp)=fp(k,iyp)+Lxyz(2)
              if (fp(k,iyp)< xyz0(2)) then
                print*, 'boundconds_particles: ERROR - particle ', ispar(k), &
                    ' was further than Ly outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ispar, xxp=', iproc, ispar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
            endif
!  yp > y1
            if (fp(k,iyp)>=xyz1(2)) then
              fp(k,iyp)=fp(k,iyp)-Lxyz(2)
              if (fp(k,iyp)>=xyz1(2)) then
                print*, 'boundconds_particles: ERROR - particle ', ispar(k), &
                    ' was further than Ly outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ispar, xxp=', iproc, ispar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
            endif
          enddo
        else
          print*, 'boundconds_particles: No such boundary condition bcpy=', bcspy
          call stop_it('boundconds_particles')
        endif
      endif
!
!  Boundary condition in the z-direction.
!
      if (nzgrid/=1) then
        if (bcspz=='p') then
          do k=1,nspar
!  zp < z0
            if (fp(k,izp)< xyz0(3)) then
              fp(k,izp)=fp(k,izp)+Lxyz(3)
              if (fp(k,izp)< xyz0(3)) then
                print*, 'boundconds_particles: ERROR - particle ', ispar(k), &
                    ' was further than Lz outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ipar, xxp=', iproc, ispar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
            endif
!  zp > z1
            if (fp(k,izp)>=xyz1(3)) then
              fp(k,izp)=fp(k,izp)-Lxyz(3)
              if (fp(k,izp)>=xyz1(3)) then
                print*, 'boundconds_particles: ERROR - particle ', ispar(k), &
                    ' was further than Lz outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ispar, xxp=', iproc, ispar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
            endif
          enddo
        else
          print*, 'boundconds_particles: No such boundary condition bcpz=', bcspz
          call stop_it('boundconds_particles')
        endif
      endif
!
!  Share all sink particles among processors (internal boundary conditions).
!
      call share_allparticles_procs(fp)
!
    endsubroutine boundconds_sink_particles
!***********************************************************************
    subroutine share_allparticles_procs(fp)
!
!  For N-body runs (few sink particles), keep particles at root
!  processor and inform other processors of positions and velocities.
!
!  27-aug-06/wlad: coded
!
      use Mpicomm
!
      real, dimension (mpar_loc,mpvar) :: fp
!
      integer :: k,nspar
!
      do k=1,nspar
        call mpibcast_real(fp(k,:),mpvar)
      enddo
!
    endsubroutine share_allparticles_procs
!***********************************************************************
    subroutine reset_center_of_mass(fp,dfp,gp,gs)
!
!  If the center of mass was accelerated, reset its position
!  to the center of the grig
!
!  27-aug-06/wlad: coded
!
      real, dimension (mpar_loc,mpvar) :: fp,dfp
      real :: gs,gp,invtotmass,vx_cm,vy_cm,vz_cm
!
      invtotmass = 1./(gp+gs)
!
      vx_cm = invtotmass * (gp*fp(1,ivpx) + gs*fp(2,ivpx))
      vy_cm = invtotmass * (gp*fp(1,ivpy) + gs*fp(2,ivpy))
      vz_cm = invtotmass * (gp*fp(1,ivpz) + gs*fp(2,ivpz))
!
      dfp(1,ixp) = dfp(1,ixp) - vx_cm
      dfp(2,ixp) = dfp(2,ixp) - vx_cm
!
      dfp(1,iyp) = dfp(1,iyp) - vy_cm
      dfp(2,iyp) = dfp(2,iyp) - vy_cm
!
      dfp(1,izp) = dfp(1,izp) - vz_cm
      dfp(2,izp) = dfp(2,izp) - vz_cm
!
     endsubroutine reset_center_of_mass
!***********************************************************************
    subroutine get_particles_interdistances(fp,rp_mn,rpcyl_mn)
!
!  Should be pencils of dimension (nx,nspar)
!  because they will be used by planet and gravity
!
!  18-jul-06/wlad: coded
!
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension (nx,nspar) :: rp_mn,rpcyl_mn
      integer :: i
!
      intent(out) :: rp_mn,rpcyl_mn
!
! more readable variable names
!
      do i=1,nspar
!
! Spherical and cylindrical distances
!
         rp_mn(:,i)    = &
              sqrt((x(l1:l2)-fp(i,ixp))**2 + (y(m)-fp(i,iyp))**2 &
              + (z(n)-fp(i,izp))**2) + tini
         rp_mn(:,i)    = &
              sqrt((x(l1:l2)-fp(i,ixp))**2 + (y(m)-fp(i,iyp))**2  &
              + (z(n)-fp(i,izp))**2) + tini
!
         rpcyl_mn(:,i) = &
              sqrt((x(l1:l2)-fp(i,ixp))**2 + (y(m)-fp(i,iyp))**2) + tini
!
      enddo
!
    endsubroutine get_particles_interdistances
!***********************************************************************
endmodule Particles_nbody
