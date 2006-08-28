! $Id: particles_nbody.f90,v 1.11 2006-08-28 20:35:03 wlyra Exp $
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
           "$Id: particles_nbody.f90,v 1.11 2006-08-28 20:35:03 wlyra Exp $")
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
! Tag the sink particles  - equal to the first ipars
!
      if (lroot) &
           print*,'particles_nbody: Initialize ispar array'
!
      do k=1,nspar
         ispar(k)=ipar(k)
         print*,k,ispar(k)
      enddo
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
    subroutine init_particles_nbody(f,fp,fsp)
!
!  Initial positions and velocities of sink particles.
!  Overwrite the position asserted by the dust module
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
      real, dimension (nspar,mpvar) :: fsp
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
      call boundconds_particles(fp,npar_loc,ipar)
      call share_sinkparticles(fp,fsp)
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
! Broadcast the velocities as well
!
      call share_sinkparticles(fp,fsp)
!
    endsubroutine init_particles_nbody
!***********************************************************************
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
    subroutine dvvp_dt_nbody(f,df,fp,dfp,fsp,dfsp,ineargrid)
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
      real, dimension (nspar,mpvar) :: fsp,dfsp
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
      intent (in) ::     f,  fp,  fsp, ineargrid
      intent (inout) :: df, dfp, dfsp
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
      !if (Omega/=0.) then
      !   if (lheader) print*,'dvvp_dt: Add Coriolis force; Omega=', Omega
      !   Omega2=2*Omega
      !   dfp(1:npar_loc,ivpx) = dfp(1:npar_loc,ivpx) + Omega2*fp(1:npar_loc,ivpy)
      !   dfp(1:npar_loc,ivpy) = dfp(1:npar_loc,ivpy) - Omega2*fp(1:npar_loc,ivpx)
!
!  With shear there is an extra term due to the background shear flow.
!
      !   if (lshear) dfp(1:npar_loc,ivpy) = &
      !        dfp(1:npar_loc,ivpy) + qshear*Omega*fp(1:npar_loc,ivpx)
      !endif
!
!  More readable variables names
!
      ax = fsp(1,ixp)  ; axs = fsp(2,ixp)
      ay = fsp(1,iyp)  ; ays = fsp(2,iyp)
      az = fsp(1,izp)  ; azs = fsp(2,izp)
!
      vx = fsp(1,ivpx) ; vxs = fsp(2,ivpx)
      vy = fsp(1,ivpy) ; vys = fsp(2,ivpy)
      vz = fsp(1,ivpz) ; vzs = fsp(2,ivpz)
!
!  Relative positions and velocities
!
      axr = ax - axs ; ayr = ay - ays ; azr = az - azs
      vxr = vx - vxs ; vyr = vy - vys ; vzr = vz - vzs
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
      dfsp(2,ivpx) = dfsp(2,ivpx) + gp_acc*r1sep*(axs-ax)
      dfsp(2,ivpy) = dfsp(2,ivpy) + gp_acc*r1sep*(ays-ay)
      dfsp(2,ivpz) = dfsp(2,ivpz) + gp_acc*r1sep*(azs-az)
!
!  Star's gravity on planet
!
      dfsp(1,ivpx) = dfsp(2,ivpx) + gs_acc*r1sep*(ax-axs)
      dfsp(1,ivpy) = dfsp(2,ivpy) + gs_acc*r1sep*(ay-ays)
      dfsp(1,ivpz) = dfsp(2,ivpz) + gs_acc*r1sep*(az-azs)
!
! Check the position of the center of mass of the sink 
! particles, and reset it to the center of the grid
!
      call reset_center_of_mass(fsp,dfsp,gp,gs)
!
!  Put it back on the dfp array
!
      do k=1,nspar
         dfp(ispar(k),:) = dfp(ispar(k),:) + dfsp(k,:)
      enddo
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
    subroutine reset_center_of_mass(fsp,dfsp,gp,gs)
!
!  If the center of mass was accelerated, reset its position
!  to the center of the grig
!
!  27-aug-06/wlad: coded
!
      real, dimension (nspar,mpvar) :: fsp,dfsp
      real :: gs,gp,invtotmass,vx_cm,vy_cm,vz_cm
!
      intent (in)    ::  fsp,gs,gp
      intent (inout) :: dfsp
!
      invtotmass = 1./(gp+gs)
!
      vx_cm = invtotmass * (gp*fsp(1,ivpx) + gs*fsp(2,ivpx))
      vy_cm = invtotmass * (gp*fsp(1,ivpy) + gs*fsp(2,ivpy))
      vz_cm = invtotmass * (gp*fsp(1,ivpz) + gs*fsp(2,ivpz))
!
      dfsp(1,ixp) = dfsp(1,ixp) - vx_cm
      dfsp(2,ixp) = dfsp(2,ixp) - vx_cm
!                        
      dfsp(1,iyp) = dfsp(1,iyp) - vy_cm
      dfsp(2,iyp) = dfsp(2,iyp) - vy_cm
!                        
      dfsp(1,izp) = dfsp(1,izp) - vz_cm
      dfsp(2,izp) = dfsp(2,izp) - vz_cm
!
     endsubroutine reset_center_of_mass
!***********************************************************************
    subroutine get_particles_interdistances(fsp,rp_mn,rpcyl_mn)
!
!  Should be pencils of dimension (nx,nspar)
!  because they will be used by planet and gravity
!
!  18-jul-06/wlad: coded
!
      real, dimension (nspar,mpvar) :: fsp
      real, dimension (nx,nspar) :: rp_mn,rpcyl_mn
      integer :: k
!
      intent(out) :: rp_mn,rpcyl_mn
!
! more readable variable names
!
      do k=1,nspar
!
! Spherical and cylindrical distances
!
         rp_mn(:,ispar(k))    = &
              sqrt((x(l1:l2)-fsp(k,ixp))**2 + (y(m)-fsp(k,iyp))**2  &
              + (z(n)-fsp(k,izp))**2) + tini
!
         rpcyl_mn(:,ispar(k)) = &
              sqrt((x(l1:l2)-fsp(k,ixp))**2 + (y(m)-fsp(k,iyp))**2) + tini
!
      enddo
!
    endsubroutine get_particles_interdistances
!***********************************************************************
endmodule Particles_nbody
