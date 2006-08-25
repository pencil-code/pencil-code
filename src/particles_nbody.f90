! $Id: particles_nbody.f90,v 1.8 2006-08-25 14:50:52 wlyra Exp $
!
!  This module takes care of everything related to sink particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 6
! CPARAM logical, parameter :: lparticles=.true.
! CPARAM logical, parameter :: lparticles_nbody=.true.
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
  integer :: idiag_xstar=0,idiag_ystar=0,idiag_zstar=0
  integer :: idiag_vxstar=0,idiag_vystar=0,idiag_vzstar=0
  integer :: idiag_vxplanet=0,idiag_vyplanet=0,idiag_vzplanet=0
  integer :: idiag_xplanet=0,idiag_yplanet=0,idiag_zplanet=0

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
           "$Id: particles_nbody.f90,v 1.8 2006-08-25 14:50:52 wlyra Exp $")
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
      select case(initxxp)

      case ('origin')
        if (lroot) print*, 'init_particles: All particles at origin'
        fp(1:npar_loc,ixp:izp)=0.
        
      case('constant') 
         if (lroot) &
              print*, 'init_particles: All particles at x,y,z=', xp0, yp0, zp0
         fp(1:npar_loc,ixp)=xp0
         fp(1:npar_loc,iyp)=yp0
         fp(1:npar_loc,izp)=zp0

      case ('fixed-cm')
!
! center of mass fixed on center of grid
!
         call get_ramped_mass(gp,gs,g0)
         fp(2,ixp)=gp 
         fp(2,iyp)=0.
         fp(2,izp)=0.
!
         fp(1,ixp) = -gs 
         fp(1,iyp) = 0.
         fp(1,izp) = 0. 
!
         if (lroot) &
              print*, 'init_particles: All particles at x=', fp(1,ixp), fp(2,ixp)

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
      if (nxgrid==1) fp(1:npar_loc,ixp)=x(nghost+1)
      if (nygrid==1) fp(1:npar_loc,iyp)=y(nghost+1)
      if (nzgrid==1) fp(1:npar_loc,izp)=z(nghost+1)
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

      case ('fixed-cm')
         call get_ramped_mass(gp,gs,g0)
!
         fp(2,ivpx)=0.
         fp(2,ivpy)=gp
         fp(2,ivpz)=0.
!
         fp(1,ivpx) = 0.
         fp(1,ivpy) = -gs
         fp(1,ivpz) = 0.
!
         if (lroot) &
              print*, 'init_particles: All particles at vy=', fp(1,ivpy), fp(2,ivpy)

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
          call interpolate_quadratic(f,iux,iuz,fp(k,ixp:izp),uup,ineargrid(k,:))
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
    subroutine pencil_criteria_particles()
! 
!  All pencils that the Particles module depends on are specified here.
! 
!  20-04-06/anders: coded
! 
    endsubroutine pencil_criteria_particles
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
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      if (NO_WARN) print*, f, p
!
    endsubroutine calc_pencils_particles
!***********************************************************************
    subroutine dxxp_dt_pencil(f,df,fp,dfp,p,ineargrid)
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
    endsubroutine dxxp_dt_pencil
!***********************************************************************
    subroutine dvvp_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of sink particles velocity due to gas gravity
!  These particle are way to big and massive to be influenced
!  by gas drag
!
!  24-aug-06/wlad: coded
!
      use Mpicomm
      use Gravity, only: g0,r0_pot
      use Planet, only: Gvalue,b_pot,lmigrate
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      real, dimension(nx) :: re,grav_gas,re1
      real :: sumx_loc,sumy_loc,sumz_loc,sumx,sumy,sumz
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
      integer :: k
      real :: r_smooth,dv
      logical :: lheader,lfirstcall=.true.    
!
! Header information
!
      lheader=lfirstcall .and. lroot
!
! Exclude the gas away from the accretion radius of the particle
! For the planet, the accretion radius is the Roche Lobe...
! And in the case of the star, what would it be?
! Shouldn't it be general? Like.... a pre-set accretion radius?
!
!
      if (lmigrate) then
         if (lheader) print*,'adding gas gravity to particles!'
!
         do k=1,npar
!
            if (k==2) then
               r_smooth = r0_pot     !star
            else if (k==1) then
               r_smooth = b_pot      !planet
            else                                                           
               call stop_it('gravity_gas - more than 2 planet particles?')
            endif
!
            re = sqrt((x(l1:l2) - fp(k,ixp))**2 &
                 +    (y(  m  ) - fp(k,iyp))**2 &
                 +    (z(  n  ) - fp(k,izp))**2)
!
            dv = dx*dy
            if (nzgrid/=1) dv = dv*dz
!                                                                      
! Gas gravity
!
            grav_gas = Gvalue*p%rho*dv*re*(re**2 + r_smooth**2)**(-1.5)
            re1=1./re
!
            where (re.le.r_smooth)
               grav_gas=0.
            endwhere
!
            sumx_loc(1) = sum(grav_gas * (x(l1:l2) - fp(k,ixp))*re1)
            sumy_loc(1) = sum(grav_gas * (y(  m  ) - fp(k,iyp))*re1)        
            sumz_loc(1) = sum(grav_gas * (z(  n  ) - fp(k,izp))*re1)        
!                                                                      
            call mpireduce_sum_scl(sumx_loc,sumx)
            call mpireduce_sum_scl(sumy_loc,sumy)
            call mpireduce_sum_scl(sumz_loc,sumz)
!
            if (lroot) then
               dfp(k,ivpx) = dfp(k,ivpx) + sumx
               dfp(k,ivpy) = dfp(k,ivpy) + sumy
               dfp(k,ivpz) = dfp(k,ivpz) + sumz
            endif
!
         enddo
      endif
!
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine dvvp_dt_pencil
!***********************************************************************
    subroutine dxxp_dt(f,df,fp,dfp,ineargrid)
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
          dfp(1:npar_loc,ixp) = dfp(1:npar_loc,ixp) + fp(1:npar_loc,ivpx)
      if (nygrid/=1) &
          dfp(1:npar_loc,iyp) = dfp(1:npar_loc,iyp) + fp(1:npar_loc,ivpy)
      if (nzgrid/=1) &
          dfp(1:npar_loc,izp) = dfp(1:npar_loc,izp) + fp(1:npar_loc,ivpz)
!
!  With shear there is an extra term due to the background shear flow.
!
      if (lshear.and.nygrid/=0) dfp(1:npar_loc,iyp) = &
          dfp(1:npar_loc,iyp) - qshear*Omega*fp(1:npar_loc,ixp)
!
!  With masses and disk gravity, the torques of the disk must be 
!  counter-balanced to keep the center of mass fixed in space
!
      !if (lmigrate) call reset_center_of_mass(fp,dfp,gp,gs)
!
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine dxxp_dt
!***********************************************************************
    subroutine dvvp_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of sink particles velocities
!
!  17-nov-05/anders+wlad: coded
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
      real, dimension (mpar_loc) :: xstar,ystar,zstar,vxstar,vystar,vzstar
      real, dimension (mpar_loc) :: xplanet,yplanet,zplanet
      real, dimension (mpar_loc) :: vxplanet,vyplanet,vzplanet
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
      if (lheader) print*,'dvvp_dt: Calculate dvvp_dt'
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
    endsubroutine dvvp_dt
!***********************************************************************
    subroutine reset_center_of_mass(fp,dfp,gp,gs)
!
! If it was accelerated, reset its position!
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
! 18-jul-06/wlad: coded
!
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension (nx,mpar_loc) :: rp_mn,rpcyl_mn
      integer :: i
!      
      intent(out) :: rp_mn,rpcyl_mn
!
! more readable variable names
!
      do i=1,mpar_loc
!
! Spherical and cylindrical distances
! 
         rp_mn(:,i)    = &
              sqrt((x(l1:l2)-fp(i,ixp))**2 + (y(m)-fp(i,iyp))**2  + (z(n)-fp(i,izp))**2) + tini
!
         rpcyl_mn(:,i) = &
              sqrt((x(l1:l2)-fp(i,ixp))**2 + (y(m)-fp(i,iyp))**2) + tini
!
      enddo
!
    endsubroutine get_particles_interdistances
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
      real, dimension (mx,my,mz,mfarray) :: f
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
    endsubroutine rprint_particles
!***********************************************************************
endmodule Particles
