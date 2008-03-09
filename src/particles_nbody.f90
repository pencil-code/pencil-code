! $Id: particles_nbody.f90,v 1.64 2008-03-09 17:11:25 wlyra Exp $
!
!  This module takes care of everything related to sink particles.
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
!!!! 08 mar 08, supposed to be the new one, apart from a break in eos
!    i still don't understand
  use Cdata
  use Messages
  use Particles_cdata
  use Particles_sub

  implicit none

  include 'particles_nbody.h'

  real, dimension(nspar,mpvar) :: fsp
  real, dimension(nspar) :: xsp0=0.0, ysp0=0.0, zsp0=0.0
  real, dimension(nspar) :: vspx0=0.0, vspy0=0.0, vspz0=0.0
  real, dimension(nspar) :: pmass=1.,r_smooth,pmass1,final_ramped_mass=1.
  real :: delta_vsp0=1.0,totmass,totmass1
  character (len=labellen) :: initxxsp='random', initvvsp='nothing'
  logical :: lcalc_orbit=.true.,lbackreaction=.false.,lnorm=.true.
  logical :: lreset_cm=.false.,lnogravz_star=.false.,lexclude_frozen=.false.
  logical :: lnoselfgrav_star=.true.
  logical, dimension(nspar) :: lcylindrical_gravity_nbody=.false.
  logical, dimension(nspar) :: lfollow_particle=.false.
  real :: GNewton=impossible,prhs_cte
  integer :: ramp_orbits=5
  logical :: lramp=.false.
  logical :: linterpolate_gravity=.false.,linterpolate_linear=.true.
  logical :: linterpolate_quadratic_spline=.false.

  integer :: iglobal_ggp=0

  namelist /particles_nbody_init_pars/ &
       initxxsp, initvvsp, xsp0, ysp0, zsp0, vspx0, vspy0, vspz0, delta_vsp0, &
       pmass, r_smooth, lcylindrical_gravity_nbody, &
       lexclude_frozen, GNewton, bcspx, bcspy, bcspz, &
       ramp_orbits,lramp,final_ramped_mass,prhs_cte,linterpolate_gravity,&
       linterpolate_quadratic_spline

  namelist /particles_nbody_run_pars/ &
       dsnap_par_minor, linterp_reality_check, lcalc_orbit, lreset_cm, &
       lnogravz_star,lfollow_particle, lbackreaction, lexclude_frozen, &
       GNewton, bcspx, bcspy, bcspz,prhs_cte,lnoselfgrav_star,&
       linterpolate_quadratic_spline

  integer, dimension(nspar,3) :: idiag_xxspar=0,idiag_vvspar=0
  integer, dimension(nspar)   :: idiag_torqint=0,idiag_torqext=0
  integer                     :: idiag_totenergy=0,idiag_totangmom=0

  logical :: ldust=.false.

  contains

!***********************************************************************
    subroutine register_particles_nbody()
!
!  Set up indices for access to the f and fsp
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
           "$Id: particles_nbody.f90,v 1.64 2008-03-09 17:11:25 wlyra Exp $")
!
!  No need to solve the N-body equations for non-N-body problems.
!
      if (nspar < 2) then
        if (lroot) write(0,*) 'nspar = ', nspar
        call fatal_error('register_particles_nbody','the number of massive'//&
             ' particles is less than 2. There is no need to use the'//&
             ' N-body code. Consider setting gravity as a global variable'//&
             ' using gravity_r.f90 instead.')
      endif
!
!  Check that we aren't registering too many auxiliary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux= ', maux
        call fatal_error('register_particles_nbody','naux > maux')
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
      use Mpicomm,only:stop_it
      use FArrayManager
!
      logical :: lstarting
!
! G_Newton. Overwrite the one set by start.in if set again here, 
! because I might want units in which both G and GM are 1. 
!
      if (GNewton == impossible) then
        GNewton=G_Newton
      endif
!
! inverse mass
!
      if (lramp) then
        pmass(1:nspar-1) = epsi
        pmass(nspar)=1-epsi*(nspar-1)
      endif
!
      pmass1=1./pmass
!
! inverse total mass
!
      totmass=sum(pmass)
      totmass1=1./totmass
!
! check for consistency between the cylindrical gravities switches
! from cdata and the one from nbody
!
      if (((lcylindrical_gravity).and.&
        (.not.lcylindrical_gravity_nbody(nspar))).or.&
             (.not.lcylindrical_gravity).and.&
             (lcylindrical_gravity_nbody(nspar))) then 
        call stop_it("initialize_particles_nbody: inconsitency "//&
             "between lcylindrical_gravity from cdata and the "//&
             "one from nbody")
      endif
!
      if (rsmooth.ne.r_smooth(nspar)) then 
        call stop_it("initialize_particles_nbody: inconsitency "//&
             "between rsmooth from cdata and the "//&
             "one from nbody")
      endif
!
! Check for consistency between the poisson and integral 
! calculations of selfgravity
!
      if (lbackreaction.and.lselfgravity) then 
        print*,'self-gravity is already taken into '//&
             'account. lbackreaction is a way of '//&
             'integrating the self-gravity only in few '//&
             'specific points of the grid'
        call fatal_error("initialize_particles_nbody","")
      endif
!
!  The presence of dust particles needs to be known
!
      if (npar > nspar) ldust=.true.
      if (linterpolate_gravity.and.(.not.ldust)) then
        if (lroot) print*,'interpolate gravity is just'//&
             ' for the dust component. No need for it if'//&
             ' you are not using dust particles'
        call fatal_error("initialize_particles_nbody","")
      endif
      if (linterpolate_gravity) then
         if (lroot) print*,'initializing global array for sink gravity'
         call farray_register_global('ggp',iglobal_ggp,vector=3)
      endif
!      
      if (linterpolate_quadratic_spline) &
           linterpolate_linear=.false.
!
      if ((.not.linterpolate_gravity).and.&
           linterpolate_quadratic_spline) then
        if (lroot) print*,'no need for linterpolate_quadratic_spline'//&
             'if linterpolate_gravity is false'
        call fatal_error("initialize_particles_nbody","")
      endif
!
    endsubroutine initialize_particles_nbody
!***********************************************************************
    subroutine pencil_criteria_par_nbody()
!
!  All pencils that the Particles_nbody module depends on are specified here.
!
!  22-sep-06/wlad: adapted
!
      lpenc_requested(i_rho)=.true.
      if (ldust) lpenc_requested(i_rhop)=.true.
      lpenc_requested(i_r_mn)=.true.
      lpenc_requested(i_rcyl_mn)=.true.
      
      
      if (idiag_totenergy/=0) then
        lpenc_diagnos(i_u2)=.true.
        lpenc_diagnos(i_rcyl_mn)=.true.
      endif
!
      if (any(idiag_torqext/=0) .or. any(idiag_torqint/=0)) then
        lpenc_diagnos(i_rcyl_mn)=.true.
      endif
!
    endsubroutine pencil_criteria_par_nbody
!***********************************************************************
    subroutine pencil_interdep_par_nbody(lpencil_in)
!
!  Interdependency among pencils provided by the Particles_nbody module
!  is specified here.
!
!  22-sep-06/wlad: adapted
!
      logical, dimension(npencils) :: lpencil_in
!
      if (NO_WARN) print*, lpencil_in
!
    endsubroutine pencil_interdep_par_nbody
!***********************************************************************
    subroutine calc_pencils_par_nbody(f,p)
!
!  Calculate sink particle pencils
!
!  22-sep-06/wlad: adapted
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
!  Overwrite the position asserted by the dust module
!
!  17-nov-05/anders+wlad: adapted
!
      use General, only: random_number_wrapper
      use Sub
      use Mpicomm
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      real, dimension(nspar) :: ang_vel,kep_vel,sma
      real, dimension(nspar,3) :: position,velocity
      real :: tmp,tmp2,rr,fac,parc
      integer :: k,ks

      intent (in) :: f
      intent (out) :: fp
!
!  Initialize particles' positions and velocities
!
      position(:,1) = xsp0 ; position(:,2) = ysp0 ; position(:,3) = zsp0
      velocity(:,1) = vspx0; velocity(:,2) = vspy0; velocity(:,3) = vspz0
!
      select case(initxxsp)

      case ('origin')
        if (lroot) then
          print*, 'init_particles_nbody: All sink particles at origin'
          fp(1:nspar,ixp:izp)=0.
        endif
!
      case('constant')
        if (lroot) then
          print*, 'init_particles_nbody: All sink particles at x,y,z=', xsp0, ysp0, zsp0
          fp(1:nspar,ixp)=xsp0
          fp(1:nspar,iyp)=ysp0
          fp(1:nspar,izp)=zsp0
        endif
!
      case ('random')
        if (lroot) print*, 'init_particles_nbody: Random particle positions'
        do ks=1,nspar
          if (nxgrid/=1) call random_number_wrapper(position(ks,ixp))
          if (nygrid/=1) call random_number_wrapper(position(ks,iyp))
          if (nzgrid/=1) call random_number_wrapper(position(ks,izp))
        enddo

        if (nxgrid/=1) &
             position(1:nspar,ixp)=xyz0_loc(1)+position(1:nspar,ixp)*Lxyz_loc(1)
        if (nygrid/=1) &
             position(1:nspar,iyp)=xyz0_loc(2)+position(1:nspar,iyp)*Lxyz_loc(2)
        if (nzgrid/=1) &
             position(1:nspar,izp)=xyz0_loc(3)+position(1:nspar,izp)*Lxyz_loc(3)
!
!  Loop through ipar to allocate the sink particles
!
        do k=1,npar_loc
          if (ipar(k) <= nspar) then
            print*,'initparticles_nbody. Slot for sink particle ',ipar(k),&
                 ' was at fp position ',k,' at processor ',iproc
!
            fp(k,ixp:izp)=position(ipar(k),1:3)
!
! Correct for non-existing dimensions (not really needed, I think)
!
            if (nxgrid==1) fp(k,ixp)=x(nghost+1)
            if (nygrid==1) fp(k,iyp)=y(nghost+1)
            if (nzgrid==1) fp(k,izp)=z(nghost+1)
!
            print*,'initparticles_nbody. Sink particle ',ipar(k),&
                 ' located at ixp=',fp(k,ixp)
          endif
        enddo
!
      case ('fixed-cm')
!
! Ok, I have the masses and the positions of all sinks except the last, 
! which will have a position determined to fix the center of mass on 
! the center of the grid
!
        if (any(ysp0/=0)) call stop_it("init_particles_nbody: not yet generalized"//&
             " for non-zero azimuthal initial position")
!
        if (any(zsp0/=0)) call stop_it("init_particles_nbody: nbody code not"//&
             " yet generalized to allow initial inclinations")
!
        if (lcylindrical_coords) then
          if (any(xsp0.lt.0)) &
               call stop_it("init_particles_nbody: in cylindrical coordinates"//&
               " all the radial positions must be positive")
        endif
!
        if (lroot) then
          print*,'fixed-cm: redefining the mass of the last sink particle'
          print*,'fixed-cm: it assumes that the sum of the mass'//&
               ' of the particles is always g0'
        endif
!
        tmp = 0.;parc=0
        do ks=1,nspar-1 
          sma(ks)=abs(position(ks,1))
          tmp=tmp+pmass(ks)
          parc = parc - sma(ks)*pmass(ks)
        enddo
        pmass(nspar)=1.- tmp;pmass1=1./pmass;totmass=1.;totmass1=1.
        parc = parc*totmass1
        if (tmp .ge. 1.) &
             call stop_it("particles_nbody,init_particles. The mass of one "//& 
             "(or more) of the particles is too big! The masses should "//&
             "never be bigger than g0. Please scale your assemble so that "//&
             "the combined mass of the (n-1) particles is less than that. "//&
             "The mass of the last particle in the pmass array will be "//&
             "reassigned to ensure that the total mass is g0")
!
        do ks=1,nspar-1 
          position(ks,1)=sign(1.,position(ks,1))* (sma(ks) + parc)
        enddo
!
! The last one (star) fixes the CM at Rcm=zero
!
        position(nspar,1)=parc
        if (lcylindrical_coords) then
          !put the star in positive coordinates, with pi for azimuth
          position(nspar,1)=abs(parc)
          position(nspar,2)=pi
        endif
!
        if (lroot) then 
          print*,'pmass =',pmass
          print*,'position (x)=',position(:,1)
          print*,'position (y)=',position(:,2)
          print*,'position (z)=',position(:,3)
        endif
!
! Loop through ipar to allocate the sink particles
!
        do k=1,npar_loc
          if (ipar(k) <= nspar) then
!
            print*,'initparticles_nbody. Slot for sink particle ',ipar(k),&
                 ' was at fp position ',k,' at processor ',iproc
!
            fp(k,ixp:izp)=position(ipar(k),1:3)
!
! Correct for non-existing dimensions (not really needed, I think)
!
            if (nxgrid==1) fp(k,ixp)=x(nghost+1)
            if (nygrid==1) fp(k,iyp)=y(nghost+1)
            if (nzgrid==1) fp(k,izp)=z(nghost+1)
!
            print*,'initparticles_nbody. Sink particle ',ipar(k),&
                 ' located at ixp=',fp(k,ixp)
          endif
        enddo
!
      case default
        if (lroot) print*,"init_particles: No such such value for"//&
             " initxxsp: ",trim(initxxsp)
        call stop_it("")
!
      endselect
!
!  Redistribute particles among processors (now that positions are determined).
!
      call boundconds_particles(fp,npar_loc,ipar)
!
!  Initial particle velocity.
!
      select case(initvvsp)

      case ('nothing')
        if (lroot) print*, 'init_particles: No particle velocity set'
      case ('zero')
        if (lroot) then
          print*, 'init_particles: Zero particle velocity'
          fp(1:nspar,ivpx:ivpz)=0.
        endif
!
      case ('constant')
         if (lroot) then
           print*, 'init_particles: Constant particle velocity'
           print*, 'init_particles: vspx0, vspy0, vspz0=', vspx0, vspy0, vspz0
           fp(1:nspar,ivpx)=vspx0
           fp(1:nspar,ivpy)=vspy0
           fp(1:nspar,ivpz)=vspz0
         endif
!
      case ('fixed-cm')
!
! Keplerian velocities for the planets
!
        parc=0.
        do ks=1,nspar-1 
          kep_vel(ks)=sqrt(1./sma(ks)) !circular velocity
          parc = parc - kep_vel(ks)*pmass(ks)
        enddo
        parc = parc*totmass
        do ks=1,nspar-1 
          if (lcartesian_coords) then
            velocity(ks,2) = sign(1.,position(ks,1))*(kep_vel(ks) + parc)
          elseif (lcylindrical_coords) then
            !positive for the planets
            velocity(ks,2) = abs(kep_vel(ks) + parc)
          endif
        enddo
!
! The last one (star) fixes the CM also with velocity zero
!
        velocity(nspar,2)=parc
        if (lcylindrical_coords) velocity(nspar,2)=-parc
!
! Loop through ipar to allocate the sink particles
!
         do k=1,npar_loc
           if (ipar(k)<=nspar) then
             print*,&
                  'initparticles_nbody. Slot for sink particle ',ipar(k),&
                  ' was at fp position ',k,&
                  ' at processor ',iproc
!
             fp(k,ivpx:ivpz) = velocity(ipar(k),1:3)
!
             print*,'initparticles_nbody. Sink particle ',ipar(k),&
                  ' has velocity y=',fp(k,ivpy)
!
           endif
         enddo
!
      case default
        if (lroot) print*, "init_particles: No such such value for"//&
             "initvvsp: ",trim(initvvsp)
        call stop_it("")
!
      endselect
!
! Make the particles known to all processors
!
      call boundconds_particles(fp,npar_loc,ipar)
      call share_sinkparticles(fp)
!
    endsubroutine init_particles_nbody
!***********************************************************************
    subroutine dvvp_dt_nbody_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Add gravity from the particles to the gas
!  and the backreaction from the gas onto the particles
!
!  Adding the gravity on the dust component via the grid
!  is less accurate, but much faster than looping through all 
!  npar_loc particle and add the gravity of all sinks to it.
!
!  07-sep-06/wlad: coded
!
      use Messages, only: fatal_error
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      real, dimension (nx,nspar) :: rp_mn,rpcyl_mn
      real, dimension (nx,3) :: ggt
      real, dimension (nx) :: pot_energy
      real, dimension (3) :: xxpar,accg
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
      integer :: ks,k
      logical :: lintegrate,lparticle_out
!
      intent (in) :: f, p, fp, ineargrid
      intent (inout) :: df,dfp
!
! Get the total gravity field. In the case of dust,
! it is already pre-calculated
!
      if (linterpolate_gravity) then
!
! Interpolate the gravity to the position of the particles
!
          if (npar_imn(imn)/=0) then
            do k=k1_imn(imn),k2_imn(imn) !loop throught in the pencil
              if (ipar(k).gt.nspar) then !for dust
                !interpolate the gravity
                if (linterpolate_linear) then
                  call interpolate_linear(f,iglobal_ggp,&
                       iglobal_ggp+2,fp(k,ixp:izp),accg,ineargrid(k,:),ipar(k))
                else if (linterpolate_quadratic_spline) then
                  !
                  ! WL: I am not sure if this interpolation
                  !     works for cylindrical coordinates, so
                  !     beware
                  !
                  call interpolate_quadratic_spline(f,iglobal_ggp,&
                       iglobal_ggp+2,fp(k,ixp:izp),accg,ineargrid(k,:),ipar(k))
                endif
                dfp(k,ivpx:ivpz)=dfp(k,ivpx:ivpz)+accg
              endif
            enddo
          endif
      endif
! 
! Add the acceleration to the gas
!
      if (lhydro) then
        if (linterpolate_gravity) then
          !already calculate, so no need to lose time
          !calculating again
          ggt=f(l1:l2,m,n,iglobal_ggp:iglobal_ggp+2)
        else
          call get_total_gravity(ggt)
        endif
        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + ggt 
      
!
! Backreaction of the gas+dust gravity onto the massive particles
! The integration is needed for these two cases:
!  
!   1. We are not solving the Poisson equation (lbackreaction)
!   2. We are, but a particle is out of the box (a star, for instance)
!      and therefore the potential cannot be interpolated.
!
        do ks=1,nspar
          lparticle_out=.false.
          if (lselfgravity) then
            if ((fsp(ks,ixp)< xyz0(1)).or.(fsp(ks,ixp) > xyz1(1)) .or. &
                (fsp(ks,iyp)< xyz0(2)).or.(fsp(ks,iyp) > xyz1(2)) .or. &
                (fsp(ks,izp)< xyz0(3)).or.(fsp(ks,izp) > xyz1(3))) then 
              !particle out of box
              lparticle_out=.true.
            endif
          endif
          lintegrate=lbackreaction.or.lparticle_out
!
! Sometimes making the star feel the selfgravity of the disk
! leads to numerical troubles as the star is too close to the 
! origin (in cylindrical coordinates).
!
          if ((ks==nspar).and.lnoselfgrav_star) &
               lintegrate=.false.
!
          if (lintegrate) then
!
! Get the acceleration particle ks suffers due to self-gravity
!
            xxpar = fsp(ks,ixp:izp)
            if (lcylindrical_gravity_nbody(ks)) then
              call integrate_selfgravity(p,rpcyl_mn(:,ks),&
                   xxpar,accg,r_smooth(ks))
            else
              call integrate_selfgravity(p,rp_mn(:,ks),&
                   xxpar,accg,r_smooth(ks))
            endif
            
!
! Add it to its dfp
!
            do k=1,npar_loc
              if (ipar(k)==ks) &
                   dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + accg(1:3)
            enddo
          endif
!
! Diagnostic
!
          if (ldiagnos) then
            if (idiag_totenergy/=0) pot_energy=0.
            call get_radial_distance(rp_mn(:,ks),rpcyl_mn(:,ks),&
                 E1_=fsp(ks,ixp),E2_=fsp(ks,iyp),E3_=fsp(ks,izp))
!
! Calculate torques for output, if needed
!
            if ((idiag_torqext(ks)/=0).or.(idiag_torqint(ks)/=0)) &
                 call calc_torque(p,rpcyl_mn(:,ks),ks)
!
! Total energy      
!
            if (idiag_totenergy/=0) then
              !potential energy
              pot_energy = pot_energy - &
                   GNewton*pmass(ks)*&
                   (rpcyl_mn(:,ks)**2+r_smooth(ks)**2)**(-0.5)
              if (ks==nspar) &
                   call sum_lim_mn_name(.5*p%rho*p%u2 + pot_energy,&
                   idiag_totenergy,p)
            endif
          endif
        enddo
      endif
!
    endsubroutine dvvp_dt_nbody_pencil
!***********************************************************************
    subroutine dxxp_dt_nbody(dfp)
!
!  If the center of mass of the sink particles was moved from the
!  center of the grid, reset it.
!
!  22-sep-06/wlad: coded
!
      real, dimension (mpar_loc,mpvar) :: dfp
!
      if (lreset_cm) &
           call reset_center_of_mass(dfp)
!
    endsubroutine dxxp_dt_nbody
!***********************************************************************
    subroutine dvvp_dt_nbody(f,df,fp,dfp,ineargrid)
!
!  Evolution of sink and dust particles velocities due to
!  particle-particle interaction only. 
!
!  Coriolis and shear already added in particles_dust
!
!  27-aug-06/wlad: coded
!
      use Cdata
      use Mpicomm, only: mpibcast_real,stop_it
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension(3) :: xxpar,accg
      integer :: k, ks, j, jpos, jvel
      logical :: lheader, lfirstcall=.true.
!
      intent (in) ::     f,  fp,  ineargrid
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
!  Evolve sink particle positions due to the gravity of other 
!  sink particles. The gravity of the sinks on the dust will be
!  added inside the pencil in dvvp_dt_nbody_pencil
!
      if (lramp) call get_ramped_mass
!
      do k=1,npar_loc
        if (linterpolate_gravity) then
          !only loop through sinks
          if (ipar(k).le.nspar) then
            xxpar=fp(k,ixp:izp)
            call loop_through_sinks(xxpar,ipar(k),accg)
            dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + accg
          endif
        else
          !for all particles
          xxpar=fp(k,ixp:izp)
          call loop_through_sinks(xxpar,ipar(k),accg)
          dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + accg
        endif
      enddo
!
!  Position and velocity diagnostics (per sink particle)
! 
     if (ldiagnos) then
        do ks=1,nspar
          if (lfollow_particle(ks)) then
            do j=1,3
              jpos=j+ixp-1 ; jvel=j+ivpx-1
              if (idiag_xxspar(ks,j)/=0) &
                   call point_par_name(fsp(ks,jpos),idiag_xxspar(ks,j))
              if (idiag_vvspar(ks,j)/=0) &
                   call point_par_name(fsp(ks,jvel),idiag_vvspar(ks,j))
            enddo
          endif
        enddo
      endif
!
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine dvvp_dt_nbody
!************************************************************
    subroutine loop_through_sinks(xxpar,kj,accg)
!
      real, dimension (3) :: xxpar,evr, accg
      real :: e1,e2,e3,e10,e20,e30
      real :: r2_ij,invr3_ij
      integer :: k, ks, kj
!
      intent(in)  :: xxpar,kj
      intent(out) :: accg 
!
      accg=0.
      do ks=1,nspar
        if (kj/=ks) then
!
          e1=xxpar(1);e10=fsp(ks,ixp)
          e2=xxpar(2);e20=fsp(ks,iyp)
          e3=xxpar(3);e30=fsp(ks,izp)
!
!  Get the distances in each ortogonal component. 
!  These are NOT (x,y,z) for all.
!  For cartesian it is (x,y,z), for cylindrical (s,phi,z)
!  for spherical (r,theta,phi)
! 
          if (lcartesian_coords) then
            evr(1) = e1 - e10  
            evr(2) = e2 - e20  
            evr(3) = e3 - e30  
          elseif (lcylindrical_coords) then
            evr(1) = e1 - e10*cos(e2-e20)  
            evr(2) = e10*sin(e2-e20)       
            evr(3) = e3 - e30              
          elseif (lspherical_coords) then
            call fatal_error("loop_through_sinks","not yet implemented "//&
                 "for spherical polars")
          endif
!
!  Particles relative distance from each other
!
!  r_ij = sqrt(ev1**2 + ev2**2 + ev3**2)
!  invr3_ij = r_ij**(-3)
!
          r2_ij = sum(evr**2)+r_smooth(ks)**2
!
          invr3_ij = r2_ij**(-1.5)
!
!  Gravitational acceleration: g=g0/|r-r0|^3 (r-r0)
!  The acceleration is in non-coordinate basis (all have dimension of length). 
!  The main dxx_dt of particle_dust takes care of 
!  transforming the linear velocities to angular changes 
!  in position.
!
          accg=accg-GNewton*pmass(ks)*invr3_ij*evr(1:3)
!
        endif
      enddo
!
    endsubroutine loop_through_sinks
!**********************************************************
    subroutine point_par_name(a,iname)
!
!  Register a, a simple scalar, as diagnostic.
!  Works for individual particle diagnostics. 
!
!  07-mar-08/wlad: adapted from sum_par_name
!
      use Cdata
!
      real ::  a
      integer :: iname
!
      if (iname/=0) then
        fname(iname)=a
      endif
!
!  There is no entry in itype_name.
!
    endsubroutine point_par_name
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
    subroutine reset_center_of_mass(dfp)
!
!  If the center of mass was accelerated, reset its position
!  to the center of the grid. Must be called by a dxxp_dt_nbody?
!
!  Assumes that the total mass of the particles is one.
!
!  27-aug-06/wlad: coded
!
      use Mpicomm,only:stop_it
!
      real, dimension(mpar_loc,mpvar),intent(inout) :: dfp
      real, dimension(3) :: vcm
      integer :: k
!
        vcm(1) = sum(pmass*fsp(:,ivpx))
        vcm(2) = sum(pmass*fsp(:,ivpy))
        vcm(3) = sum(pmass*fsp(:,ivpz))
!
      do k=1,npar_loc
        if (ipar(k)<=nspar) then
          dfp(k,ixp:izp) = dfp(k,ixp:izp) - vcm*totmass1
        endif
      enddo
!
    endsubroutine reset_center_of_mass
!***********************************************************************
    subroutine integrate_selfgravity(p,rrp,xxpar,accg,rp0)
!
! Calculates acceleration on the point (x,y,z)=xxpar
! due to the gravity of the gas+dust.
!
! 15-sep-06/wlad : coded
!
      use Mpicomm
!
      real, dimension(nx,3) :: dist
      real, dimension(nx) :: rrp,selfgrav,density
      real, dimension(nx) :: dv,jac,dqy,tmp
      real :: dqx,dqz,rp0,fac
      real, dimension(3) :: xxpar,accg,sum_loc
      integer :: k,j
      type (pencil_case) :: p      
!
      intent(out) :: accg
!
      if (headtt) &
           print*,'Adding gas+dust gravity to the massive particles'
!
      if (coord_system=='cartesian') then
        jac=1.;dqx=dx;dqy=dy;dqz=dz
        dist(:,1)=x(l1:l2)-xxpar(1)
        dist(:,2)=y(  m  )-xxpar(2)
        dist(:,3)=z(  n  )-xxpar(3)
      elseif (coord_system=='cylindric') then
        jac=x(l1:l2);dqx=dx;dqy=x(l1:l2)*dy;dqz=dz
        dist(:,1)=x(l1:l2)-xxpar(1)*cos(y(m)-xxpar(2))
        dist(:,2)=         xxpar(2)*sin(y(m)-xxpar(2))
        dist(:,3)=z(  n  )-xxpar(3)
      elseif (coord_system=='spherical') then
        call stop_it("integrate_selfgravity: "//&
             " not yet implemented for spherical polars")
      endif

      if (nzgrid==1) then
        dv=dqx*dqy
      else
        dv=dqx*dqy*dqz
      endif
!
! The gravity of every single cell - should exclude inner and outer radii...
!
! selfgrav = G*((rho+rhop)*dv)*mass*r*(r**2 + r0**2)**(-1.5)
! gx = selfgrav * r\hat dot x\hat
! -> gx = selfgrav * (x-x0)/r = G*((rho+rhop)*dv)*mass*(r**2+r0**2)**(-1.5) * (x-x0)
!
      density=0.
      if (npar>nspar) then !npar>nspar equals dust is being used
        density=density+p%rhop
      else
        density=density+p%rho
      endif
      selfgrav = prhs_cte*density*jac*dv*(rrp**2 + rp0**2)**(-1.5)
!
! Everything inside the accretion radius of the particle should
! not exert gravity (numerical problems otherwise)
!
      where (rrp.le.rp0)
        selfgrav = 0
      endwhere
!
! Exclude the frozen zones
!
      if (lexclude_frozen) then
        if (lcylinder_in_a_box) then
          where ((p%rcyl_mn.le.r_int).or.(p%rcyl_mn.ge.r_ext))
            selfgrav = 0
          endwhere
        else
          where ((p%r_mn.le.r_int).or.(p%r_mn.ge.r_ext))
            selfgrav = 0
          endwhere
        endif
      endif
!
! Integrate the accelerations on this processor
! And sum over processors with mpireduce
!
      do j=1,3
        tmp=selfgrav*dist(:,j)
        !take proper care of the trapezoidal rule
        !in the case of non-periodic boundaries
        fac = 1.
        if ((m==m1.and.ipy==0).or.(m==m2.and.ipy==nprocy-1)) then
          if (.not.lperi(2)) fac = .5*fac
        endif
!
        if (lperi(1)) then
          sum_loc(j) = fac*sum(tmp)
        else
          sum_loc(j) = fac*(sum(tmp(2:nx-1))+.5*(tmp(1)+tmp(nx)))
        endif
        call mpireduce_sum_scl(sum_loc(j),accg(j))
      enddo
!
!  Broadcast particle acceleration
!
      call mpibcast_real(accg,3)
!
    endsubroutine integrate_selfgravity
!***********************************************************************
    subroutine share_sinkparticles(fp)
!
! Broadcast sink particles across processors
! The non-root processors, if they find a sink particle in their
! fp array, they:
!    send it to root with a true logical
!    else send a false logical
!
! The root processor receives all the logicals, and if
! they are true, then receives the value, broadcasting it
!
! 07-sep-06/wlad: coded
!
      use Mpicomm
!
      real, dimension(mpar_loc,mpvar) :: fp
      logical, dimension(nspar) :: lsink
      integer :: ks,k,tagsend,j,tagrecv
!
      if (lmpicomm) then
!
! Loop through the sink particles
!
        do ks=1,nspar
!
! Set the logical to false initially
!
          lsink(ks)=.false.
!
! Loop through the particles on this processor
!
          do k=1,npar_loc
            if (ks==ipar(k)) then
!
! A sink was found here. Turn the logical true and copy fp to fsp
!
              lsink(ks) = .true.
              fsp(ks,:) = fp(k,:)
!
! Send it to root. As there will be just nspar calls to
! mpisend, the tag can be ipar itself
!
              if (.not.lroot) then
                call mpisend_real(fsp(ks,:),mpvar,root,ks)
                if (ip<=6) print*,'logical for particle ',ks,&
                     ' set to true on processor ',iproc, &
                     ' with tag=',ks
              endif
            endif
          enddo
!
! Send the logicals from each processor. As all processors send nspar calls,
! the tags are now in base nspar. It assures that two logicals will not have
! the same tag.
!
          if (.not.lroot) then
            tagsend = nspar*iproc + ks
            call mpisend_logical(lsink(ks),1,root,tagsend)
          else
!
! The root receives all logicals. Same tag.
!
            do j=1,ncpus-1
              tagrecv = nspar*j + ks
              call mpirecv_logical(lsink(ks),1,j,tagrecv)
!
! Test the received logicals
!
              if (lsink(ks)) then
!
! Found a sink particle. Get the value of fsp
!
                call mpirecv_real(fsp(ks,:),mpvar,j,ks)
                if (ip<=6) print*,'logical for particle ',ks,&
                     ' is true on processor ',j, &
                     ' with tag=',ks,' on root'
              endif
            enddo
          endif
!
! Broadcast the received fsp
!
          call mpibcast_real(fsp(ks,:),mpvar)
!
! Print the result in all processors
!
          if (ip<=8)  print*,'share_sinkparticles: finished loop. '//&
               'sink particles in proc ',iproc,&
               ' are fsp(ks,:)=',fsp(ks,:)
        enddo
      else
!
! Non-parallel. Just copy fp to fsp
!
        do k=1,npar_loc
          if (ipar(k)<=nspar) fsp(ipar(k),:) = fp(k,:)
        enddo
!
      endif
!
    endsubroutine share_sinkparticles
!***********************************************************************
    subroutine get_totalmass(tmass)
!
! called from global_shear to set the initial keplerian field
!
      real :: tmass
!
      if (lnorm) then 
        tmass=1.
      else
        tmass=sum(pmass)
      endif
!   
    endsubroutine get_totalmass
!***********************************************************************
    subroutine get_gravity_field_nbody(grr,gg,ks)

      use Cdata,only: coord_system,x,y,z,l1,l2,&
                      m,n,nx,lcylindrical_gravity
      use Mpicomm,only:stop_it
!
      real, dimension(nx),intent(in) :: grr
      real, dimension(nx,3),intent(out) :: gg
      real :: x0,y0,z0
      integer :: ks
!
      x0=fsp(ks,ixp);y0=fsp(ks,iyp);z0=fsp(ks,izp)
!
      if (coord_system=='cartesian') then
        gg(:,1) = (x(l1:l2)-x0)*grr
        gg(:,2) = (y(  m  )-y0)*grr
        gg(:,3) = (z(  n  )-z0)*grr
      elseif (coord_system=='cylindric') then
        gg(:,1) = (x(l1:l2)-x0*cos(y(m)-y0))*grr
        gg(:,2) =           x0*sin(y(m)-y0) *grr
        gg(:,3) = (z(  n  )-z0)             *grr
      elseif (coord_system=='spherical') then
        call stop_it("get_gravity_field_nbody: off-center gravity"//&
             " field not yet implemented for spherical polars")
      endif
!
    endsubroutine get_gravity_field_nbody
!***********************************************************************
    subroutine calc_torque(p,dist,ks)
!
!  Output torque diagnostic for sink particle ks
!
!  05-nov-05/wlad : coded
!
      use Sub
      use Mpicomm, only: stop_it
!
      type (pencil_case) :: p
      real, dimension(nx) :: torque,torqint,torqext
      real, dimension(nx) :: dist,rpre
      real :: rr,w2,smap,hills
      integer :: ks,i
!
      if (.not.lcartesian_coords) &
           call stop_it("particles_nbody: calc_torque not yet "//&
           "implemented for curvilinear coordinates.")
!
      if (ks==nspar) &
           call stop_it("Nonsense to calculate torques for the star")

      rr    = sqrt(fsp(ks,ixp)**2 + fsp(ks,iyp)**2 + fsp(ks,izp)**2)
      w2    = fsp(ks,ivpx)**2 + fsp(ks,ivpy)**2 + fsp(ks,ivpz)**2
      smap  = 1./(2./rr - w2)
      hills = smap*(pmass(ks)*pmass1(nspar)/3.)**(1./3.)
!
      rpre  = fsp(ks,ixp)*y(m) - fsp(ks,iyp)*x(l1:l2)
      torque = GNewton*pmass(ks)*p%rho*rpre*&
           (dist**2 + r_smooth(ks)**2)**(-1.5)
!
      torqext=0.
      torqint=0.
!
      do i=1,nx
!
! Exclude material from inside the Roche Lobe
!
        if (dist(i).ge.hills) then
!
! External torque
!
          if (p%rcyl_mn(i).ge.rr) torqext(i) = torque(i)
!
! Internal torque
!
          if (p%rcyl_mn(i).le.rr) torqint(i) = torque(i)
!
        endif
!
      enddo
!
      call sum_lim_mn_name(torqext,idiag_torqext(ks),p)
      call sum_lim_mn_name(torqint,idiag_torqint(ks),p)
!
    endsubroutine calc_torque
!***********************************************************************
    subroutine get_ramped_mass
!     
      real :: ramping_period,tmp
      integer :: ks
!
      ramping_period=2*pi*ramp_orbits
      if (t .le. ramping_period) then
        !sin ((pi/2)*(t/(ramp_orbits*2*pi))
        tmp=0.
        do ks=1,nspar-1
          pmass(ks)= max(tini,&
               final_ramped_mass(ks)*(sin((.5*pi)*(t/ramping_period))**2))
          tmp=tmp+pmass(ks)
        enddo
        pmass(nspar)= 1-tmp
      else
        pmass=final_ramped_mass
      endif
!
    endsubroutine get_ramped_mass
!***********************************************************************
    subroutine calc_nbodygravity_particles(f)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(nx,3) :: ggt
!
      if (linterpolate_gravity) then
!
        if (lramp) call get_ramped_mass
!
! Calculate grid - sink particles distances
!
        do n=n1,n2
          do m=m1,m2
            call get_total_gravity(ggt)
            f(l1:l2,m,n,iglobal_ggp:iglobal_ggp+2)=ggt
          enddo
        enddo
!
!  else do nothing
!
      endif
!
    endsubroutine calc_nbodygravity_particles
!***********************************************************************
      subroutine get_total_gravity(ggt)
!  
        use Sub
!
      real, dimension (nx,nspar) :: rp_mn,rpcyl_mn
      real, dimension (nx,3)     :: ggp,ggt
      real, dimension (nx)       :: grav_particle,rrp
      integer                    :: ks
!
      intent(out) :: ggt
!
      ggt=0.
      do ks=1,nspar
!
! Spherical and cylindrical distances
!
        call get_radial_distance(rp_mn(:,ks),rpcyl_mn(:,ks),&
             e1_=fsp(ks,ixp),e2_=fsp(ks,iyp),e3_=fsp(ks,izp))
!
! Check which particle has cylindrical gravity switched on
!
        if (lcylindrical_gravity_nbody(ks)) then
          rrp = rpcyl_mn(:,ks)
        else
          rrp = rp_mn(:,ks)
        endif
!
! Gravity field from the particle ks
!
        grav_particle =-GNewton*pmass(ks)*(rrp**2+r_smooth(ks)**2)**(-1.5)
        call get_gravity_field_nbody(grav_particle,ggp,ks)
!
        if ((ks==nspar).and.lnogravz_star) &
             ggp(:,3) = 0.
!
! Sum up the accelerations of the sinks
!
        ggt=ggt+ggp
!
      enddo
!
    endsubroutine get_total_gravity
!***********************************************************************
    subroutine rprint_particles_nbody(lreset,lwrite)
!
!  Read and register print parameters relevant for sink particles.
!
!  17-nov-05/anders+wlad: adapted
!
      use Cdata
      use Sub, only: parse_name
      use General, only: chn
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      integer :: iname,ks,j
      character :: str
      character (len=5) :: sks
!
!  Write information to index.pro
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset
!
      if (lreset) then
        idiag_xxspar=0;idiag_vvspar=0
        idiag_torqint=0;idiag_torqext=0
        idiag_totenergy=0;idiag_totangmom=0
      endif
!
!  Run through all possible names that may be listed in print.in
!
      if (lroot .and. ip<14) print*,'rprint_particles: run through parse list'
!
!  Now check diagnostics for specific particles
!
      do ks=1,nspar
        call chn(ks,sks)
        do j=1,3
          if (j==1) str='x';if (j==2) str='y';if (j==3)  str='z'
          do iname=1,nname
            call parse_name(iname,cname(iname),cform(iname),&
                 trim(str)//'par'//trim(sks),idiag_xxspar(ks,j))
            call parse_name(iname,cname(iname),cform(iname),&
                 'v'//trim(str)//'par'//trim(sks),idiag_vvspar(ks,j))
          enddo
!
!  Run through parse list again
!
          if (lwr) then
            write(3,*) ' i_'//trim(str)//'par'//trim(sks)//'=',&
                 idiag_xxspar(ks,j)
            write(3,*) 'i_v'//trim(str)//'par'//trim(sks)//'=',&
                 idiag_vvspar(ks,j)
          endif
!
        enddo
!
        do iname=1,nname
          call parse_name(iname,cname(iname),cform(iname),&
               'torqint_'//trim(sks),idiag_torqint(ks))
          call parse_name(iname,cname(iname),cform(iname),&
               'torqext_'//trim(sks),idiag_torqext(ks))
        enddo
!
        if (lwr) then
          write(3,*) 'i_torqint_'//trim(sks)//'=',idiag_torqint(ks)
          write(3,*) 'i_torqext_'//trim(sks)//'=',idiag_torqext(ks)
        endif
      enddo
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),&
             'totenergy',idiag_totenergy)
        call parse_name(iname,cname(iname),cform(iname),&
             'totangmom',idiag_totangmom)
      enddo
!
     if (lwr) then
       write(3,*) 'i_totenergy=',idiag_totenergy
       write(3,*) 'i_totangmom=',idiag_totangmom
     endif
!
   endsubroutine rprint_particles_nbody
!***********************************************************************
 endmodule Particles_nbody
