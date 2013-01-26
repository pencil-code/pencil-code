! $Id$
!
!  This module takes care of direct N-body gravity between particles.
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
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_sub
!
  implicit none
!
  include 'particles_nbody.h'
!
  real, dimension(nspar,mspvar) :: fsp
  real, dimension(nspar) :: xsp0=0.0, ysp0=0.0, zsp0=0.0
  real, dimension(nspar) :: vspx0=0.0, vspy0=0.0, vspz0=0.0
  real, dimension(nspar) :: pmass=0.0, r_smooth=0.0, pmass1
  real, dimension(nspar) :: accrete_hills_frac=0.2, final_ramped_mass=0.0
  real :: delta_vsp0=1.0, totmass, totmass1
  real :: create_jeans_constant=0.25, GNewton1
  real :: GNewton=impossible, prhs_cte
  real :: cdtpnbody=0.1
  real, pointer :: rhs_poisson_const, tstart_selfgrav
  integer :: ramp_orbits=5, mspar_orig=1
  integer :: iglobal_ggp=0, istar=1, imass=0
  integer :: maxsink=10*nspar, icreate=100
  logical, dimension(nspar) :: lcylindrical_gravity_nbody=.false.
  logical, dimension(nspar) :: lfollow_particle=.false., laccretion=.false.
  logical, dimension(nspar) :: ladd_mass=.false.
  logical :: lcalc_orbit=.true., lbackreaction=.false., lnorm=.true.
  logical :: lreset_cm=.false., lnogravz_star=.false., lexclude_frozen=.false.
  logical :: lnoselfgrav_star=.true.
  logical :: lramp=.false., lcreate_sinks=.false., lcreate_gas=.true.
  logical :: ldt_nbody=.false., lcreate_dust=.true.
  logical :: linterpolate_gravity=.false., linterpolate_linear=.true.
  logical :: linterpolate_quadratic_spline=.false.
  logical :: laccrete_when_create=.true.
  logical :: ldust=.false.
  character (len=labellen) :: initxxsp='random', initvvsp='nothing'
!
  namelist /particles_nbody_init_pars/ &
      initxxsp, initvvsp, xsp0, ysp0, zsp0, vspx0, vspy0, vspz0, delta_vsp0, &
      pmass, r_smooth, lcylindrical_gravity_nbody, lexclude_frozen, GNewton, &
      bcspx, bcspy, bcspz, ramp_orbits, lramp, final_ramped_mass, prhs_cte, &
      linterpolate_gravity, linterpolate_quadratic_spline, laccretion, &
      accrete_hills_frac, istar, maxsink, lcreate_sinks, icreate, lcreate_gas, &
      lcreate_dust, ladd_mass, laccrete_when_create, ldt_nbody, cdtpnbody
!
  namelist /particles_nbody_run_pars/ &
      dsnap_par_minor, linterp_reality_check, lcalc_orbit, lreset_cm, &
      lnogravz_star, lfollow_particle, lbackreaction, lexclude_frozen, &
      GNewton, bcspx, bcspy, bcspz,prhs_cte, lnoselfgrav_star, &
      linterpolate_quadratic_spline, laccretion, accrete_hills_frac, istar, &
      maxsink, lcreate_sinks, icreate, lcreate_gas, lcreate_dust, ladd_mass, &
      laccrete_when_create, ldt_nbody, cdtpnbody
!
  integer, dimension(nspar,3) :: idiag_xxspar=0,idiag_vvspar=0
  integer, dimension(nspar)   :: idiag_torqint=0,idiag_torqext=0
  integer                     :: idiag_totenergy=0,idiag_totangmom=0
!
  contains
!***********************************************************************
    subroutine register_particles_nbody()
!
!  Set up indices for access to the f and fsp.
!
!  27-aug-06/wlad: adapted
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  Set up mass as particle index. Plus seven, since the other 6 are
!  used by positions and velocities.
!
      imass=npvar+1
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
      if (npar < nspar) then
        if (lroot) write(0,*) 'npar, nspar = ', npar, nspar
        call fatal_error('register_particles_nbody','the number of massive'//&
             ' particles (nspar) is less than the allocated number of particles'//&
             ' (npar). Increase npar to the minimum number (npar=npsar) needed'//&
             ' in cparam.local and recompile')
      endif
!
    endsubroutine register_particles_nbody
!***********************************************************************
    subroutine initialize_particles_nbody(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  27-aug-06/wlad: adapted
!
      use FArrayManager
      use SharedVariables
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      integer :: ierr,ks
!
!  Look for initialized masses.
!
      mspar_orig=0
      do ks=1,nspar
        if (pmass(ks)/=0.0) then
          mspar_orig=mspar_orig+1
          ipar_nbody(mspar_orig)=ks
        endif
      enddo
!
!  Just needs to do this starting, otherwise mspar will be overwritten after
!  reading the snapshot
!
      if (lstarting) then
        mspar=mspar_orig
      else
        !else read pmass from fsp
        pmass(1:mspar)=fsp(1:mspar,imass)
      endif
!
!  When first called, mspar was zero, so no diagnostic index was written to
!  index.pro
!
      if (lroot) open(3, file=trim(datadir)//'/index.pro', &
          STATUS='old', POSITION='append')
      call rprint_particles_nbody(.false.,LWRITE=lroot)
      if (lroot) close(3)
!
!  G_Newton. Overwrite the one set by start.in if set again here,
!  because I might want units in which both G and GM are 1.
!  If we are solving for selfgravity, get that one instead.
!
      if (GNewton == impossible) then
        !ONLY REASSIGN THE GNEWTON
        !IF IT IS NOT SET IN START.IN!!!!!
        if (lselfgravity) then
          call get_shared_variable('rhs_poisson_const',rhs_poisson_const,ierr)
          if (ierr/=0) then
            if (lroot) print*, 'initialize_particles_nbody: '// &
                'there was a problem when getting rhs_poisson_const!'
            call fatal_error('initialize_particles_nbody','')
          endif
          GNewton=rhs_poisson_const/(4*pi)
        else
          GNewton=G_Newton
        endif
      endif
!
      GNewton1=1./GNewton
!
!  Get the variable tstart_selfgrav from selfgrav, so we know when to create
!  sinks.
!
      if (lselfgravity) then
        call get_shared_variable('tstart_selfgrav',tstart_selfgrav,ierr)
        if (ierr/=0) then
          if (lroot) print*, 'initialize_particles_nbody: '// &
              'there was a problem when getting tstart_selfgrav!'
          call fatal_error('initialize_particles_nbody','')
        endif
      endif
!
!  inverse mass
!
      if (lstarting) then
        if (lramp) then
          do ks=1,mspar
            if (ks/=istar) pmass(ks) = epsi
          enddo
          pmass(istar)=1-epsi*(mspar-1)
        endif
      else
        !read imass from the snapshot
        pmass(1:mspar)=fsp(1:mspar,imass)
      endif
!
      pmass1=1./max(pmass,tini)
!
!  inverse total mass
!
      totmass=sum(pmass)
      totmass1=1./max(totmass, tini)
!
!  Check for consistency between the cylindrical gravities switches
!  from cdata and the one from nbody.
!
      if (((lcylindrical_gravity).and.&
        (.not.lcylindrical_gravity_nbody(istar))).or.&
             (.not.lcylindrical_gravity).and.&
             (lcylindrical_gravity_nbody(istar))) then
        call fatal_error('initialize_particles_nbody','inconsitency '//&
            'between lcylindrical_gravity from cdata and the '//&
            'one from nbody')
      endif
!
      if (rsmooth/=r_smooth(istar)) then
        print*,'rsmooth from cdata=',rsmooth
        print*,'r_smooth(istar)=',r_smooth(istar)
        call fatal_error('initialize_particles_nbody','inconsitency '//&
            'between rsmooth from cdata and the '//&
            'one from nbody')
      endif
!
!  Check for consistency between the poisson and integral
!  calculations of selfgravity.
!
      if (lbackreaction.and.lselfgravity) then
        print*,'self-gravity is already taken into '//&
            'account. lbackreaction is a way of '//&
            'integrating the self-gravity only in few '//&
            'specific points of the grid'
        call fatal_error('initialize_particles_nbody','')
      endif
!
!  The presence of dust particles needs to be known.
!
      if (npar > nspar) ldust=.true.
      if (linterpolate_gravity.and.(.not.ldust)) then
        if (lroot) print*,'interpolate gravity is just'//&
            ' for the dust component. No need for it if'//&
            ' you are not using dust particles'
        call fatal_error('initialize_particles_nbody','')
      endif
      if (any(laccretion).and.linterpolate_gravity) then
        if (lroot) print*,'interpolate gravity  not yet '//&
            'implemented in connection with accretion'
        call fatal_error('initialize_particles_nbody','')
      endif
!
      if (linterpolate_gravity) then
         if (lroot) print*,'initializing global array for nbody gravity'
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
        call fatal_error('initialize_particles_nbody','')
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_nbody
!***********************************************************************
    subroutine pencil_criteria_par_nbody()
!
!  All pencils that the Particles_nbody module depends on are specified here.
!
!  22-sep-06/wlad: adapted
!
      if (ldensity) lpenc_requested(i_rho)=.true.
      if (ldust.and.lbackreaction) &
          lpenc_requested(i_rhop)=.true.
!
      if (lexclude_frozen) then
        if (lcylinder_in_a_box) then
          lpenc_requested(i_rcyl_mn)=.true.
        else
          lpenc_requested(i_r_mn)=.true.
        endif
      endif
!
      if (idiag_totenergy/=0) then
        lpenc_diagnos(i_u2)=.true.
        lpenc_diagnos(i_rho)=.true.
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
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_par_nbody
!***********************************************************************
    subroutine calc_pencils_par_nbody(f,p)
!
!  Calculate nbody particle pencils
!
!  22-sep-06/wlad: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_par_nbody
!***********************************************************************
    subroutine init_particles_nbody(f,fp)
!
!  Initial positions and velocities of nbody particles.
!  Overwrite the position asserted by the dust module
!
!  17-nov-05/anders+wlad: adapted
!
      use General, only: random_number_wrapper
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension(mspar) :: kep_vel,sma
      real, dimension(mspar,3) :: position,velocity
      real :: tmp,parc
      integer :: k,ks
!
      intent (in) :: f
      intent (out) :: fp
!
!  Initialize particles' positions and velocities.
!
      if (mspar < nspar) then
        call warning("init_particles_nbody", 'mspar < nspar:')
        print*, 'mspar = ', mspar, ', nspar = ', nspar
        print*, 'nspar is the number of allocated n-body particles'
        print*, 'mspar is the number of those particles with'
        print*, ' non-zero mass in this processor.'
        print*, ''
        print*, 'you should set pmass to non-zero values in start.in'
      endif
      if (mspar > 0) then
        position(1:nspar,1) = xsp0
        position(:,2) = ysp0
        position(:,3) = zsp0
        velocity(:,1) = vspx0
        velocity(:,2) = vspy0
        velocity(:,3) = vspz0
      endif
!
      select case (initxxsp)
!
      case ('nothing')
        if (lroot) print*, 'init_particles_nbody: nothing'
!
      case ('origin')
        if (lroot) then
          print*, 'init_particles_nbody: All nbody particles at origin'
          fp(1:mspar,ixp:izp)=0.0
        endif
!
      case ('constant')
        if (lroot) &
            print*, 'init_particles_nbody: Place nbody particles at x,y,z=', &
            xsp0, ysp0, zsp0
        do k=1,npar_loc
          if (ipar(k)<=mspar) fp(k,ixp:izp)=position(ipar(k),1:3)
        enddo
!
      case ('random')
        if (lroot) print*, 'init_particles_nbody: Random particle positions'
        do ks=1,mspar
          if (nxgrid/=1) call random_number_wrapper(position(ks,ixp))
          if (nygrid/=1) call random_number_wrapper(position(ks,iyp))
          if (nzgrid/=1) call random_number_wrapper(position(ks,izp))
        enddo
!
        if (nxgrid/=1) &
             position(1:mspar,ixp)=xyz0_loc(1)+position(1:mspar,ixp)*Lxyz_loc(1)
        if (nygrid/=1) &
             position(1:mspar,iyp)=xyz0_loc(2)+position(1:mspar,iyp)*Lxyz_loc(2)
        if (nzgrid/=1) &
             position(1:mspar,izp)=xyz0_loc(3)+position(1:mspar,izp)*Lxyz_loc(3)
!
!  Loop through ipar to allocate the nbody particles.
!
        do k=1,npar_loc
          if (ipar(k) <= mspar) then
            if (ldebug) &
                 print*,'initparticles_nbody. Slot for nbody particle ',&
                 ipar(k),' was at fp position ',k,' at processor ',iproc
!
            fp(k,ixp:izp)=position(ipar(k),1:3)
!
!  Correct for non-existing dimensions (not really needed, I think).
!
            if (nxgrid==1) fp(k,ixp)=x(nghost+1)
            if (nygrid==1) fp(k,iyp)=y(nghost+1)
            if (nzgrid==1) fp(k,izp)=z(nghost+1)
!
            if (ldebug) &
                 print*,'initparticles_nbody. Nbody particle ',&
                 ipar(k),' located at ixp=',fp(k,ixp)
          endif
        enddo
!
      case ('fixed-cm')
        if (lgrav) then
          print*,'a gravity module is being used. Are you using '//&
               'both a fixed central gravity and nbody gravity? '//&
               'better stop and check'
          call fatal_error('init_particles_nbody','')
        endif
!
!  Ok, I have the masses and the positions of all massive particles
!  except the last, which will have a position determined to fix the
!  center of mass on the center of the grid.
!
        if (any(ysp0/=0)) then
          if (lspherical_coords) then
            call fatal_error('init_particles_nbody','not yet generalized'//&
                 ' for non-zero initial inclinations')
          else
            call fatal_error('init_particles_nbody','not yet generalized'//&
                 ' for non-zero azimuthal initial position')
          endif
        endif
        if (any(zsp0/=0)) then
          if (lspherical_coords) then
            call fatal_error('init_particles_nbody','not yet generalized'//&
                 ' for non-zero azimuthal initial position')
          else
            call fatal_error('init_particles_nbody','nbody code not'//&
                 ' yet generalized to allow initial inclinations')
          endif
        endif
        if (lcylindrical_coords.or.lspherical_coords) then
          if (any(xsp0<0)) &
              call fatal_error('init_particles_nbody', &
              'in cylindrical coordinates '//&
              'all the radial positions must be positive')
        endif
!
        if (lroot) then
          print*,'init_particles_nbody: fixed-cm - mass and position arranged'
          print*,'                      so that the center of mass is at rest'
          print*,'                      at the center of the grid.'
        endif
!
        if (lspherical_coords) then
          if (lroot) print*,'put all particles in the midplane'
          position(1:mspar,iyp)=pi/2
        endif
!
        tmp = 0.;parc=0
        do ks=1,mspar
          if (ks/=istar) then
            sma(ks)=abs(position(ks,1))
            tmp=tmp+pmass(ks)
            parc = parc - sma(ks)*pmass(ks)
          endif
        enddo
!
!  Fixed-cm assumes that the total mass is always one. The mass of the
!  star is adjusted to ensure this.
!
        pmass(istar)=1.- tmp;pmass1=1./max(pmass,tini);totmass=1.;totmass1=1.
!
        parc = parc*totmass1
        if (tmp >= 1.0) &
            call fatal_error('init_particles_nbody', &
            'The mass of one '//&
            '(or more) of the particles is too big! The masses should '//&
            'never be bigger than g0. Please scale your assemble so that '//&
            'the combined mass of the (n-1) particles is less than that. '//&
            'The mass of the last particle in the pmass array will be '//&
            'reassigned to ensure that the total mass is g0')
!
        do ks=1,mspar
          if (ks/=istar) &
              position(ks,1)=sign(1.,position(ks,1))* (sma(ks) + parc)
        enddo
!
!  The last one (star) fixes the CM at Rcm=zero
!
        if (lcartesian_coords) then
          position(istar,1)=parc
        elseif (lcylindrical_coords) then
          !put the star in positive coordinates, with pi for azimuth
          position(istar,1)=abs(parc)
          position(istar,2)=pi
        elseif (lspherical_coords) then
          position(istar,1)=abs(parc)
          position(istar,3)=pi
        endif
!
        if (ldebug) then
          print*,'pmass =',pmass
          print*,'position (x)=',position(:,1)
          print*,'position (y)=',position(:,2)
          print*,'position (z)=',position(:,3)
        endif
!
!  Loop through ipar to allocate the nbody particles
!
        do k=1,npar_loc
          if (ipar(k) <= mspar) then
!
            if (ldebug) &
                 print*,'initparticles_nbody. Slot for nbody particle ',&
                 ipar(k),' was at fp position ',k,' at processor ',iproc
!
!  Here I substitute the first mspar dust particles by massive ones,
!  since the first ipars are less than mspar
!
            fp(k,ixp:izp)=position(ipar(k),1:3)
!
!  Correct for non-existing dimensions (not really needed, I think)
!
            if (nxgrid==1) fp(k,ixp)=x(nghost+1)
            if (nygrid==1) fp(k,iyp)=y(nghost+1)
            if (nzgrid==1) fp(k,izp)=z(nghost+1)
!
            if (ldebug) & 
                 print*,'init_particles_nbody. Nbody particle ',ipar(k),&
                 ' located at ixp=',fp(k,ixp)
          endif
        enddo
!
      case default
        if (lroot) print*,'init_particles_nbody: No such such value for'//&
            ' initxxsp: ',trim(initxxsp)
        call fatal_error('init_particles_nbody','')
!
      endselect
!
!  Redistribute particles among processors (now that positions are determined).
!
      call boundconds_particles(fp,ipar)
!
!  Initial particle velocity.
!
      select case (initvvsp)
!
      case ('nothing')
        if (lroot) print*, 'init_particles_nbody: No particle velocity set'
!
      case ('zero')
        if (lroot) then
          print*, 'init_particles_nbody: Zero particle velocity'
          fp(1:mspar,ivpx:ivpz)=0.
        endif
!
      case ('constant')
        if (lroot) then
          print*, 'init_particles_nbody: Constant particle velocity'
          print*, 'init_particles_nbody: vspx0, vspy0, vspz0=', vspx0, vspy0, vspz0
        endif
        do k=1,npar_loc
          if (ipar(k) <= mspar) &
              fp(k,ivpx:ivpz)=velocity(ipar(k),1:3)
        enddo
!
      case ('fixed-cm')
!
!  Keplerian velocities for the planets
!
        parc=0.
        do ks=1,mspar
          if (ks/=istar) then
            kep_vel(ks)=sqrt(1./sma(ks)) !circular velocity
            parc = parc - kep_vel(ks)*pmass(ks)
          endif
        enddo
        parc = parc*totmass
        do ks=1,mspar
          if (ks/=istar) then
            if (lcartesian_coords) then
              velocity(ks,2) = sign(1.,position(ks,1))*(kep_vel(ks) + parc)
            elseif (lcylindrical_coords) then
              !positive for the planets
              velocity(ks,2) = abs(kep_vel(ks) + parc)
            elseif (lspherical_coords) then
              velocity(ks,3) = abs(kep_vel(ks) + parc)
            endif
          endif
        enddo
!
!  The last one (star) fixes the CM also with velocity zero
!
        if (lcartesian_coords) then
          velocity(istar,2)=parc
        elseif (lcylindrical_coords) then
          velocity(istar,2)=-parc
        elseif (lspherical_coords) then
          velocity(istar,3)=-parc
        endif
!
!  Loop through ipar to allocate the nbody particles
!
         do k=1,npar_loc
           if (ipar(k)<=mspar) then
             if (ldebug) & 
                  print*,&
                  'init_particles_nbody. Slot for nbody particle ',ipar(k),&
                  ' was at fp position ', k, ' at processor ', iproc
!
             fp(k,ivpx:ivpz) = velocity(ipar(k),1:3)
!
             if (ldebug) &
                  print*,'init_particles_nbody. Nbody particle ',ipar(k),&
                  ' has velocity y=',fp(k,ivpy)
!
           endif
         enddo
!
      case default
        if (lroot) print*, 'init_particles_nbody: No such such value for'//&
             'initvvsp: ',trim(initvvsp)
        call fatal_error('init_particles_nbody','')
!
      endselect
!
!  Make the particles known to all processors
!
      call boundconds_particles(fp,ipar)
      call bcast_nbodyarray(fp)
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_particles_nbody
!***********************************************************************
    subroutine dvvp_dt_nbody_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Add gravity from the particles to the gas and the backreaction from the
!  gas onto the particles.
!
!  Adding the gravity on the dust component via the grid is less accurate,
!  but much faster than looping through all npar_loc particle and add the
!  gravity of all massive particles to it.
!
!  07-sep-06/wlad: coded
!
      use Diagnostics
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (nx,mspar) :: rp_mn, rpcyl_mn
      real, dimension (mx,3) :: ggt
      real, dimension (nx) :: pot_energy
      real, dimension (3) :: xxpar,accg
      integer :: ks, k
      logical :: lintegrate, lparticle_out
!
      intent (in) :: f, p, fp, ineargrid
      intent (inout) :: df, dfp
!
!  Get the total gravity field. In the case of dust, it is already
!  pre-calculated
!
      if (linterpolate_gravity) then
!
!  Interpolate the gravity to the position of the particles.
!
        if (npar_imn(imn)/=0) then
          do k=k1_imn(imn),k2_imn(imn) !loop through the pencil
            if (ipar(k)>mspar) then !for dust
              !interpolate the gravity
              if (linterpolate_linear) then
                call interpolate_linear(f,iglobal_ggp,&
                    iglobal_ggp+2,fp(k,ixp:izp),accg,ineargrid(k,:),0,ipar(k))
              else if (linterpolate_quadratic_spline) then
!
!  This interpolation works also for cylindrical coordinates.
!
                call interpolate_quadratic_spline(f,iglobal_ggp,&
                    iglobal_ggp+2,fp(k,ixp:izp),accg,ineargrid(k,:),0,ipar(k))
              endif
              dfp(k,ivpx:ivpz)=dfp(k,ivpx:ivpz)+accg
            endif
          enddo
        endif
      endif
!
!  Add the acceleration to the gas.
!
      if (lhydro) then
        if (linterpolate_gravity) then
          !already calculated, so no need to lose time
          !calculating again
          ggt=f(:,m,n,iglobal_ggp:iglobal_ggp+2)
        else
          call get_total_gravity(ggt)
        endif
        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + ggt(l1:l2,:)
!
!  Backreaction of the gas+dust gravity onto the massive particles
!  The integration is needed for these two cases:
!
!   1. We are not solving the Poisson equation (lbackreaction)
!   2. We are, but a particle is out of the box (a star, for instance)
!      and therefore the potential cannot be interpolated.
!
        do ks=1,mspar
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
!  Sometimes making the star feel the selfgravity of the disk leads to
!  numerical troubles as the star is too close to the origin (in cylindrical
!  coordinates).
!
          if ((ks==istar).and.lnoselfgrav_star) lintegrate=.false.
!
          if (lintegrate) then
!
!  Get the acceleration particle ks suffers due to self-gravity.
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
!  Add it to its dfp
!
            do k=1,npar_loc
              if (ipar(k)==ks) &
                  dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + accg(1:3)
            enddo
          endif
!
!  Diagnostic
!
          if (ldiagnos) then
            if (idiag_totenergy/=0) pot_energy=0.0
            call get_radial_distance(rp_mn(:,ks),rpcyl_mn(:,ks),&
                E1_=fsp(ks,ixp),E2_=fsp(ks,iyp),E3_=fsp(ks,izp))
!
!  Calculate torques for output, if needed
!
            if ((idiag_torqext(ks)/=0).or.(idiag_torqint(ks)/=0)) &
                call calc_torque(p,rpcyl_mn(:,ks),ks)
!
!  Total energy
!
            if (idiag_totenergy/=0) then
              !potential energy
              pot_energy = pot_energy - &
                  GNewton*pmass(ks)*&
                  (rpcyl_mn(:,ks)**2+r_smooth(ks)**2)**(-0.5)
              if (ks==mspar) &
                  call sum_lim_mn_name(.5*p%rho*p%u2 + pot_energy,&
                  idiag_totenergy,p)
            endif
          endif
        enddo
      endif !if hydro
!
    endsubroutine dvvp_dt_nbody_pencil
!***********************************************************************
    subroutine dxxp_dt_nbody(dfp)
!
!  If the center of mass of the nbody particles was moved from the
!  center of the grid, reset it.
!
!  22-sep-06/wlad: coded
!
      real, dimension (mpar_loc,mpvar) :: dfp
!
      if (lreset_cm) call reset_center_of_mass(dfp)
!
    endsubroutine dxxp_dt_nbody
!***********************************************************************
    subroutine dvvp_dt_nbody(f,df,fp,dfp,ineargrid)
!
!  Evolution of nbody and dust particles velocities due to particle-particle
!  interaction only.
!
!  Coriolis and shear are already added in particles_dust.
!
!  27-aug-06/wlad: coded
!
      use Mpicomm, only: mpibcast_real
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension(mspar) :: sq_hills
      real :: rr, w2, sma2
      integer :: k, ks, j, jpos, jvel
      logical :: lheader, lfirstcall=.true.
!
      intent (in) :: f, ineargrid
      intent (inout) :: fp, df, dfp
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
!
!  Identify module.
!
      if (lheader) print*, 'dvvp_dt_nbody: Calculate dvvp_dt_nbody'
!
!  Evolve massive particle positions due to the gravity of other massive
!  particles. The gravity of the massive particles on the dust will be added
!  inside the pencil in dvvp_dt_nbody_pencil.
!
      if (lramp) call get_ramped_mass
!
!  Calculate Hills radius if laccretion is switched on.
!
      do ks=1,mspar
        if (laccretion(ks).and.(ks/=istar)) then
          if (lcartesian_coords) then
            rr=sqrt(fsp(ks,ixp)**2 + fsp(ks,iyp)**2 + fsp(ks,izp)**2)
          elseif (lcylindrical_coords) then
            rr=fsp(ks,ixp)
            if (nzgrid/=1) rr=sqrt(fsp(ks,ixp)**2+fsp(ks,izp)**2)
          elseif (lspherical_coords) then
            rr=fsp(ks,ixp)
          else
            call fatal_error('dvvp_dt_nbody','wrong coord system')
            rr=0.0
          endif
          !particle velocities are non-coordinate (linear)
          w2    = fsp(ks,ivpx)**2 + fsp(ks,ivpy)**2 + fsp(ks,ivpz)**2
          !squared semi major axis - assumes GM=1, so beware...
          sma2  = (rr/(2-rr*w2))**2
          !squared hills radius
          sq_hills(ks)=sma2*(pmass(ks)*pmass1(istar)/3)**(2./3.)
        else
          sq_hills(ks)=0.0
        endif
      enddo
!
!  Add the gravity from all N-body particles.
!
      do k=npar_loc,1,-1
!
        if (linterpolate_gravity) then
          !only loop through massive particles
          if (ipar(k)<=mspar) then
            call loop_through_nbodies(fp,dfp,k,sq_hills,ineargrid)
          endif
        else
          !for all particles
          call loop_through_nbodies(fp,dfp,k,sq_hills,ineargrid)
        endif
      enddo
!
!  Position and velocity diagnostics (per nbody particle).
!
      if (ldiagnos) then
        do ks=1,mspar
          if (lfollow_particle(ks)) then
            do j=1,3
              jpos=j+ixp-1 ; jvel=j+ivpx-1
              if (idiag_xxspar(ks,j)/=0) then
                call point_par_name(fsp(ks,jpos),idiag_xxspar(ks,j))
              endif
              if (idiag_vvspar(ks,j)/=0) &
                  call point_par_name(fsp(ks,jvel),idiag_vvspar(ks,j))
            enddo
          endif
        enddo
      endif
!
      if (lheader) print*,'dvvp_dt_nbody: Finished dvvp_dt_nbody'
      if (lfirstcall) lfirstcall=.false.
!
      call keep_compiler_quiet(f,df)
!
    endsubroutine dvvp_dt_nbody
!************************************************************
    subroutine loop_through_nbodies(fp,dfp,k,sq_hills,ineargrid)
!
!  Subroutine that adds the gravity from all massive particles.
!
!  07-mar-08/wlad:coded
!
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer :: k
      real, dimension (mspar) :: sq_hills
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: evr
      real :: r2_ij, rs2, invr3_ij, v_ij
      integer :: ks
!
      intent(inout) :: fp, dfp
      intent(in)  :: k
!
      do ks=1,mspar
        if (ipar(k)/=ks) then
!
          call get_particles_interdistance(fp(k,ixp:izp),fsp(ks,ixp:izp), &
              VECTOR=evr,DISTANCE=r2_ij,lsquare=.true.)
!
!  Particles relative distance from each other:
!
!  r_ij = sqrt(ev1**2 + ev2**2 + ev3**2)
!
!  invr3_ij = r_ij**(-3)
!
          rs2=(accrete_hills_frac(ks)**2)*sq_hills(ks)
!
!  If there is accretion, remove the accreted particles from the simulation.
!
          if (laccretion(ks) .and. (r2_ij<=rs2)) then
            call remove_particle(fp,ipar,k,dfp,ineargrid,ks)
            !add mass of the removed particle to the accreting particle
            if (ladd_mass(ks)) pmass(ks)=pmass(ks)+mp_swarm
            goto 99
          else
            r2_ij=r2_ij+r_smooth(ks)**2
            if (r2_ij > 0) then
              invr3_ij = r2_ij**(-1.5)
            else                ! can happen during pencil_check
              invr3_ij = 0.0
            endif
!
!  Gravitational acceleration: g=g0/|r-r0|^3 (r-r0)
!
!  The acceleration is in non-coordinate basis (all have dimension of length).
!  The main dxx_dt of particle_dust takes care of transforming the linear
!  velocities to angular changes in position.
!
            dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) - &
                GNewton*pmass(ks)*invr3_ij*evr(1:3)
!
!  Time-step constraint from N-body particles. We use both the criterion
!  that the distance to the N-body particle must not change too much in
!  one time-step and additionally we use the free-fall time-scale.
!
            if (lfirst.and.ldt.and.ldt_nbody) then
              v_ij=sqrt(sum((fp(k,ivpx:ivpz)-fp(ks,ivpx:ivpz))**2))
              dt1_max(ineargrid(k,1)-nghost)= &
                  max(dt1_max(ineargrid(k,1)-nghost),v_ij/sqrt(r2_ij)/cdtpnbody)
              dt1_max(ineargrid(k,1)-nghost)= &
                  max(dt1_max(ineargrid(k,1)-nghost), &
                  sqrt(GNewton*pmass(ks)*invr3_ij)/cdtpnbody)
            endif
!
          endif !if accretion
!
        endif !if (ipar(k)/=ks)
!
      enddo !nbody loop
!
99    continue
!
    endsubroutine loop_through_nbodies
!**********************************************************
    subroutine point_par_name(a,iname)
!
!  Register a, a simple scalar, as diagnostic.
!
!  Works for individual particle diagnostics.
!
!  07-mar-08/wlad: adapted from sum_par_name
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
!  18-mar-08/wlad: cylindrical and spherical corrections
!
      real, dimension(mpar_loc,mpvar),intent(inout) :: dfp
      real, dimension(mspar,mspvar) :: ftmp
      real, dimension(3) :: vcm
      real :: xcm,ycm,zcm,thtcm,phicm,vxcm,vycm,vzcm
      integer :: k
!
      ftmp=fsp(1:mspar,:)
!
      if (lcartesian_coords) then
        vcm(1) = sum(ftmp(:,imass)*ftmp(:,ivpx))
        vcm(2) = sum(ftmp(:,imass)*ftmp(:,ivpy))
        vcm(3) = sum(ftmp(:,imass)*ftmp(:,ivpz))
      else if (lcylindrical_coords) then
!
!  This is not really functional for other grids than Cartesian.
!
        xcm=sum(ftmp(:,imass)*(ftmp(:,ipx)*cos(ftmp(:,iyp))))
        ycm=sum(ftmp(:,imass)*(ftmp(:,ipx)*sin(ftmp(:,iyp))))
        phicm=atan2(ycm,xcm)
        vxcm=sum(ftmp(:,imass)*(&
            ftmp(:,ivpx)*cos(ftmp(:,iyp))-ftmp(:,ivpy)*sin(ftmp(:,iyp))))
        vycm=sum(ftmp(:,imass)*(&
            ftmp(:,ivpx)*sin(ftmp(:,iyp))+ftmp(:,ivpy)*cos(ftmp(:,iyp))))
        vcm(1)= vxcm*cos(phicm) + vycm*sin(phicm)
        vcm(2)=-vxcm*sin(phicm) + vycm*cos(phicm)
        vcm(3) = sum(ftmp(:,imass)*ftmp(:,ivpz))
      else if (lspherical_coords) then
        vxcm=sum(ftmp(:,imass)*( &
              ftmp(:,ivpx)*sin(ftmp(:,iyp))*cos(ftmp(:,izp))&
             +ftmp(:,ivpy)*cos(ftmp(:,iyp))*cos(ftmp(:,izp))&
             -ftmp(:,ivpz)*sin(ftmp(:,izp))                ))
        vycm=sum(ftmp(:,imass)*( &
              ftmp(:,ivpx)*sin(ftmp(:,iyp))*sin(ftmp(:,izp))&
             +ftmp(:,ivpy)*cos(ftmp(:,iyp))*sin(ftmp(:,izp))&
             +ftmp(:,ivpz)*cos(ftmp(:,izp))                ))
        vzcm=sum(ftmp(:,imass)*(&
             ftmp(:,ivpx)*cos(ftmp(:,iyp))-ftmp(:,ivpy)*sin(ftmp(:,iyp))))
!
        xcm=sum(ftmp(:,imass)*(ftmp(:,ixp)*sin(ftmp(:,iyp))*cos(ftmp(:,izp))))
        ycm=sum(ftmp(:,imass)*(ftmp(:,ixp)*sin(ftmp(:,iyp))*sin(ftmp(:,izp))))
        zcm=sum(ftmp(:,imass)*(ftmp(:,ixp)*cos(ftmp(:,iyp))))
!
        thtcm=atan2(sqrt(xcm**2+ycm**2),zcm)
        phicm=atan2(ycm,xcm)
!
        vcm(1)= vxcm*sin(thtcm)*cos(phicm) + &
            vycm*sin(thtcm)*sin(phicm) + vzcm*cos(thtcm)
        vcm(2)= vxcm*cos(thtcm)*cos(phicm) + &
            vycm*cos(thtcm)*sin(phicm) - vzcm*sin(thtcm)
        vcm(3)=-vxcm*sin(phicm)            + &
            vycm*cos(phicm)
      endif
!
      do k=1,npar_loc
        if (ipar(k)<=mspar) then
          dfp(k,ixp:izp) = dfp(k,ixp:izp) - vcm*totmass1
        endif
      enddo
!
    endsubroutine reset_center_of_mass
!***********************************************************************
    subroutine integrate_selfgravity(p,rrp,xxpar,accg,rp0)
!
!  Calculates acceleration on the point (x,y,z)=xxpar
!  due to the gravity of the gas+dust.
!
!  15-sep-06/wlad : coded
!
      use Mpicomm
!
      real, dimension(nx,3) :: dist
      real, dimension(nx) :: rrp,selfgrav,density
      real, dimension(nx) :: dv,jac,dqy,tmp
      real :: dqx,dqz,rp0,fac
      real, dimension(3) :: xxpar,accg,sum_loc
      integer :: j
      type (pencil_case) :: p
      logical :: lfirstcall=.true.
!
      intent(out) :: accg
!
      if (lfirstcall) &
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
        call fatal_error('integrate_selfgravity', &
            ' not yet implemented for spherical polars')
        dqx=0.;dqy=0.;dqz=0.
      else
        call fatal_error('integrate_selfgravity','wrong coord_system')
        dqx=0.;dqy=0.;dqz=0.
      endif
!
      if (nzgrid==1) then
        dv=dqx*dqy
      else
        dv=dqx*dqy*dqz
      endif
!
!  The gravity of every single cell - should exclude inner and outer radii...
!
!  selfgrav = G*((rho+rhop)*dv)*mass*r*(r**2 + r0**2)**(-1.5)
!  gx = selfgrav * r\hat dot x\hat
!  -> gx = selfgrav * (x-x0)/r = G*((rho+rhop)*dv)*mass*(r**2+r0**2)**(-1.5) * (x-x0)
!
      density=0.
      if (ldust) then !npar>mspar equals dust is being used
        density=density+p%rhop
      else
        density=density+p%rho
      endif
      selfgrav = prhs_cte*density*jac*dv*(rrp**2 + rp0**2)**(-1.5)
!
!  Everything inside the accretion radius of the particle should
!  not exert gravity (numerical problems otherwise)
!
      where (rrp<=rp0)
        selfgrav = 0
      endwhere
!
!  Exclude the frozen zones
!
      if (lexclude_frozen) then
        if (lcylinder_in_a_box) then
          where ((p%rcyl_mn<=r_int).or.(p%rcyl_mn>=r_ext))
            selfgrav = 0
          endwhere
        else
          where ((p%r_mn<=r_int).or.(p%r_mn>=r_ext))
            selfgrav = 0
          endwhere
        endif
      endif
!
!  Integrate the accelerations on this processor
!  And sum over processors with mpireduce
!
      do j=1,3
        tmp=selfgrav*dist(:,j)
        !take proper care of the trapezoidal rule
        !in the case of non-periodic boundaries
        fac = 1.
        if ((m==m1.and.lfirst_proc_y).or.(m==m2.and.llast_proc_y)) then
          if (.not.lperi(2)) fac = .5*fac
        endif
!
        if (lperi(1)) then
          sum_loc(j) = fac*sum(tmp)
        else
          sum_loc(j) = fac*(sum(tmp(2:nx-1))+.5*(tmp(1)+tmp(nx)))
        endif
        call mpireduce_sum(sum_loc(j),accg(j))
      enddo
!
!  Broadcast particle acceleration
!
      call mpibcast_real(accg,3)
!
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine integrate_selfgravity
!***********************************************************************
    subroutine bcast_nbodyarray(fp)
!
!  Broadcast nbody particles across processors
!  The non-root processors, if they find a nbody particle in their
!  fp array, they:
!     send it to root with a true logical
!     else send a false logical
!
!  The root processor receives all the logicals, and if
!  they are true, then receives the value, broadcasting it
!
!  07-sep-06/wlad: coded
!
      use Mpicomm
!
      real, dimension(mpar_loc,mpvar) :: fp
      logical, dimension(mspar) :: lnbody
      integer :: ks,k,tagsend,j,tagrecv
!
      if (lmpicomm) then
!
!  Loop through the nbody particles
!
        do ks=1,mspar
!
!  Set the logical to false initially
!
          lnbody(ks)=.false.
!
!  Loop through the particles on this processor
!
          do k=1,npar_loc
            if (ipar_nbody(ks)==ipar(k)) then
!
!  A nbody was found here. Turn the logical true and copy fp to fsp
!
              lnbody(ks) = .true.
              fsp(ks,ixp:izp)   = fp(k,ixp:izp)
              fsp(ks,ivpx:ivpz) = fp(k,ivpx:ivpz)
              fsp(ks,imass)      = pmass(ks)
!
!  Send it to root. As there will be just mspar calls to
!  mpisend, the tag can be ipar itself
!
              if (.not.lroot) then
                call mpisend_real(fsp(ks,:),mspvar,root,ks)
                if (ip<=6) print*,'logical for particle ',ks,&
                     ' set to true on processor ',iproc, &
                     ' with tag=',ks
              endif
            endif
          enddo
!
!  Send the logicals from each processor. As all processors send mspar calls,
!  the tags are now in base mspar. It assures that two logicals will not have
!  the same tag.
!
          if (.not.lroot) then
            tagsend = mspar*iproc + ks
            call mpisend_logical(lnbody(ks),1,root,tagsend)
          else
!
!  The root receives all logicals. Same tag.
!
            do j=1,ncpus-1
              tagrecv = mspar*j + ks
              call mpirecv_logical(lnbody(ks),1,j,tagrecv)
!
!  Test the received logicals
!
              if (lnbody(ks)) then
!
!  Found a nbody particle. Get the value of fsp
!
                call mpirecv_real(fsp(ks,:),mspvar,j,ks)
                if (ip<=6) print*,'logical for particle ',ks,&
                     ' is true on processor ',j, &
                     ' with tag=',ks,' on root'
              endif
            enddo
          endif
!
!  Broadcast the received fsp
!
          call mpibcast_real(fsp(ks,:),mspvar)
!
!  Print the result in all processors
!
          if (ip<=8)  print*,'bcast_nbodyarray: finished loop. '//&
               'nbody particles in proc ',iproc,&
               ' are fsp(ks,:)=',fsp(ks,:)
        enddo
      else
!
!  Non-parallel. Just copy fp to fsp
!
        do ks=1,mspar
          do k=1,npar_loc
            if (ipar_nbody(ks)==ipar(k)) then
              fsp(ks,ixp:izp)   = fp(k,ixp:izp)
              fsp(ks,ivpx:ivpz) = fp(k,ivpx:ivpz)
            endif
          enddo
          fsp(ks,imass)    = pmass(ks)
        enddo
!
      endif
!
      if (ldebug) print*,'bcast_nbodyarray finished'
!
    endsubroutine bcast_nbodyarray
!***********************************************************************
    subroutine particles_nbody_special
!
!  Fetch fsp array to special module
!
!  01-mar-08/wlad: coded
!
      use Special, only: special_calc_particles_nbody
!
      call special_calc_particles_nbody(fsp)
!
    endsubroutine particles_nbody_special
!***********************************************************************
    subroutine get_totalmass(tmass)
!
!  called from global_shear to set the initial keplerian field
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
!
!  Subroutine that returns the gravity field a particle
!  exterts in a pencil, respecting the coordinate system used
!
!  07-mar-08/wlad: coded
!
      real, dimension(mx),intent(in) :: grr
      real, dimension(mx,3),intent(out) :: gg
      real :: x0,y0,z0
      integer :: ks
!
      x0=fsp(ks,ixp);y0=fsp(ks,iyp);z0=fsp(ks,izp)
!
      if (coord_system=='cartesian') then
        gg(:,1) = (x   -x0)*grr
        gg(:,2) = (y(m)-y0)*grr
        gg(:,3) = (z(n)-z0)*grr
      elseif (coord_system=='cylindric') then
        gg(:,1) = (x   -x0*cos(y(m)-y0))*grr
        gg(:,2) =       x0*sin(y(m)-y0) *grr
        gg(:,3) = (z(n)-z0             )*grr
      elseif (coord_system=='spherical') then
        gg(:,1)  = (x-x0*sinth(m)*sin(y0)*cos(z(n)-z0))*grr
        gg(:,2)  = (x0*(sinth(m)*cos(y0)-&
             costh(m)*sin(y0)*cos(z(n)-z0)))*grr
        gg(:,3)  = (x0*sin(y0)*sin(z(n)-z0))*grr
      endif
!
    endsubroutine get_gravity_field_nbody
!***********************************************************************
    subroutine calc_torque(p,dist,ks)
!
!  Output torque diagnostic for nbody particle ks
!
!  05-nov-05/wlad : coded
!
      use Diagnostics
!
      type (pencil_case) :: p
      real, dimension(nx) :: torque,torqint,torqext
      real, dimension(nx) :: dist,rpre
      real :: rr,w2,smap,hills
      integer :: ks,i
!
      real, dimension(nx) :: rr_mn,dist2,tmp,tempering
      real :: rp,phip,phi,pcut
!
      if (ks==istar) call fatal_error('calc_torque', &
          'Nonsense to calculate torques for the star')
!
      if (lcartesian_coords) then
!
        rr    = sqrt(fsp(ks,ixp)**2 + fsp(ks,iyp)**2 + fsp(ks,izp)**2)
        w2    = fsp(ks,ivpx)**2 + fsp(ks,ivpy)**2 + fsp(ks,ivpz)**2
        smap  = 1./(2./rr - w2)
        hills = smap*(pmass(ks)*pmass1(istar)/3.)**(1./3.)
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
!  Exclude material from inside the Roche Lobe. Compute internal
!  and external torques.
!
          if (dist(i)>=hills) then
            if (p%rcyl_mn(i)>=rr) torqext(i) = torque(i)
            if (p%rcyl_mn(i)<=rr) torqint(i) = torque(i)
          endif
!
        enddo
!
        call sum_lim_mn_name(torqext,idiag_torqext(ks),p)
        call sum_lim_mn_name(torqint,idiag_torqint(ks),p)
!
      else if (lcylindrical_coords) then
        rp = fsp(ks,ixp)
        hills = rp*(pmass(ks)*pmass1(istar)/3.)**(1./3.)
!
        rr_mn = x(l1:l2)     ; phi  = y(m)
        rp    = fsp(ks,ixp)  ; phip = fsp(ks,iyp)
        dist2 = rr_mn**2 + rp**2 - 2*rr_mn*rp*cos(phi-phip)
        rpre  = rr_mn*rp*sin(phi-phip)
        tmp   = -GNewton*pmass(ks)*p%rho*rpre*&
             (dist2 + r_smooth(ks)**2)**(-1.5)*rr_mn
!
        pcut=0.8 !*hills
        tempering = 1./(exp(-(sqrt(dist2)/hills - pcut)/(.1*pcut))+1.)
!
        torque = tmp * tempering
!
        torqext=0.
        torqint=0.
!
        do i=1,nx
!
!  Exclude material from inside the Roche Lobe. Compute internal
!  and external torques.
!
        !  if (dist2(i)>=hills**2) then
            if (p%rcyl_mn(i)>=rp) torqext(i) = torque(i)
            if (p%rcyl_mn(i)<=rp) torqint(i) = torque(i)
        !  endif
!
        enddo
!
        call sum_mn_name(torqext,idiag_torqext(ks))
        call sum_mn_name(torqint,idiag_torqint(ks))
!
      else
        call fatal_error('calc_torque','calc_torque not yet '//&
             'implemented for spherical coordinates.')
      endif
!
    endsubroutine calc_torque
!***********************************************************************
    subroutine get_ramped_mass
!
      real :: ramping_period,tmp
      integer :: ks
!
      ramping_period=2*pi*ramp_orbits
      if (t <= ramping_period) then
        !sin ((pi/2)*(t/(ramp_orbits*2*pi))
        tmp=0.
        !just need to do that for the original masses
        do ks=1,mspar_orig
          if (ks/=istar) then
            pmass(ks)= max(dble(tini),&
                 final_ramped_mass(ks)*(sin((.5*pi)*(t/ramping_period))**2))
            tmp=tmp+pmass(ks)
          endif
        enddo
        pmass(istar)= 1-tmp
      else
        pmass(1:mspar_orig)=final_ramped_mass(1:mspar_orig)
      endif
!
    endsubroutine get_ramped_mass
!***********************************************************************
    subroutine calc_nbodygravity_particles(f)
!
!  For the case of interpolated gravity, add the gravity
!  of all n-body particles to slots of the f-array.
!
!  08-mar-08/wlad: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,3) :: ggt
!
      if (linterpolate_gravity) then
!
        if (lramp) call get_ramped_mass
!
!  Calculate grid - nbody particles distances
!
        do n=1,mz
        do m=1,my
          call get_total_gravity(ggt)
          f(:,m,n,iglobal_ggp:iglobal_ggp+2)=ggt
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
!  Sum the gravities of all massive particles
!
!  08-mar-08/wlad: coded
!
      use Sub
!
      real, dimension (mx,mspar) :: rp_mn,rpcyl_mn
      real, dimension (mx,3)     :: ggp,ggt
      real, dimension (mx)       :: grav_particle,rrp
      integer                    :: ks
!
      intent(out) :: ggt
!
      ggt=0.
      do ks=1,mspar
!
!  Spherical and cylindrical distances
!
        call get_radial_distance(rp_mn(:,ks),rpcyl_mn(:,ks),&
             e1_=fsp(ks,ixp),e2_=fsp(ks,iyp),e3_=fsp(ks,izp))
!
!  Check which particle has cylindrical gravity switched on
!
        if (lcylindrical_gravity_nbody(ks)) then
          rrp = rpcyl_mn(:,ks)
        else
          rrp = rp_mn(:,ks)
        endif
!
!  Gravity field from the particle ks
!
        grav_particle =-GNewton*pmass(ks)*(rrp**2+r_smooth(ks)**2)**(-1.5)
        call get_gravity_field_nbody(grav_particle,ggp,ks)
!
        if ((ks==istar).and.lnogravz_star) &
            ggp(:,3) = 0.
!
!  Sum up the accelerations of the massive particles
!
        ggt=ggt+ggp
!
      enddo
!
    endsubroutine get_total_gravity
!***********************************************************************
    subroutine create_particles_sink_nbody(f,fp,dfp,ineargrid)
!
!  If the flow in any place of the grid has gone gravitationally
!  unstable, substitute the local flow by a sink particle. By now,
!  the only criterion is the local Jeans mass. The Jeans length is
!
!    lambda_J = sqrt(pi*cs2/(G*rho))
!
!  We substitute lambda_J by the biggest resolution element and write
!  the condition in terms of density
!
!     rho_J = pi*cs2/G*Delta_x^2
!
!  The constant is a free parameter of the module.
!
!  13-mar-08/wlad: coded
!
      use EquationOfState, only: cs20
      use FArrayManager,   only: farray_use_global
      use Mpicomm
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mpvar) :: fp, dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      real, dimension(maxsink,mspvar) :: fcsp,fcsp_loc
      real, dimension(maxsink,mspvar) :: fcsp_proc
      real, dimension(nx,3) :: vvpm
      real, dimension(nx) :: rho_jeans,rho_jeans_dust,vpm2
      real, dimension(nx) :: Delta1,dvolume
      integer, dimension(nx,mpar_loc) :: pik
      integer :: i,k,kn,inx0,nc,ks,nc_loc,j
      integer::nc_proc
      integer, pointer :: iglobal_cs2
!
      logical :: ltime_to_create
!
      real, dimension(nx,3)::puu
      real, dimension(nx)  ::prho,prhop,pcs2
      integer, dimension(nx) :: pnp,npik
!
!  Just activate this routine if we want sinks to be created and self-
!  gravity is used. If one does not add the t>=tstart_selfgrav flag, the
!  particles initially have zero velocity dispersion and collapse
!  immediately.
!
      if (lcreate_sinks.and.t>=tstart_selfgrav) then
!
!  just do this for some specific timesteps, otherwise it takes too long!
!
        ltime_to_create = (mod(it-1,icreate) == 0)
        if (ltime_to_create.and.llast) then
!
          if (ldebug) print*,'Entered create_sink_particles_nbody'
!
          do i=1,nx
            if (lcartesian_coords) then
              Delta1(i)=max(dx_1(i),dy_1(mpoint),dz_1(npoint))
            elseif (lcylindrical_coords) then
              Delta1(i)=max(dx_1(i),rcyl_mn1(i)*dy_1(mpoint),dz_1(npoint))
            elseif (lspherical_coords) then
              call fatal_error('create_sink_particle_nbody', &
                  'not yet implemented for spherical polars')
            endif
          enddo
!
!  start number of created particles
!
          fcsp_loc=0.
          nc_loc=0
!
          do imn=1,ny*nz
            n=nn(imn)
            m=mm(imn)
!
!  define the "pencils"
!
            if (ldensity_nolog) then
              prho=f(l1:l2,m,n,irho)
            else
              prho=exp(f(l1:l2,m,n,ilnrho))
            endif
!
            if (ldust)  then
              prhop=f(l1:l2,m,n,irhop)
              pnp  =int(f(l1:l2,m,n,inp))
            endif
!
            if (lhydro) puu=f(l1:l2,m,n,iux:iuz)
!
            if (llocal_iso) then
              call farray_use_global('cs2',iglobal_cs2)
              pcs2=f(l1:l2,m,n,iglobal_cs2)
            else
              pcs2=cs20
            endif
!
            dvolume = dVol_x(l1:l2)*dVol_y(m)*dVol_z(n)
!
!  Jeans analysis of the Gas
!
!  The Jeans length is lambda_J = sqrt(pi*cs2/(G*rho))
!  We substitute lambda_J by the biggest resolution
!  element and write the condition in terms of density
!
!     rho_J = pi*cs2/G*Delta_x^2
!
!  The constant is a free parameter of the module
!
            if (lcreate_gas) then !test purposes
!
              if (ldebug) print*,"Entered create sink from gas"
!
              if (nzgrid/=1) then
                rho_jeans = create_jeans_constant*3.*pi*pcs2*GNewton1*Delta1**2/32
              else
                rho_jeans = create_jeans_constant*   pi*pcs2*GNewton1*Delta1   /8
              endif
!
!  start particle counter for dust particles
!
              do i=1,nx
                if (prho(i)>rho_jeans(i)) then
!
!  create a new particle at the center of the grid
!
                  nc_loc=nc_loc+1
!
!  check if we are not creating too many particles
!
                  if (nc_loc>maxsink) then
                    print*,'too many gas sink particles were created. '//&
                         'Stop and allocated more'
                    print*,'number of created partciles (nc)=',nc_loc
                    print*,'maximum allowed number of sinks '//&
                         'before merging (maxsink)=',maxsink
                    call fatal_error('create_sink_particles_nbody','')
                  endif
!
!  store these particles in a temporary array fcsp - "F array of Created Sink Particles"
!
                  fcsp_loc(nc_loc,ixp) = x(i+nghost)
                  fcsp_loc(nc_loc,iyp) = y(m)
                  fcsp_loc(nc_loc,izp) = z(n)
!
!  give it the velocity of the grid cell
!
                  fcsp_loc(nc_loc,ivpx:ivpz) = puu(i,:)
!
!  the mass of the new particle is the
!  collapsed mass m=[rho - rho_jeans]*dV
!
                  fcsp_loc(nc_loc,imass)   = (prho(i)-rho_jeans(i))*dvolume(i)
!
!  and that amount was lost by the grid, so only rho_jeans remains
!
                  if (ldensity_nolog) then
                    f(i+nghost,m,n,irho)   = rho_jeans(i)
                  else
                    f(i+nghost,m,n,ilnrho) = log(rho_jeans(i))
                  endif
!
                endif
              enddo
            endif!if lcreate_gas
!
!  Jeans analysis of the dust
!
            if (ldust.and.lparticles_selfgravity.and.lcreate_dust) then
!
              if (ldebug) print*,"Entered create sink from dust"
!
!  k1s,k2s and ineargrids are already defined. Substitute sound speed
!  for vpm2=<(vvp-<vvp>)^2>, the particle's velocity dispersion
!
              vvpm=0.0; vpm2=0.0
              do k=k1_imn(imn),k2_imn(imn)
                if (.not.(any(ipar(k)==ipar_nbody))) then
                  inx0=ineargrid(k,1)-nghost
                  vvpm(inx0,:) = vvpm(inx0,:) + fp(k,ivpx:ivpz)
                endif
              enddo
              do i=1,nx
                if (pnp(i)>1.0) vvpm(i,:)=vvpm(i,:)/pnp(i)
              enddo
!  vpm2
              do k=k1_imn(imn),k2_imn(imn)
                if (.not.(any(ipar(k)==ipar_nbody))) then
                  inx0=ineargrid(k,1)-nghost
                  vpm2(inx0) = vpm2(inx0) + (fp(k,ivpx)-vvpm(inx0,1))**2 + &
                       (fp(k,ivpy)-vvpm(inx0,2))**2 + &
                       (fp(k,ivpz)-vvpm(inx0,3))**2
                  endif
              enddo
              do i=1,nx
                if (pnp(i)>1.0) vpm2(i)=vpm2(i)/pnp(i)
              enddo
!
              if (nzgrid/=1) then
                rho_jeans_dust = create_jeans_constant*3*pi*vpm2*GNewton1*Delta1**2/32
              else
                rho_jeans_dust = create_jeans_constant*  pi*vpm2*GNewton1*Delta1   /8
              endif
!
!  Now fill up an array with the identification index of the particles
!  present in each grid cell
!
              npik=0.
              do k=k1_imn(imn),k2_imn(imn)
                if (.not.(any(ipar(k)==ipar_nbody))) then
                  inx0=ineargrid(k,1)-nghost
                  npik(inx0)=npik(inx0)+1
                  pik(inx0,npik(inx0)) = k
                endif
              enddo
!
!  Now check for unstable particle concentrations within this pencil
!
              do i=1,nx
                if (prhop(i)>rho_jeans_dust(i)) then
!
!  This cell is unstable. Remove all particles from it
!
                  if (pnp(i) > 1.0) then
                    !removing must always be done backwards
                    do kn=pnp(i),1,-1
                      if (.not.(any(ipar(k)==ipar_nbody))) &
                           call remove_particle(fp,ipar,pik(i,kn))
                    enddo
!
!  create a new particle at the center of the cell
!
                    nc_loc=nc_loc+1
!
!  check if we are not making too many particles
!
                    if (nc_loc>maxsink) then
                      print*,'too many dust sink particles were created. '//&
                           'Stop and allocated more'
                      print*,'number of created partciles (nc)=',nc_loc
                      print*,'maximum allowed number of sinks '//&
                           'before merging (maxsink)=',maxsink
                      call fatal_error('create_sink_particle_nbody','')
                    endif
!
                    fcsp_loc(nc_loc,ixp) = x(i+nghost)
                    fcsp_loc(nc_loc,iyp) = y(m)
                    fcsp_loc(nc_loc,izp) = z(n)
!
!  give it the mean velocity of the group of particles that
!  was accreted
!
                    fcsp_loc(nc_loc,ivpx:ivpz) = vvpm(i,:)
!
!  the mass of the new particle is the
!  collapsed mass M=m_particle*np
!
                    fcsp_loc(nc_loc,imass)    = mp_swarm*pnp(i)
!
                  endif
                endif
              enddo
!
            endif !if (ldust)
          enddo
          if ((ip<=8).and.nc_loc/=0) &
               print*,'I, processor ',iproc,&
               ' created ',nc_loc,' particles'
!
!  Finished creating them. Now merge them across the processors
!  into a single array
!
          fcsp_proc=0.
          nc=0
          if (lmpicomm) then
            if (.not.lroot) then
              if (ldebug) print*,'processor ',iproc,'sending ',nc_loc,' particles'
              call mpisend_int(nc_loc,1,root,222)
              if (nc_loc/=0) then
                print*,'sending',fcsp_loc(1:nc_loc,:)
                call mpisend_real(fcsp_loc(1:nc_loc,:),(/nc_loc,mspvar/),root,111)
              endif
            else
              !root
              nc_proc=nc_loc
              if (nc_loc/=0) &
                   fcsp(nc+1:nc+nc_proc,:)=fcsp_loc(1:nc_proc,:)
              nc=nc+nc_loc
              !get from the other processors
              do j=1,ncpus-1
                call mpirecv_int(nc_proc,1,j,222)
                if (ldebug) print*,'received ',nc_proc,'particles from processor,',j
!                call fatal_error("","")
                if (nc_proc/=0) then
                  if (ldebug) print*,'receiving ',nc_proc,' particles from processor ',j
                  call mpirecv_real(fcsp_proc(1:nc_proc,:),(/nc_proc,mspvar/),j,111)
                  if (ldebug) print*,'particle received=',fcsp_proc(1:nc_proc,:)
                  fcsp(nc+1:nc+nc_proc,:)=fcsp_proc(1:nc_proc,:)
                  nc=nc+nc_proc
                endif
              enddo
            endif
          else
            !serial, just copy it
            nc =nc_loc
            fcsp=fcsp_loc
          endif
!
!  Call merge only if particles were created
!
          call mpibcast_int(nc,1)
!
          if (nc/=0) then
            call merge_and_share(fcsp,nc,fp)
            call bcast_nbodyarray(fp)
          endif
!
        endif
!
!  print to the screen the positions and velocities of the
!  newly generated sinks
!
        if (ldebug) then
          print*,'---------------------------'
          print*,'the massive particles present are:'
          do ks=1,mspar
            print*,'ks=',ks
            print*,'positions=',fsp(ks,ixp:izp)
            print*,'velocities=',fsp(ks,ivpx:ivpz)
            print*,'mass=',fsp(ks,imass)
            print*,'ipar_nbody=',ipar_nbody
            print*,''
          enddo
          print*,'---------------------------'
        endif
!
      endif
!
      call keep_compiler_quiet(dfp)
!
    endsubroutine create_particles_sink_nbody
!***********************************************************************
    subroutine remove_particles_sink_nbody(f,fp,dfp,ineargrid)
!
!  Subroutine for taking particles out of the simulation due to their
!  proximity to a sink particle or sink point.
!
!  Just a dummy routine for now.
!
!  25-sep-08/anders: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mpvar) :: fp, dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine remove_particles_sink_nbody
!**************************************************************************
    subroutine merge_and_share(fcsp,nc,fp)
!
!  Subroutine to merge a gravitationally unstable
!  cluster of particles into a sink particle
!
!  13-mar-08/wlad: coded
!
      use Mpicomm
!
      real, dimension(mpar_loc,mpvar) :: fp
      real, dimension(maxsink,mspvar) :: fcsp
      real, dimension(nspar,mspvar) :: fleft
!
      integer :: nc,nf,kc,j
!
      real, dimension(0:ncpus-1,nspar,mpvar) :: fcsp_mig
      integer, dimension(0:ncpus-1) :: nsmig
      integer :: iz0,iy0,ipz_rec,ipy_rec,iproc_rec,npar_tot
      integer:: ns,i
      double precision :: dy1,dz1
!
!  How many particles are present in the whole grid? Needed to set
!  the new ipars of the particles that will be created.
!
      call mpiallreduce_sum_int(npar_loc,npar_tot)
!
!  The root processor is the only one that has the right fcsp (the
!  array of particles that were created). Merge them with a friends of
!  friends algorithm
!
      if (lroot) then
!
        call friends_of_friends(fcsp,fleft,nc,nf)
!
!  Friends of friends merged the created particles into few massive
!  clusters. These new clusters are truly sinks. Now we have nf new
!  particles.
!
        pmass(mspar+1:mspar+nf)=fleft(1:nf,imass)
!
!  Allocate the new ipar_nbodies
!
        do i=1,nf
          ipar_nbody(mspar+i)=mspar+i
!  Activate accretion for the newly created sinks
          if (laccrete_when_create) then
            ladd_mass(mspar+i)=.true.
            laccretion(mspar+i)=.true.
          endif
        enddo
        mspar=mspar+nf
!
!  But check if we did not end up with too many particles
!
        if (mspar>nspar) then
          print*,'after friends_of_friends, we still have '//&
               'too many nbody particles.'//&
               'Stop and allocated more'
          print*,'the total number of massive particles (mspar)=',mspar
          print*,'is bigger than the maximum number of '//&
               'allowed massive particles (nspar)=',nspar
          call fatal_error("merge and share","")
        endif
      endif
!
!  Broadcast mspar, ipar_nbody and pmass
!
      call mpibcast_int(mspar,1)
      call mpibcast_int(ipar_nbody(1:mspar),mspar)
      call mpibcast_real(pmass,mspar)
      call mpibcast_logical(ladd_mass,mspar)
      call mpibcast_logical(laccretion,mspar)
!
!  Migrate the particles to their respective processors
!
      if (lmpicomm) then
!
        dy1=1/dy; dz1=1/dz
          !y0_mig=0.5*(y(m1)+y(m1-1)); y1_mig=0.5*(y(m2)+y(m2+1))
          !z0_mig=0.5*(z(n1)+z(n1-1)); z1_mig=0.5*(z(n2)+z(n2+1))
!
!  Only the root processor knows where the particle is (fleft)
!
        nsmig(:)=0
!
        if (lroot) then
!
          do kc=1,nf
!  y-index
            ipy_rec=ipy
            iy0=nint((fleft(kc,iyp)-y(m1))*dy1+nygrid)-nygrid+1
            do while (iy0>ny)
              ipy_rec=ipy_rec+1
              iy0=iy0-ny
            enddo
!  z-index
            ipz_rec=ipz
            iz0=nint((fleft(kc,izp)-z(n1))*dz1+nzgrid)-nzgrid+1
            do while (iz0>nz)
              ipz_rec=ipz_rec+1
              iz0=iz0-nz
            enddo
!  Calculate serial index of receiving processor.
            iproc_rec=ipy_rec+nprocy*ipz_rec
!  Check that the processor exists
            if (iproc_rec>=ncpus .or. iproc_rec<0) then
              call warning('merge_and_share','',iproc)
              print*, 'merge_and_share: receiving proc does not exist'
              print*, 'merge_and_share: iproc, iproc_rec=', &
                   iproc, iproc_rec
              call fatal_error_local("","")
            endif
!
!  Fill up the migration arrays if the particles are not at the root processor
!
            if (iproc_rec/=root) then
              nsmig(iproc_rec)=nsmig(iproc_rec)+1
              fcsp_mig(iproc_rec,nsmig(iproc_rec),:)=fleft(kc,1:mpvar)
            else
!  The particle is at the root processor, so create it here
              npar_loc=npar_loc+1
              fp(npar_loc,:)=fleft(kc,1:mpvar)
              ipar(npar_loc)=npar_tot+1
            endif
          enddo !loop over particles
!
!  Send them now
!
          do j=1,ncpus-1
            ns=nsmig(j)
            call mpisend_int(ns, 1, j, 111)
            if (ns/=0) then
               call mpisend_real(fcsp_mig(j,1:ns,:), (/ns,mpvar/), j, 222)
            endif
          enddo
!
        else !not lroot
!
!  The other processors receive the particles
!
          call mpirecv_int(nsmig(iproc), 1, root, 111)
          ns=nsmig(iproc)
          if (ns/=0) then
            call mpirecv_real(fp(npar_loc+1:npar_loc+ns,:),(/ns,mpvar/),root,222)
!  update the relevant quantities
            do kc=1,ns
              npar_loc=npar_loc+1
              ipar(npar_loc)=npar_tot+kc
            enddo
!
          endif
!
        endif !if lroot
!
      else !serial version. Just copy it all
!
        do kc=1,nf
          npar_loc=npar_loc+1
          fp(npar_loc,:)=fleft(kc,1:mpvar)
          ipar(npar_loc)=npar_tot+1
        enddo
!
      endif
!
      if (ldebug) then
        print*,'merge_and_share finished. '
        print*,'We have now ',mspar,' nbody particles, '
        print*,'located at '
        print*,ipar_nbody(1:mspar)
      endif
!
    endsubroutine merge_and_share
!***********************************************************************
    subroutine friends_of_friends(fcsp,fleft,nc,nfinal)
!
!  Algorithm that cluster particles together based on a proximity
!  threshold. It takes a particle and find all its neighbours (=friends).
!  Then loops through the neighbours finding the neighbours of the
!  neighbours (friends of friends), and so on.
!
!  13-mar-08/wlad: coded
!
      integer, intent(in)  :: nc
      real,    dimension(maxsink,mspvar) :: fcsp
      real,    dimension(nspar,mspvar)   :: fleft
      integer, dimension(nc,nc)          :: clusters
      integer, dimension(nc)             :: cluster
      logical, dimension(nc)             :: lchecked
      real, dimension(mspvar)            :: particle
!
      integer :: ks,kf,nclusters,nf,i,nfinal
      integer :: min_number_members,kpar
!
      intent(out) :: nfinal,fleft
!
      min_number_members=3
!
!  All particles are initially unchecked
!
      lchecked(1:nc)=.false.
!
!  No clusters yet, counter is zero
!
      nclusters=0
!
!  And the final number of particles is the same as the initial
!
      kf=0
!
!  Loop through the particles to cluster them
!
      do ks=1,nc
        if (.not.lchecked(ks)) then
!
!  If a particle is unchecked, start a new cluster
!
          call make_cluster(ks,lchecked,cluster,nc,nf,fcsp)
!
!  Cluster done, check if the number of particles is enough
!
          if (nf >= min_number_members) then
!
!  Okay, it is a cluster, raise the counter
!
            nclusters=nclusters+1
            clusters(1:nf,nclusters)=cluster(1:nf)
!
!  Collapse the cluster into a single sink particle
!
            call collapse_cluster(cluster,nf,fcsp,particle)
!
!  this is all a single particle
!
            kf=kf+1
            fleft(kf,:)=particle
!
          else
            !this (or these) are isolated sinks
            do i=1,nf
              kf=kf+1
              kpar=cluster(i)
              fleft(kf,:)=fcsp(kpar,:)
            enddo
          endif
!
        endif
      enddo
!
      nfinal=kf
!
    endsubroutine friends_of_friends
!***********************************************************************
    subroutine make_cluster(ks,lchecked,cluster,nc,nf,fcsp)
!
!  This subroutine starts a cluster of particles, by checking the
!  first one. Once a particle is checked, it will also look for its
!  neighbours, checking them. In the end, all friends and friends
!  of friends are checked into the cluster.
!
!  13-mar-08/wlad: coded
!
      integer, intent(in) :: nc
      real,    dimension(maxsink,mspvar)   :: fcsp
      integer, intent(inout), dimension(nc) :: cluster
      logical, intent(inout), dimension(nc) :: lchecked
      integer, intent(in)  :: ks
      integer, intent(out) :: nf
!
      integer :: ic
!
!  Start the cluster counter
!
      ic=0
!
!  Tag the particles as checked. The subroutine will also
!  loop through its friends and friends of friends
!
      call check_particle(ks,lchecked,cluster,ic,nc,fcsp)
!
!  When it ends this loop, the cluster is done. Nf is the final
!  number of members in the cluster
!
      nf=ic
!
    endsubroutine make_cluster
!***********************************************************************
    subroutine check_particle(k,lchecked,cluster,ic,nc,fcsp)
!
!  Check (tag) a particle and its neighbours
!
!  13-mar-08/wlad : coded
!
      integer, intent(in) :: nc
      real,    dimension(maxsink,mspvar)   :: fcsp
      integer, intent(inout), dimension(nc) :: cluster
      logical, intent(inout), dimension(nc) :: lchecked
      integer, intent(inout) :: ic
      integer, intent(in) :: k
!
!  Add the particle to the cluster
!
      ic=ic+1
      cluster(ic)=k
!
!  Tag it as visited and add its friends to the cluster
!  as well
!
      lchecked(k)=.true.
      call add_friends(k,lchecked,cluster,ic,nc,fcsp)
!
    endsubroutine check_particle
!***********************************************************************
    subroutine add_friends(par,lchecked,cluster,ic,nc,fcsp)
!
!  Add all neighbours of a particle to the cluster. This subroutine
!  is called by check_particle, but also calls it. The procedure is
!  therefore recursive: the loop will be over when all friends of
!  friends are added to the cluster.
!
!  13-mar-08/wlad: coded
!
      integer, intent(in)                   :: nc
      real,    dimension(maxsink,mspvar)   :: fcsp
      logical, intent(inout), dimension(nc) :: lchecked
      integer, intent(inout), dimension(nc) :: cluster
      integer, intent(inout)                :: ic
      integer, intent(in)                   :: par
!
      real    :: dist,link_length
      integer :: k
!
!  linking length
!
      link_length=4.*dx
!
!  Loop through the particle
!
      do k=1,nc
        if (.not.lchecked(k)) then
          call get_particles_interdistance(&
               fcsp(k,ixp:izp),fcsp(par,ixp:izp),&
               DISTANCE=dist)
!
          if (dist<=link_length) then
!
!  if so, add it to the cluster, tag it and call its friends
!
            call check_particle(k,lchecked,cluster,ic,nc,fcsp)
          endif
!  done
        endif
      enddo
!
    endsubroutine add_friends
!***********************************************************************
    subroutine collapse_cluster(cluster,nf,fcsp,particle)
!
!  Collapse the cluster into a single particle
!  with the total mass, momentum and energy of the
!  collapsed stuff
!
!  13-mar-08/wlad: coded
!
      integer, intent(in) :: nf
      real,    dimension(maxsink,mspvar)   :: fcsp
      integer, intent(inout), dimension(nf) :: cluster
      real, dimension(mspvar) :: particle
      integer :: ic,k
!
      real :: mass,xcm,ycm,vxcm,vycm
      real, dimension(3) :: position,velocity,xxp,vvp
!
      intent(out):: particle
!
      mass=0;position=0;velocity=0
!
      do ic=1,nf
        k=cluster(ic)
        mass=mass+fcsp(k,imass)
        if (lcartesian_coords) then
          xxp(1)=fcsp(k,ixp) ; vvp(1)=fcsp(k,ivpx)
          xxp(2)=fcsp(k,iyp) ; vvp(2)=fcsp(k,ivpy)
          xxp(3)=fcsp(k,izp) ; vvp(3)=fcsp(k,ivpz)
        else if (lcylindrical_coords) then
          xxp(1)=fcsp(k,ixp)*cos(fcsp(k,iyp))
          xxp(2)=fcsp(k,ixp)*sin(fcsp(k,iyp))
          xxp(3)=fcsp(k,izp)
!
          vvp(1)=fcsp(k,ivpx)*cos(fcsp(k,iyp))-fcsp(k,ivpy)*sin(fcsp(k,iyp))
          vvp(2)=fcsp(k,ivpx)*sin(fcsp(k,iyp))+fcsp(k,ivpy)*cos(fcsp(k,iyp))
          vvp(3)=fcsp(k,ivpz)
        else if (lspherical_coords) then
          call fatal_error("","")
        endif
        position=position+fcsp(k,imass)*xxp
        velocity=velocity+fcsp(k,imass)*vvp
      enddo
!
      position=position/mass
      if (lcylindrical_coords) then
        xcm=position(1);vxcm=velocity(1)
        ycm=position(2);vycm=velocity(2)
        position(1)=sqrt(xcm**2+ycm**2)
        position(2)=atan2(ycm,xcm)
!
        velocity(1)= vxcm*cos(position(2))+vycm*sin(position(2))
        velocity(2)=-vxcm*sin(position(2))+vycm*cos(position(2))
      endif
!
!  put it onto the output array
!
      particle(ixp : izp)=position
      particle(ivpx:ivpz)=velocity
      particle(imass)    =mass
!
    endsubroutine collapse_cluster
!***********************************************************************
    subroutine particles_nbody_read_snapshot(filename)
!
!  Read nbody particle info
!
!  01-apr-08/wlad: dummy
!
      use Mpicomm, only:mpibcast_int,mpibcast_real
!
      character (len=*) :: filename
!
      if (lroot) then
        open(1,FILE=filename,FORM='unformatted')
        read(1) mspar
        if (mspar/=0) read(1) fsp(1:mspar,:),ipar_nbody(1:mspar)
        if (ip<=8) print*, 'read snapshot', filename
        close(1)
      endif
!
      call mpibcast_int(mspar,1)
      call mpibcast_real(fsp(1:mspar,:),(/mspar,mspvar/))
      call mpibcast_int(ipar_nbody(1:mspar),mspar)
!
      if (ldebug) then
        print*,'particles_nbody_read_snapshot'
        print*,'mspar=',mspar
        print*,'fsp(1:mspar)=',fsp(1:mspar,:)
        print*,'ipar_nbody(1:mspar)=',ipar_nbody
        print*,''
      endif
!
    endsubroutine particles_nbody_read_snapshot
!***********************************************************************
    subroutine particles_nbody_write_snapshot(snapbase,enum,flist)
!
      use General, only:safe_character_assign
      use Sub, only: update_snaptime, read_snaptime
!
!  Input and output of information about the massive particles
!
!  01-apr-08/wlad: coded
!
      logical, save :: lfirst_call=.true.
      integer, save :: nsnap
      real, save :: tsnap
      character (len=*) :: snapbase,flist
      character (len=fnlen) :: snapname, filename_diag
      logical :: enum,lsnap
      character (len=intlen) :: nsnap_ch
      optional :: flist
!
      if (enum) then
        call safe_character_assign(filename_diag,trim(datadir)//'/tsnap.dat')
        if (lfirst_call) then
          call read_snaptime(filename_diag,tsnap,nsnap,dsnap,t)
          lfirst_call=.false.
        endif
        call update_snaptime(filename_diag,tsnap,nsnap,dsnap,t,lsnap,nsnap_ch)
        if (lsnap) then
          snapname=snapbase//nsnap_ch
!
!  Write number of massive particles and their data
!
          open(lun_output,FILE=snapname,FORM='unformatted')
          write(lun_output) mspar
          if (mspar/=0) write(lun_output) fsp(1:mspar,:),ipar_nbody(1:mspar)
          close(lun_output)
          if (ip<=10 .and. lroot) &
               print*,'written snapshot ', snapname
        endif
      else
!
!  Write snapshot without label
!
        snapname=snapbase
        open(lun_output,FILE=snapname,FORM='unformatted')
        write(lun_output) mspar
        if (mspar/=0) write(lun_output) fsp(1:mspar,:),ipar_nbody(1:mspar)
        close(lun_output)
        if (ip<=10 .and. lroot) &
             print*,'written snapshot ', snapname
      endif
!
      call keep_compiler_quiet(flist)
!
    endsubroutine particles_nbody_write_snapshot
!***********************************************************************
    subroutine particles_nbody_write_spdim(filename)
!
!  Write nspar and mspvar to file.
!
!  01-apr-08/wlad: coded
!
      character (len=*) :: filename
!
      open(1,file=filename)
      write(1,'(3i9)') nspar, mspvar
      close(1)
!
    endsubroutine particles_nbody_write_spdim
!***********************************************************************
    subroutine rprint_particles_nbody(lreset,lwrite)
!
!  Read and register print parameters relevant for nbody particles.
!
!  17-nov-05/anders+wlad: adapted
!
      use Diagnostics
      use General, only: itoa
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      integer :: iname,ks,j
      character :: str
      character (len=intlen) :: sks
!
!  Write information to index.pro
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
      if (lwr) then
        write(3,*) 'imass=', imass
      endif
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
        sks=itoa(ks)
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
!  Diagnostic related to quantities summed over all particles
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),&
             'totenergy',idiag_totenergy)
        call parse_name(iname,cname(iname),cform(iname),&
             'totangmom',idiag_totangmom)
      enddo
!
      if (lwr) then
!
      endif
!
    endsubroutine rprint_particles_nbody
!***********************************************************************
 endmodule Particles_nbody
