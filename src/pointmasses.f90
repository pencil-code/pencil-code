! $Id$
!
!  This module takes care of direct N-body gravity between point masses.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MQVAR CONTRIBUTION 7
! MQAUX CONTRIBUTION 0
! CPARAM logical, parameter :: lpointmasses=.true.
!
!***************************************************************
module PointMasses
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Cparam
!
  implicit none
!
  include 'pointmasses.h'
!
  character (len=10), dimension(mqarray) :: qvarname
  real, dimension(nqpar,mqarray) :: fq
  real, dimension(nqpar,mqvar) :: dfq
  real, dimension(nqpar) :: xq0=0.0, yq0=0.0, zq0=0.0
  real, dimension(nqpar) :: vxq0=0.0, vyq0=0.0, vzq0=0.0
  real, dimension(nqpar) :: pmass=0.0, r_smooth=impossible, pmass1
  real, dimension(nqpar) :: frac_smooth=0.3
  real, dimension(nqpar) :: accrete_hills_frac=0.2, final_ramped_mass=0.0
  real :: eccentricity=0.0, semimajor_axis=1.0
  real :: totmass, totmass1
  real :: GNewton1, GNewton=impossible, density_scale=0.001
  real :: cdtq=0.1
  real :: hills_tempering_fraction=0.8
  real, pointer :: rhs_poisson_const, tstart_selfgrav
  integer :: ramp_orbits=5
  integer :: iglobal_ggp=0, istar=1, iplanet=2
!
  integer :: imass=0, ixq=0, iyq=0, izq=0, ivxq=0, ivyq=0, ivzq=0
  integer :: nqvar=0, nqaux=0
!
  logical, dimension(nqpar) :: lcylindrical_gravity_nbody=.false.
  logical, dimension(nqpar) :: lfollow_particle=.false., laccretion=.false.
  logical, dimension(nqpar) :: ladd_mass=.false.
  logical :: lcalc_orbit=.true., lbackreaction=.false., lnorm=.true.
  logical :: lreset_cm=.false., lnogravz_star=.false., lexclude_frozen=.false.
  logical :: lnoselfgrav_star=.true.
  logical :: lramp=.false.
  logical :: ldt_pointmasses=.true.
  logical :: linterpolate_gravity=.false., linterpolate_linear=.true.
  logical :: linterpolate_quadratic_spline=.false.
  logical :: ldust=.false.
  logical :: ltempering=.false.
  logical :: lretrograde=.false.
  logical :: lgas_gravity=.false.,ldust_gravity=.false.
  logical :: linertial_frame=.true.
!
  character (len=labellen) :: initxxq='random', initvvq='nothing'
  character (len=2*bclen+1) :: bcqx='p', bcqy='p', bcqz='p'
!
  type IndexDustParticles
    integer :: ixw=0,izw=0,ivxw=0,ivzw=0
  endtype IndexDustParticles
  type (IndexDustParticles) :: index
!
  namelist /pointmasses_init_pars/ &
      initxxq, initvvq, xq0, yq0, zq0, vxq0, vyq0, vzq0,  &
      pmass, r_smooth, lcylindrical_gravity_nbody, lexclude_frozen, GNewton, &
      bcqx, bcqy, bcqz, ramp_orbits, lramp, final_ramped_mass, density_scale, &
      linterpolate_gravity, linterpolate_quadratic_spline, laccretion, &
      accrete_hills_frac, istar, &
      ladd_mass, ldt_pointmasses, cdtq, lretrograde, &
      linertial_frame, eccentricity, semimajor_axis
!
  namelist /pointmasses_run_pars/ &
      lcalc_orbit, lreset_cm, &
      lnogravz_star, lfollow_particle, lbackreaction, lexclude_frozen, &
      GNewton, bcqx, bcqy, bcqz, density_scale, lnoselfgrav_star, &
      linterpolate_quadratic_spline, laccretion, accrete_hills_frac, istar, &
      ladd_mass, ldt_pointmasses, cdtq, hills_tempering_fraction, &
      ltempering, lgas_gravity, ldust_gravity, linertial_frame
!
  integer, dimension(nqpar,3) :: idiag_xxq=0,idiag_vvq=0
  integer, dimension(nqpar)   :: idiag_torqint=0,idiag_torqext=0
  integer, dimension(nqpar)   :: idiag_torqext_gas=0,idiag_torqext_par=0
  integer, dimension(nqpar)   :: idiag_torqint_gas=0,idiag_torqint_par=0
  integer                     :: idiag_totenergy=0
!
  contains
!***********************************************************************
    subroutine register_pointmasses()
!
!  Set up indices for access to the f and fq.
!
!  27-aug-06/wlad: adapted
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  No need to solve the N-body equations for non-N-body problems.
!
      if (nqpar < 2) then
        if (lroot) write(0,*) 'nqpar = ', nqpar
        call fatal_error('register_pointmasses','the number of massive'//&
             ' particles is less than 2. There is no need to use the'//&
             ' N-body code. Consider setting gravity as a global variable'//&
             ' using gravity_r.f90 instead.')
      endif
!
      !if (npar < nqpar) then
      !  if (lroot) write(0,*) 'npar, nqpar = ', npar, nqpar
      !  call fatal_error('register_pointmasses','the number of massive'//&
      !       ' particles (nqpar) is less than the allocated number of particles'//&
      !       ' (npar). Increase npar to the minimum number (npar=npsar) needed'//&
      !       ' in cparam.local and recompile')
      !endif
!
!  Auxiliary variables for polar coordinates
!
      ixq = nqvar+1
      qvarname(nqvar+1)='ixq'
      iyq = nqvar+2
      qvarname(nqvar+1)='iyq'
      izq = nqvar+3
      qvarname(nqvar+1)='izq'
      nqvar=nqvar+3
!
      ivxq = nqvar+1
      qvarname(nqvar+1)='ivxq'
      ivyq = nqvar+2
      qvarname(nqvar+1)='ivyq'
      ivzq = nqvar+3
      qvarname(nqvar+1)='ivzq'
      nqvar=nqvar+3
!
!  Set up mass as particle index. Plus seven, since the other 6 are
!  used by positions and velocities.
!
      imass=nqvar+1
      qvarname(nqvar+1)='imass'
      nqvar=nqvar+1

      !ivxq_cart = mqvar+1
      !qvarname(mqvar+1)='ivxq_cart'
      !ivyq_cart = mqvar+2
      !qvarname(mqvar+1)='ivyq_cart'
      !ivzq_cart = mqvar+3
      !qvarname(mqvar+1)='ivzq_cart'
      !mqvar=mqvar+3
!
!  Check that the fq and dfq arrays are big enough.
!
      if (nqvar > mqvar) then
        if (lroot) write(0,*) 'nqvar = ', nqvar, ', mqvar = ', mqvar
        call fatal_error('register_pointmasses','nqvar > mqvar')
      endif
!
    endsubroutine register_pointmasses
!***********************************************************************
    subroutine initialize_pointmasses(f)
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
!
      integer :: ks
!
!  Look for initialized masses.
!
      !nqpar_orig=0
      do ks=1,nqpar
        if (pmass(ks)/=0.0) then
          fq(ks,imass)=pmass(ks)
        else
          call fatal_error("initialize_pointmasses",&
                "one of the bodies has zero mass")  
        endif
      enddo
!
!  Just needs to do this starting, otherwise mqar will be overwritten after
!  reading the snapshot
!
      if (.not.lstart) pmass(1:nqpar)=fq(1:nqpar,imass)

!
!  When first called, nqpar was zero, so no diagnostic index was written to
!  index.pro
!
      if (lroot) open(3, file=trim(datadir)//'/index.pro', &
          STATUS='old', POSITION='append')
      call rprint_pointmasses(.false.,LWRITE=lroot)
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
          call get_shared_variable('rhs_poisson_const',rhs_poisson_const,caller='initialize_pointmasses')
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
      if (lselfgravity) &
        call get_shared_variable('tstart_selfgrav',tstart_selfgrav,caller='initialize_pointmasses')
!
!  inverse mass
!
      if (lstart) then
        if (lramp) then
          do ks=1,nqpar
            if (ks/=istar) pmass(ks) = epsi
          enddo
          pmass(istar)=1-epsi*(nqpar-1)
        endif
      else
        !read imass from the snapshot
        pmass(1:nqpar)=fq(1:nqpar,imass)
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
!  from cdata and the one from point mass.
!
      if (((lcylindrical_gravity).and.&
        (.not.lcylindrical_gravity_nbody(istar))).or.&
             (.not.lcylindrical_gravity).and.&
             (lcylindrical_gravity_nbody(istar))) then
        call fatal_error('initialize_pointmasses','inconsitency '//&
            'between lcylindrical_gravity from cdata and the '//&
            'one from point mass')
      endif
!
!  Smoothing radius 
!
      if (any(r_smooth == impossible)) then 
        do ks=1,nqpar
          if (ks/=istar) then
            r_smooth(ks) = frac_smooth(ks) * xq0(ks) * (pmass(ks)/3.)**(1./3)
          else
            r_smooth(ks) = rsmooth
          endif
        enddo
      endif
!
      if (rsmooth/=r_smooth(istar)) then
        print*,'rsmooth from cdata=',rsmooth
        print*,'r_smooth(istar)=',r_smooth(istar)
        call fatal_error('initialize_pointmasses','inconsistency '//&
            'between rsmooth from cdata and the '//&
            'one from point mass')
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
        call fatal_error('initialize_pointmasses','')
      endif
!
!  The presence of dust particles needs to be known.
!
      if (npar > nqpar) ldust=.true.
      if (linterpolate_gravity.and.(.not.ldust)) then
        if (lroot) print*,'interpolate gravity is just'//&
            ' for the dust component. No need for it if'//&
            ' you are not using dust particles'
        call fatal_error('initialize_pointmasses','')
      endif
      if (any(laccretion).and.linterpolate_gravity) then
        if (lroot) print*,'interpolate gravity  not yet '//&
            'implemented in connection with accretion'
        call fatal_error('initialize_pointmasses','')
      endif
!
      if (linterpolate_gravity) then
         if (lroot) print*,'initializing global array for point mass gravity'
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
        call fatal_error('initialize_pointmasses','')
      endif
!
      if ((.not.linertial_frame).and.((.not.lcartesian_coords).or.nqpar>2))  &
           call fatal_error('initialize_pointmasses','Fixed star only for Cartesian and nqpar=2')
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_pointmasses
!***********************************************************************
    subroutine pencil_criteria_pointmasses()
!
!  All pencils that the pointmasses module depends on are specified here.
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
      !if (any(idiag_torqext/=0) .or. any(idiag_torqint/=0)) then
      !  lpenc_diagnos(i_rcyl_mn)=.true.
      !endif
!
      lpenc_diagnos(i_rcyl_mn)=.true.
      lpenc_diagnos(i_rhop)=.true.
!
    endsubroutine pencil_criteria_pointmasses
!***********************************************************************
    subroutine pencil_interdep_pointmasses(lpencil_in)
!
!  Interdependency among pencils provided by the pointmasses module
!  is specified here.
!
!  22-sep-06/wlad: adapted
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_pointmasses
!***********************************************************************
    subroutine calc_pencils_pointmasses(f,p)
!
!  Calculate point mass particle pencils
!
!  22-sep-06/wlad: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_pointmasses
!***********************************************************************
    subroutine init_pointmasses(f)!,fp)
!
!  Initial positions and velocities of point mass particles.
!  Overwrite the position asserted by the dust module
!
!  17-nov-05/anders+wlad: adapted
!
      use General, only: random_number_wrapper
      use Sub
      use Mpicomm, only: mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f
      !real, dimension (mpar_loc,mparray) :: fp
      real, dimension(nqpar) :: kep_vel,sma
      real, dimension(nqpar,3) :: velocity
      real, dimension(nqpar,3) :: positions
      real :: tmp,parc
      integer :: k,ks
!
      intent (in) :: f
      !intent (out) :: fp
!
!  Shortcuts
!
      positions(1:nqpar,1) = xq0 ; velocity(1:nqpar,1) = vxq0
      positions(1:nqpar,2) = yq0 ; velocity(1:nqpar,2) = vyq0
      positions(1:nqpar,3) = zq0 ; velocity(1:nqpar,3) = vzq0
!
!  Initialize particles' positions.
!
      select case (initxxq)
!
      case ('nothing')
        if (lroot) print*, 'init_pointmasses: nothing'
!
      case ('origin')
        if (lroot) then
          print*, 'init_pointmasses: All point mass particles at origin'
          fq(1:nqpar,ixq:izq)=0.0
        endif
!
      case ('constant')
        if (lroot) &
            print*, 'init_pointmasses: Place point mass particles at x,y,z=', &
            xq0, yq0, zq0
        do k=1,nqpar
           fq(k,ixq:izq)=positions(k,1:3)
        enddo
!
      case ('random')
        if (lroot) print*, 'init_pointmasses: Random particle positions'
        do ks=1,nqpar
          if (nxgrid/=1) call random_number_wrapper(positions(ks,ixq))
          if (nygrid/=1) call random_number_wrapper(positions(ks,iyq))
          if (nzgrid/=1) call random_number_wrapper(positions(ks,izq))
        enddo
!
        if (nxgrid/=1) &
             positions(1:nqpar,ixq)=xyz0_loc(1)+positions(1:nqpar,ixq)*Lxyz_loc(1)
        if (nygrid/=1) &
             positions(1:nqpar,iyq)=xyz0_loc(2)+positions(1:nqpar,iyq)*Lxyz_loc(2)
        if (nzgrid/=1) &
             positions(1:nqpar,izq)=xyz0_loc(3)+positions(1:nqpar,izq)*Lxyz_loc(3)
!
        do k=1,nqpar
!
            fq(k,ixq:izq)=positions(k,1:3)
!
!  Correct for non-existing dimensions (not really needed, I think).
!
            if (nxgrid==1) fq(k,ixq)=x(nghost+1)
            if (nygrid==1) fq(k,iyq)=y(nghost+1)
            if (nzgrid==1) fq(k,izq)=z(nghost+1)
!
        enddo
!
      case ('fixed-cm')
        if (lgrav) then
          print*,'a gravity module is being used. Are you using '//&
               'both a fixed central gravity and point mass gravity? '//&
               'better stop and check'
          call fatal_error('init_pointmasses','')
        endif
!
!  Ok, I have the masses and the positions of all massive particles
!  except the last, which will have a position determined to fix the
!  center of mass on the center of the grid.
!
        if (any(yq0/=0)) then
          if (lspherical_coords) then
            call fatal_error('init_pointmasses','not yet generalized'//&
                 ' for non-zero initial inclinations')
          else
            call fatal_error('init_pointmasses','not yet generalized'//&
                 ' for non-zero azimuthal initial position')
          endif
        endif
        if (any(zq0/=0)) then
          if (lspherical_coords) then
            call fatal_error('init_pointmasses','not yet generalized'//&
                 ' for non-zero azimuthal initial position')
          else
            call fatal_error('init_pointmasses','point mass code not'//&
                 ' yet generalized to allow initial inclinations')
          endif
        endif
        if (lcylindrical_coords.or.lspherical_coords) then
          if (any(xq0<0)) &
              call fatal_error('init_pointmasses', &
              'in cylindrical coordinates '//&
              'all the radial positions must be positive')
        endif
!
        if (lroot) then
          print*,'init_pointmasses: fixed-cm - mass and position arranged'
          print*,'                      so that the center of mass is at rest'
          print*,'                      at the center of the grid.'
        endif
!
        if (lspherical_coords) then
          if (lroot) print*,'put all particles in the midplane'
          positions(1:nqpar,iyq)=pi/2
        endif
!
        tmp = 0.;parc=0
        do ks=1,nqpar
          if (ks/=istar) then
            sma(ks)=abs(positions(ks,1))
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
            call fatal_error('init_pointmasses', &
            'The mass of one '//&
            '(or more) of the particles is too big! The masses should '//&
            'never be bigger than g0. Please scale your assemble so that '//&
            'the combined mass of the (n-1) particles is less than that. '//&
            'The mass of the last particle in the pmass array will be '//&
            'reassigned to ensure that the total mass is g0')
!
        do ks=1,nqpar
          if (ks/=istar) &
              positions(ks,1)=sign(1.,positions(ks,1))* (sma(ks) + parc)
        enddo
!
!  The last one (star) fixes the CM at Rcm=zero
!
        if (lcartesian_coords) then
          positions(istar,1)=parc
        elseif (lcylindrical_coords) then
          !put the star in positive coordinates, with pi for azimuth
          positions(istar,1)=abs(parc)
          positions(istar,2)=pi
        elseif (lspherical_coords) then
          positions(istar,1)=abs(parc)
          positions(istar,3)=pi
        endif
!
        if (ldebug) then
          print*,'pmass =',pmass
          print*,'position (x)=',positions(:,1)
          print*,'position (y)=',positions(:,2)
          print*,'position (z)=',positions(:,3)
        endif
!
        do k=1,nqpar
!
!  Here I substitute the first nqpar dust particles by massive ones,
!  since the first ipars are less than nqpar
!
            fq(k,ixq:izq)=positions(k,1:3)
!
!  Correct for non-existing dimensions (not really needed, I think)
!
            if (nxgrid==1) fq(k,ixq)=x(nghost+1)
            if (nygrid==1) fq(k,iyq)=y(nghost+1)
            if (nzgrid==1) fq(k,izq)=z(nghost+1)
!
        enddo
!
      case ('eccentric')
!
!  Coded only for 2 bodies
!
        if (nqpar /= 2) call fatal_error("init_pointmasses",&
             "This initial condition is currently coded for 2 massive particles only.")
!
!  Define iplanet. istar=1 and iplanet=2 is default
!
        if (istar == 2) iplanet=1
!
!  Radial position at barycentric coordinates. Start both at apocenter,
!
!     r_i=(1+e)*a_i, where a_i = sma * m_j /(mi+mj)
!
!  See, i.e., Murray & Dermott, p.45, barycentric orbits.
!
        positions(iplanet,1)=(1+eccentricity) * semimajor_axis * pmass(  istar)/totmass
        positions(  istar,1)=(1+eccentricity) * semimajor_axis * pmass(iplanet)/totmass
!
!  Azimuthal position. Planet and star phased by pi.
!
        positions(iplanet,2)=0
        positions(  istar,2)=pi
!
        do k=1,nqpar
          !if (ipar(k) <= nqpar) then
            fq(k,ixq:izq) = positions(k,1:3)
          !endif
        enddo
!
      case default
        if (lroot) print*,'init_pointmasses: No such such value for'//&
            ' initxxq: ',trim(initxxq)
        call fatal_error('init_pointmasses','')
!
      endselect
!
!  Redistribute particles among processors (now that positions are determined).
!
      call boundconds_pointmasses
      if (lroot) then
        print*,'init_pointmasses: particle positions'
        do ks=1,nqpar
          print*,'particle',ks, 'position x,y,z:',fq(ks,ixq:izq)
        enddo
      endif
!
!  Initial particle velocity.
!
      select case (initvvq)
!
      case ('nothing')
        if (lroot) print*, 'init_pointmasses: No particle velocity set'
!
      case ('zero')
        if (lroot) then
          print*, 'init_pointmasses: Zero particle velocity'
          fq(1:nqpar,ivxq:ivzq)=0.
        endif
!
      case ('constant')
        if (lroot) then
          print*, 'init_pointmasses: Constant particle velocity'
          print*, 'init_pointmasses: vxq0, vyq0, vzq0=', vxq0, vyq0, vzq0
        endif
        do k=1,nqpar
           !if (ipar(k) <= nqpar) &
              fq(k,ivxq:ivzq)=velocity(k,1:3)
        enddo
!
      case ('fixed-cm')
!
!  Keplerian velocities for the planets
!
        parc=0.
        do ks=1,nqpar
          if (ks/=istar) then
            kep_vel(ks)=sqrt(1./sma(ks)) !circular velocity
            parc = parc - kep_vel(ks)*pmass(ks)
          endif
        enddo
        parc = parc*totmass
        do ks=1,nqpar
          if (ks/=istar) then
            if (lcartesian_coords) then
              velocity(ks,2) = sign(1.,positions(ks,1))*(kep_vel(ks) + parc)
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
!  Revert all velocities if retrograde
!
        if (lretrograde) velocity=-velocity
!
!  Allocate the point mass particles
!
         do k=1,nqpar
!
             fq(k,ivxq:ivzq) = velocity(k,1:3)
!
         enddo
!
      case ('eccentric')
!
!  Coded only for 2 bodies
!
        if (nqpar /= 2) call fatal_error("init_pointmasses",&
             "This initial condition is currently coded for 2 massive particles only.")
!
!  Define iplanet. istar=1 and iplanet=2 is default.
!
        if (istar == 2) iplanet=1
        velocity(iplanet,2) = sqrt((1-eccentricity)/(1+eccentricity) * GNewton/semimajor_axis) * pmass(  istar)/totmass
        velocity(  istar,2) = sqrt((1-eccentricity)/(1+eccentricity) * GNewton/semimajor_axis) * pmass(iplanet)/totmass
!
!  Revert all velocities if retrograde.
!
        if (lretrograde) velocity=-velocity
!
!  Loop through particles to allocate the velocities.
!
        do k=1,nqpar
          fq(k,ivxq:ivzq) = velocity(k,1:3)
        enddo
!
      case default
        if (lroot) print*, 'init_pointmasses: No such such value for'//&
             'initvvq: ',trim(initvvq)
        call fatal_error('init_pointmasses','')
!
      endselect
!
!  Make the particles known to all processors
!
      call boundconds_pointmasses()!fq)!,ipar)
      call mpibcast_real(fq,(/nqpar,mqarray/))
!
      if (lroot) then
        print*,'init_pointmasses: particle velocities'
        do ks=1,nqpar
          print*,'particle',ks, 'velocities x,y,z:',fq(ks,ivxq:ivzq)
        enddo
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_pointmasses
!***********************************************************************
    subroutine pointmasses_pde_pencil(f,df,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      !call dxxpq_dt_pointmasses_pencil(f,df,fp,dfp,p,ineargrid)
      call dvvq_dt_pointmasses_pencil(f,df,p)
!
    endsubroutine pointmasses_pde_pencil
!***********************************************************************         
    subroutine dvvq_dt_pointmasses_pencil(f,df,p)!,fp,dfp,p,ineargrid)
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
      type (pencil_case) :: p
!
      real, dimension (nx,nqpar) :: rp_mn, rpcyl_mn
      real, dimension (mx,3) :: ggt
      real, dimension (nx) :: pot_energy
      real, dimension (3) :: xxpar,accg
      integer :: ks, k
      logical :: lintegrate, lparticle_out
!
      intent (in) :: f, p
      intent (inout) :: df
!
!  Get the total gravity field. In the case of dust, it is already
!  pre-calculated
!
      if (lhydro) then
        call get_total_gravity(ggt)
        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + ggt(l1:l2,:)
!
!  Diagnostic
!
        if (ldiagnos) then
          do ks=1,nqpar
            if (idiag_totenergy/=0) pot_energy=0.0
            call get_radial_distance(rp_mn(:,ks),rpcyl_mn(:,ks),&
                 E1_=fq(ks,ixq),E2_=fq(ks,iyq),E3_=fq(ks,izq))
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
                   GNewton*pmass(ks)*(rpcyl_mn(:,ks)**2+r_smooth(ks)**2)**(-0.5)
              if (ks==nqpar) &
                   call sum_lim_mn_name(.5*p%rho*p%u2 + pot_energy,idiag_totenergy,p)
            endif
          enddo
        endif
      endif
!
    endsubroutine dvvq_dt_pointmasses_pencil
!***********************************************************************
    subroutine pointmasses_pde(f,df)
!
      use Mpicomm, only: mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
!
      real, dimension(nqpar) :: hill_radius_square
!
      if (lroot) then
        call calc_hill_radius(hill_radius_square)
!        
        call dxxq_dt_pointmasses
        call dvvq_dt_pointmasses(hill_radius_square)
!
        call correct_curvilinear
      endif
!
      if (lparticles) then
        call mpibcast_real(hill_radius_square,nqpar)
        call dvvp_dt_dustparticles(hill_radius_square)
      endif
!
      call keep_compiler_quiet(f,df)
!    
    endsubroutine  pointmasses_pde
!***********************************************************************
    subroutine calc_hill_radius(hill_radius_square)

      real, dimension(nqpar), intent(out) :: hill_radius_square
      integer :: ks
      real :: rr,w2,sma2
!
      do ks=1,nqpar
        if (laccretion(ks).and.(ks/=istar)) then
          if (lcartesian_coords) then
            rr=sqrt(fq(ks,ixq)**2 + fq(ks,iyq)**2 + fq(ks,izq)**2)
          elseif (lcylindrical_coords) then
            rr=fq(ks,ixq)
            if (nzgrid/=1) rr=sqrt(fq(ks,ixq)**2+fq(ks,izq)**2)
          elseif (lspherical_coords) then
            rr=fq(ks,ixq)
          else
            call fatal_error('calc_hill_radius','wrong coord system')
            rr=0.0
          endif
!
! particle velocities are non-coordinate (linear)
!
          w2    = fq(ks,ivxq)**2 + fq(ks,ivyq)**2 + fq(ks,ivzq)**2
!
! squared semi major axis - assumes GM=1, so beware...
!
          if (totmass/=1) call fatal_error('calc_hill_radius',&
               'can only calculate semimajor axis for normalized total mass')
          sma2  = (rr/(2-rr*w2))**2
!
! squared hills radius
!
          hill_radius_square(ks)=sma2*(pmass(ks)*pmass1(istar)/3)**(2./3.)
        else
          hill_radius_square(ks)=0.0
        endif
      enddo

      
    endsubroutine calc_hill_radius
!***********************************************************************
    subroutine dxxq_dt_pointmasses
!
!  If the center of mass of the point masses was moved from the
!  center of the grid, reset it.
!
!  22-sep-06/wlad: coded
!
      if (lcartesian_coords) then
!
        if (nxgrid/=1) dfq(:,ixq) = dfq(:,ixq) + fq(:,ivxq)
        if (nygrid/=1) dfq(:,iyq) = dfq(:,iyq) + fq(:,ivyq)
        if (nzgrid/=1) dfq(:,izq) = dfq(:,izq) + fq(:,ivzq)
!
      else if (lcylindrical_coords) then
!
        if (nxgrid/=1) dfq(:,ixq) = dfq(:,ixq) + fq(:,ivxq)
        if (nygrid/=1) dfq(:,iyq) = dfq(:,iyq) + fq(:,ivyq)/max(fq(:,ixq),tini)
        if (nzgrid/=1) dfq(:,izq) = dfq(:,izq) + fq(:,ivzq)

      else if (lspherical_coords) then
!
        if (nxgrid/=1) dfq(:,ixq) = dfq(:,ixq) + fq(:,ivxq)
        if (nygrid/=1) dfq(:,iyq) = dfq(:,iyq) + fq(:,ivyq)/ max(fq(:,ixq),tini)
        if (nzgrid/=1) dfq(:,izq) = dfq(:,izq) + fq(:,ivzq)/(max(fq(:,ixq),tini)*sin(fq(:,iyq)))
!
      endif
!
      if (lreset_cm) call reset_center_of_mass
!
    endsubroutine dxxq_dt_pointmasses
!***********************************************************************
    subroutine dvvq_dt_pointmasses(hill_radius_square)
!
!  Evolution of point masses and dust particles velocities due to particle-particle
!  interaction only.
!
!  Coriolis and shear are already added in particles_dust.
!
!  27-aug-06/wlad: coded
!
      use Sub
!
      real, dimension(nqpar) :: hill_radius_square
      integer :: k, ks, j, jpos, jvel
      logical :: lheader, lfirstcall=.true.
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
!
!  Identify module.
!
      if (lheader) print*, 'dvvp_dt_pointmasses: Calculate dvvp_dt_pointmasses'
!
!  Evolve massive particle positions due to the gravity of other massive
!  particles. The gravity of the massive particles on the dust will be added
!  inside the pencil in dvvp_dt_pointmasses_pencil.
!
      if (lramp) call get_ramped_mass
!
!  To the N^2 problem of adding gravity of each massive particle to another. 
!
      do k=1,nqpar
        call loop_through_pointmasses(k,hill_radius_square,.true.)
      enddo
!
!  Position and velocity diagnostics (per massive particle).
!
      if (ldiagnos) then
        do ks=1,nqpar
          if (lfollow_particle(ks)) then
            do j=1,3
              jpos=j+ixq-1 ; jvel=j+ivxq-1
              if (idiag_xxq(ks,j)/=0) then
                call point_par_name(fq(ks,jpos),idiag_xxq(ks,j))
              endif
              if (idiag_vvq(ks,j)/=0) &
                  call point_par_name(fq(ks,jvel),idiag_vvq(ks,j))
            enddo
          endif
        enddo
      endif
!
      if (lheader) print*,'dvvp_dt_pointmasses: Finished dvvp_dt_pointmasses'
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine dvvq_dt_pointmasses
!************************************************************
    subroutine dvvp_dt_dustparticles(hill_radius_square)

      use Particles_main, only: fetch_fp_array,return_fp_array

      real, dimension (mpar_loc,mparray) :: fp_aux
      real, dimension (mpar_loc,mpvar) :: dfp_aux
      integer :: np_aux,k
      logical, dimension(mpar_loc) :: flag=.false.
      real, dimension(nqpar) :: hill_radius_square
!
      call fetch_fp_array(fp_aux,dfp_aux,np_aux,&
           index%ixw,index%izw,index%ivxw,index%ivzw)
      if (np_aux/=0) then
        do k=1,np_aux
          call loop_through_pointmasses(k,hill_radius_square,.false.,&
               fp_aux(k,:),dfp_aux(k,:),flag(k))
        enddo
      endif
      call return_fp_array(fp_aux,dfp_aux,flag)

    endsubroutine dvvp_dt_dustparticles
!************************************************************
    subroutine loop_through_pointmasses(k,hill_radius_square,&
         lcallpointmass,fp_pt,dfp_pt,flag_pt)
!      
      real, dimension(nqpar) :: hill_radius_square
      logical, intent(in) :: lcallpointmass
      integer, intent(in) :: k
!
      real, dimension (mparray), optional :: fp_pt
      real, dimension (mpvar), optional :: dfp_pt
      logical, optional :: flag_pt
!
      !
      !if (linertial_frame) then
      if (lcallpointmass) then 
        call loop_through_pointmasses_inertial(k,hill_radius_square,lcallpointmass)
      else
        call loop_through_pointmasses_inertial(k,hill_radius_square,lcallpointmass,&
             fp_pt,dfp_pt,flag_pt)
      endif
      !else
      !  call loop_through_nbodies_fixstar(fp,dfp,k,sq_hills,ineargrid)
      !endif
!
      endsubroutine loop_through_pointmasses
!************************************************************
      subroutine loop_through_pointmasses_inertial(k,hill_radius_square,&
           lcallpointmass,fp_pt,dfp_pt,flag_pt)
!
!  Subroutine that adds the gravity from all massive particles.
!  Particle gravity has always to be added in Cartesian, for better
!  conservation of the Jacobi constant.
!
!  07-mar-08/wlad:coded
!
      integer, intent(in) :: k
      real, dimension (nqpar) :: hill_radius_square
!
      real :: r2_ij, rs2, invr3_ij, v_ij,tmp1,tmp2
      integer :: ks,ip1,ip3,iv1,iv3
!
      real, dimension (3) :: evr_cart,acc_cart
!
      logical, intent(in) :: lcallpointmass
      real, dimension(3) :: positions,acceleration
!
      real, dimension (mparray), optional :: fp_pt
      real, dimension (mpvar), optional :: dfp_pt
      logical, optional :: flag_pt
!
      if (lcallpointmass) then
        positions    =  fq(k,ixq:izq)
        acceleration = dfq(k,ivxq:ivzq)  
      else
        !if (lparticles) then
        !  ip1=1;ip3=3
        !  iv1=4;iv3=6
        !else
        !  ip1=0;ip3=0
        !  iv1=0;iv3=0
        !endif
        positions    =  fp_pt(index%ixw:index%izw)
        acceleration = dfp_pt(index%ivxw:index%ivzw)
      endif
!
      do ks=1,nqpar
        if (.not.(lcallpointmass.and.k==ks)) then !prevent self-acceleration
!
          call get_evr(positions,fq(ks,ixq:izq),evr_cart)
!
!  Particles relative distance from each other:
!
!  r_ij = sqrt(ev1**2 + ev2**2 + ev3**2)
!
          tmp1=sum(evr_cart**2)
          tmp2=r_smooth(ks)**2
          r2_ij=max(tmp1,tmp2)
!
!  If there is accretion, remove the accreted particles from the simulation, if any.
!
          if (.not.(lcallpointmass).and.laccretion(ks)) then
            rs2=0.01 !(accrete_hills_frac(ks)**2)*hill_radius_square(ks)
            if (r2_ij<=rs2) then
              !flag particle for removal 
               flag_pt=.true.
              !add mass of the removed particle to the accreting particle
              !if (ladd_mass(ks)) pmass(ks)=pmass(ks)+mp_swarm
              return
            endif
          endif
!
!  Shortcut: invr3_ij = r_ij**(-3)
!
          if (r2_ij > 0) then
            invr3_ij = r2_ij**(-1.5)
          else                ! can happen during pencil_check
            invr3_ij = 0.0
          endif
!
!  Gravitational acceleration: g=-g0/|r-r0|^3 (r-r0)
!
!  The acceleration is in non-coordinate basis (all have dimension of length).
!  The main dxx_dt of particle_dust takes care of transforming the linear
!  velocities to angular changes in position.
!
          acc_cart = - GNewton*pmass(ks)*invr3_ij*evr_cart(1:3)
          !if (lcartesian_coords) then
            acceleration =  acceleration + acc_cart
          !else
            ! separate this N-body acceleration from other, added elsewhere in the code
          !  dfq(k,ivxq_cart:ivzq_cart) =  dfq(k,ivxq_cart:ivzq_cart) + acc_cart
          !endif
!
!  Time-step constraint from N-body particles. We use both the criterion
!  that the distance to the N-body particle must not change too much in
!  one time-step and additionally we use the free-fall time-scale.
!
          if (lfirst.and.ldt.and.ldt_pointmasses) then
            if (.not.lcallpointmass) then     
              v_ij=sqrt(sum((fp_pt(ivxq:ivzq)-fq(ks,ivxq:ivzq))**2))
              dt1_max(l1-nghost)= &
                   max(dt1_max(l1-nghost),v_ij/sqrt(r2_ij)/cdtq)
            endif
            dt1_max(l1-nghost)= &
                 max(dt1_max(l1-nghost), &
                 sqrt(GNewton*pmass(ks)*invr3_ij)/cdtq)
          endif
!
        endif !if (ipar(k)/=ks)
!
      enddo !nbody loop
!
      if (lcallpointmass) then
        dfq(k,ivxq:ivzq) = acceleration
      else
        dfp_pt(index%ivxw:index%ivzw) = acceleration
      endif
!
    endsubroutine loop_through_pointmasses_inertial
!**********************************************************
    subroutine get_evr(xxp,xxq,evr_cart)
!
!  Point-to-point vector distance, in different coordinate systems.
!  Return always in Cartesian.
!
!  14-feb-14/wlad: coded
!
      real, dimension(3), intent(in) :: xxp,xxq
      real, dimension(3), intent(out) :: evr_cart
      real :: x1,y1,x2,y2,z1,z2
!
      if (lcartesian_coords) then
        x1=xxp(1) ; x2=xxq(1)
        y1=xxp(2) ; y2=xxq(2)
        z1=xxp(3) ; z2=xxq(3)
        evr_cart(1)=x1-x2
        evr_cart(2)=y1-y2
        evr_cart(3)=z1-z2
      elseif (lcylindrical_coords) then
        !x1=xxp(1)*cos(xxp(2))
        !y1=xxp(1)*sin(xxp(2))
        z1=xxp(3)
!
        !x2=xxq(1)*cos(xxq(2))
        !y2=xxq(1)*sin(xxq(2))
        z2=xxq(3)

         evr_cart(1)=xxp(1) - xxq(1)*cos(xxp(2)-xxq(2))
         evr_cart(2)=xxq(1)*sin(xxp(2)-xxq(2))
         evr_cart(3)=z1-z2

      elseif (lspherical_coords) then
        x1=xxp(1)*sin(xxp(2))*cos(xxp(3))
        y1=xxp(1)*sin(xxp(2))*sin(xxp(3))
        z1=xxp(1)*cos(xxp(2))
!
        x2=xxq(1)*sin(xxq(2))*cos(xxq(3))
        y2=xxq(1)*sin(xxq(2))*sin(xxq(3))
        z2=xxq(1)*cos(xxq(2))
      endif
!
      !evr_cart(1)=x1-x2
      !evr_cart(2)=y1-y2
      !evr_cart(3)=z1-z2
!
    endsubroutine get_evr
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
    subroutine read_pointmasses_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=pointmasses_init_pars, IOSTAT=iostat)
!
    endsubroutine read_pointmasses_init_pars
!***********************************************************************
    subroutine write_pointmasses_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=pointmasses_init_pars)
!
    endsubroutine write_pointmasses_init_pars
!***********************************************************************
    subroutine read_pointmasses_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=pointmasses_run_pars, IOSTAT=iostat)
!
    endsubroutine read_pointmasses_run_pars
!***********************************************************************
    subroutine write_pointmasses_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=pointmasses_run_pars)
!
    endsubroutine write_pointmasses_run_pars
!***********************************************************************
    subroutine reset_center_of_mass()
!
!  If the center of mass was accelerated, reset its position
!  to the center of the grid.
!
!  Assumes that the total mass of the particles is one.
!
!  27-aug-06/wlad: coded
!  18-mar-08/wlad: cylindrical and spherical corrections
!
      !real, dimension(nqpar,mqvar),intent(inout) :: dfq
      real, dimension(nqpar,mqarray) :: ftmp
      real, dimension(3) :: vcm
      real :: xcm,ycm,zcm,thtcm,phicm,vxcm,vycm,vzcm
      integer :: ks
!
      ftmp=fq(1:nqpar,:)
!
      if (lcartesian_coords) then
        vcm(1) = sum(ftmp(:,imass)*ftmp(:,ivxq))
        vcm(2) = sum(ftmp(:,imass)*ftmp(:,ivyq))
        vcm(3) = sum(ftmp(:,imass)*ftmp(:,ivzq))
      else if (lcylindrical_coords) then
        xcm=sum(ftmp(:,imass)*(ftmp(:,ixq)*cos(ftmp(:,iyq))))
        ycm=sum(ftmp(:,imass)*(ftmp(:,ixq)*sin(ftmp(:,iyq))))
        phicm=atan2(ycm,xcm)
!
        vxcm=sum(ftmp(:,imass)*(&
            ftmp(:,ivxq)*cos(ftmp(:,iyq))-ftmp(:,ivyq)*sin(ftmp(:,iyq))))
        vycm=sum(ftmp(:,imass)*(&
            ftmp(:,ivxq)*sin(ftmp(:,iyq))+ftmp(:,ivyq)*cos(ftmp(:,iyq))))
!
        vcm(1)= vxcm*cos(phicm) + vycm*sin(phicm)
        vcm(2)=-vxcm*sin(phicm) + vycm*cos(phicm)
        vcm(3) = sum(ftmp(:,imass)*ftmp(:,ivzq))
!
      else if (lspherical_coords) then
        vxcm=sum(ftmp(:,imass)*( &
              ftmp(:,ivxq)*sin(ftmp(:,iyq))*cos(ftmp(:,izq))&
             +ftmp(:,ivyq)*cos(ftmp(:,iyq))*cos(ftmp(:,izq))&
             -ftmp(:,ivzq)*sin(ftmp(:,izq))                ))
        vycm=sum(ftmp(:,imass)*( &
              ftmp(:,ivxq)*sin(ftmp(:,iyq))*sin(ftmp(:,izq))&
             +ftmp(:,ivyq)*cos(ftmp(:,iyq))*sin(ftmp(:,izq))&
             +ftmp(:,ivzq)*cos(ftmp(:,izq))                ))
        vzcm=sum(ftmp(:,imass)*(&
             ftmp(:,ivxq)*cos(ftmp(:,iyq))-ftmp(:,ivyq)*sin(ftmp(:,iyq))))
!
        xcm=sum(ftmp(:,imass)*(ftmp(:,ixq)*sin(ftmp(:,iyq))*cos(ftmp(:,izq))))
        ycm=sum(ftmp(:,imass)*(ftmp(:,ixq)*sin(ftmp(:,iyq))*sin(ftmp(:,izq))))
        zcm=sum(ftmp(:,imass)*(ftmp(:,ixq)*cos(ftmp(:,iyq))))
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
      do ks=1,nqpar
        dfq(ks,ixq:izq) = dfq(ks,ixq:izq) - vcm*totmass1
      enddo
!
    endsubroutine reset_center_of_mass
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
    subroutine get_gravity_field_pointmasses(grr,gg,ks)
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
      x0=fq(ks,ixq);y0=fq(ks,iyq);z0=fq(ks,izq)
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
    endsubroutine get_gravity_field_pointmasses
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
      real, dimension(nx) :: torque_gas,torqint_gas,torqext_gas
      real, dimension(nx) :: torque_par,torqint_par,torqext_par
      real, dimension(nx) :: dist,rpre,rr_mn,tempering
      real :: rr,w2,smap,hills,phip,phi,pcut
      integer :: ks,i
!
      if (ks==istar) call fatal_error('calc_torque', &
          'Nonsense to calculate torques for the star')
!
      if (lcartesian_coords) then
        rr    = sqrt(fq(ks,ixq)**2 + fq(ks,iyq)**2 + fq(ks,izq)**2)
        rpre  = fq(ks,ixq)*y(m) - fq(ks,iyq)*x(l1:l2)
      elseif (lcylindrical_coords) then
        rr_mn = x(l1:l2)     ; phi  = y(m)
        rr    = fq(ks,ixq)  ; phip = fq(ks,iyq)
        rpre  = rr_mn*rr*sin(phi-phip)
      elseif (lspherical_coords) then
        call fatal_error("calc_torque",&
             "not yet implemented for spherical coordinates")
      else
        call fatal_error("calc_torque",&
             "the world is flat and we should never gotten here")
      endif
!
      w2    = fq(ks,ivxq)**2 + fq(ks,ivyq)**2 + fq(ks,ivzq)**2
      smap  = 1./(2./rr - w2)
      hills = smap*(pmass(ks)*pmass1(istar)/3.)**(1./3.)
!
!  Define separate torques for gas and dust/particles
!
      torque_gas = GNewton*pmass(ks)*p%rho*rpre*&
           (dist**2 + r_smooth(ks)**2)**(-1.5)
!
      if (ldust) then
        torque_par = GNewton*pmass(ks)*p%rhop*rpre*&
             (dist**2 + r_smooth(ks)**2)**(-1.5)
      else
        torque_par=0.
      endif
!
      if (ldustdensity) &
           call fatal_error("calc_torque",&
           "not implemented for the dust fluid approximation")
!
!  Zero torque outside r_int and r_ext in Cartesian coordinates
!
      if (lcartesian_coords) then
        do i=1,nx
          if (p%rcyl_mn(i) > r_ext .or. p%rcyl_mn(i) < r_int) then
            torque_gas(i)=0.
            torque_par(i)=0.
          endif
        enddo
      endif
!
!  Exclude region inside a fraction (hills_tempering_fraction) of the Hill sphere.
!
      if (ltempering) then
        pcut=hills_tempering_fraction*hills
        tempering = 1./(exp(-(sqrt(dist**2)/hills - pcut)/(.1*pcut))+1.)
        torque_gas = torque_gas * tempering
        torque_par = torque_par * tempering
      else
        do i=1,nx
          if (dist(i)<hills) then
            torque_gas(i)=0.
            torque_par(i)=0.
          endif
        enddo
      endif
!
!  Separate internal and external torques
!
      do i=1,nx
        if (p%rcyl_mn(i)>=rr) then
          torqext_gas(i) = torque_gas(i)
          torqext_par(i) = torque_par(i)
        else
          torqext_gas(i)=0.;torqext_par(i)=0.
        endif
        if (p%rcyl_mn(i)<=rr) then
          torqint_gas(i) = torque_gas(i)
          torqint_par(i) = torque_par(i)
        else
          torqint_gas(i)=0.;torqint_par(i)=0.
        endif
      enddo
!
!  Sum the different torque contributions.
!
      if (lcartesian_coords) then
        call sum_lim_mn_name(torqext_gas,idiag_torqext_gas(ks),p)
        call sum_lim_mn_name(torqext_par,idiag_torqext_par(ks),p)
        call sum_lim_mn_name(torqint_gas,idiag_torqint_gas(ks),p)
        call sum_lim_mn_name(torqint_par,idiag_torqint_par(ks),p)
!
!  Backward compatibility
!
        call sum_lim_mn_name(torqext_gas+torqext_par,idiag_torqext(ks),p)
        call sum_lim_mn_name(torqint_gas+torqint_par,idiag_torqint(ks),p)
      else
        !
        ! Hack for non-cartesian coordinates. sum_lim_mn_name is lagging
        ! behind sum_mn_name, and whould be brought up to date.
        !
        call integrate_mn_name(torqext_gas,idiag_torqext_gas(ks))
        call integrate_mn_name(torqext_par,idiag_torqext_par(ks))
        call integrate_mn_name(torqint_gas,idiag_torqint_gas(ks))
        call integrate_mn_name(torqint_par,idiag_torqint_par(ks))
        call integrate_mn_name(torqext_gas+torqext_par,idiag_torqext(ks))
        call integrate_mn_name(torqint_gas+torqint_par,idiag_torqint(ks))
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
        do ks=1,nqpar !_orig
          if (ks/=istar) then
            pmass(ks)= max(dble(tini),&
                 final_ramped_mass(ks)*(sin((.5*pi)*(t/ramping_period))**2))
            tmp=tmp+pmass(ks)
          endif
        enddo
        pmass(istar)= 1-tmp
      else
         !pmass(1:nqpar_orig)=final_ramped_mass(1:nqpar_orig)
         pmass(1:nqpar)=final_ramped_mass(1:nqpar)
      endif
!
    endsubroutine get_ramped_mass
!***********************************************************************
    subroutine get_total_gravity(ggt)
!
!  Sum the gravities of all massive particles
!
!  08-mar-08/wlad: coded
!
      use Sub
!
      real, dimension (mx,nqpar) :: rp_mn,rpcyl_mn
      real, dimension (mx,3)     :: ggp,ggt
      real, dimension (mx)       :: grav_particle,rrp
      integer                    :: ks
!
      intent(out) :: ggt
!
      ggt=0.
      do ks=1,nqpar
!
!  Spherical and cylindrical distances
!
        call get_radial_distance(rp_mn(:,ks),rpcyl_mn(:,ks),&
             e1_=fq(ks,ixq),e2_=fq(ks,iyq),e3_=fq(ks,izq))
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
        call get_gravity_field_pointmasses(grav_particle,ggp,ks)
!
        if ((ks==istar).and.lnogravz_star) &
            ggp(:,3) = 0.
!
!  Sum up the accelerations of the massive particles
!
        ggt=ggt+ggp
!
!  Add indirect term if the star is fixed at the center
!
        if (.not.linertial_frame) call add_indirect_term(ks,ggt)
!
      enddo
!
    endsubroutine get_total_gravity
!***********************************************************************
    subroutine add_indirect_term(ks,ggt)
!
!  Add the terms due to the reference frame being away from the baricenter.
!  So far, only for one perturber (two massive bodies), and in Cartesian coordinates.
!
!  23-jun-14/wlad: coded
!
      real, dimension(mx,3) :: ggt
      real :: rr1_3
      integer, intent(in) :: ks
!
      if (ks/=istar) then
        rr1_3=(fq(ks,ixq)**2 + fq(ks,iyq)**2 + fq(ks,izq)**2)**(-1.5)
        ggt(:,1) = ggt(:,1) - GNewton*pmass(ks)*fq(ks,ixq)*rr1_3
        ggt(:,2) = ggt(:,2) - GNewton*pmass(ks)*fq(ks,iyq)*rr1_3
        ggt(:,3) = ggt(:,3) - GNewton*pmass(ks)*fq(ks,izq)*rr1_3
      endif
!
    endsubroutine add_indirect_term
!***********************************************************************
    subroutine pointmasses_read_snapshot(filename)
!
!  Read nbody particle info
!
!  01-apr-08/wlad: dummy
!
      use Mpicomm, only: mpibcast_real
!
      character (len=*) :: filename
      integer :: nqpar_read
!
      if (lroot) then
        open(1,FILE=filename,FORM='unformatted')
        print*,'opened file'
        read(1) nqpar_read
        if (nqpar_read /= nqpar) call fatal_error("","")
        if (nqpar_read/=0) read(1) fq
        if (ip<=8) print*, 'read snapshot', filename
        close(1)
      endif
!
      call mpibcast_real(fq,(/nqpar,mqarray/))
!
      if (ldebug) then
        print*,'pointmasses_read_snapshot'
        print*,'nqpar=',nqpar
        print*,'fq =',fq
        print*,''
      endif
!
    endsubroutine pointmasses_read_snapshot
!***********************************************************************
    subroutine pointmasses_write_snapshot(snapbase,enum,flist)
!
      use General, only:safe_character_assign
      use Sub, only: update_snaptime, read_snaptime
      use IO, only: lun_output
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
          write(lun_output) nqpar
          if (nqpar/=0) write(lun_output) fq
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
        write(lun_output) nqpar
        if (nqpar/=0) write(lun_output) fq 
        close(lun_output)
        if (ip<=10 .and. lroot) &
             print*,'written snapshot ', snapname
      endif
!
      call keep_compiler_quiet(flist)
!
    endsubroutine pointmasses_write_snapshot
!***********************************************************************
    subroutine pointmasses_write_qdim(filename)
!
!  Write nqpar and mqvar to file.
!
!  01-apr-08/wlad: coded
!
      character (len=*) :: filename
!
      open(1,file=filename)
      write(1,'(3i9)') nqpar, mqvar
      close(1)
!
    endsubroutine pointmasses_write_qdim
!***********************************************************************
    subroutine rprint_pointmasses(lreset,lwrite)
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
        write(3,*) 'ixq=',ixq 
        write(3,*) 'iyq=',iyq 
        write(3,*) 'izq=',izq
        write(3,*) 'ivxq=',ivxq 
        write(3,*) 'ivyq=',ivyq 
        write(3,*) 'ivzq=',ivzq
        write(3,*) 'imass=', imass
      endif
!
!  Reset everything in case of reset
!
      if (lreset) then
        idiag_xxq=0;idiag_vvq=0
        idiag_torqint=0;idiag_torqext=0
        idiag_totenergy=0
        idiag_torqext_gas=0;idiag_torqext_par=0
        idiag_torqint_gas=0;idiag_torqint_par=0
      endif
!
!  Run through all possible names that may be listed in print.in
!
      if (lroot .and. ip<14) print*,'rprint_pointmasses: run through parse list'
!
!  Now check diagnostics for specific particles
!
      do ks=1,nqpar
        sks=itoa(ks)
        do j=1,3
          if (j==1) str='x';if (j==2) str='y';if (j==3)  str='z'
          do iname=1,nname
            call parse_name(iname,cname(iname),cform(iname),&
                 trim(str)//'q'//trim(sks),idiag_xxq(ks,j))
            call parse_name(iname,cname(iname),cform(iname),&
                 'v'//trim(str)//'q'//trim(sks),idiag_vvq(ks,j))
          enddo
!
!  Run through parse list again
!
          if (lwr) then
            write(3,*) ' i_'//trim(str)//'q'//trim(sks)//'=',&
                 idiag_xxq(ks,j)
            write(3,*) 'i_v'//trim(str)//'q'//trim(sks)//'=',&
                 idiag_vvq(ks,j)
          endif
!
        enddo
!
        do iname=1,nname
          call parse_name(iname,cname(iname),cform(iname),&
               'torqint_'//trim(sks),idiag_torqint(ks))
          call parse_name(iname,cname(iname),cform(iname),&
               'torqext_'//trim(sks),idiag_torqext(ks))
          call parse_name(iname,cname(iname),cform(iname),&
               'torqext_gas_'//trim(sks),idiag_torqext_gas(ks))
          call parse_name(iname,cname(iname),cform(iname),&
               'torqext_par_'//trim(sks),idiag_torqext_par(ks))
          call parse_name(iname,cname(iname),cform(iname),&
               'torqint_gas_'//trim(sks),idiag_torqint_gas(ks))
          call parse_name(iname,cname(iname),cform(iname),&
               'torqint_par_'//trim(sks),idiag_torqint_par(ks))
        enddo
!
        if (lwr) then
          write(3,*) 'i_torqint_'//trim(sks)//'=',idiag_torqint(ks)
          write(3,*) 'i_torqext_'//trim(sks)//'=',idiag_torqext(ks)
          write(3,*) 'i_torqint_gas'//trim(sks)//'=',idiag_torqint(ks)
          write(3,*) 'i_torqext_gas'//trim(sks)//'=',idiag_torqext(ks)
          write(3,*) 'i_torqint_par'//trim(sks)//'=',idiag_torqint(ks)
          write(3,*) 'i_torqext_par'//trim(sks)//'=',idiag_torqext(ks)
        endif
      enddo
!
!  Diagnostic related to quantities summed over all point masses
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'totenergy',idiag_totenergy)
      enddo
!
       if (lwr) then
         write(3,*) 'i_totenergy=',idiag_totenergy
       endif
!
    endsubroutine rprint_pointmasses
!***********************************************************************
    subroutine boundconds_pointmasses!(f)
!
!  Global boundary conditions for particles.
!
!  30-dec-04/anders: coded
!
      use Mpicomm
      use General, only: random_number_wrapper
      !use Particles_mpicomm
      use SharedVariables, only: get_shared_variable
!
      !real, dimension (mpar_loc,mparray) :: fp
      !integer, dimension (mpar_loc) :: ipar
      !real, dimension (mpar_loc,mpvar), optional :: dfp
      !logical, optional :: linsert
!
      real :: xold, yold, rad, r1old, OO, tmp
      integer :: k, ik, k1, k2
      character (len=2*bclen+1) :: boundx, boundy, boundz
!
      !intent (inout) :: fp, ipar, dfp
!
      do k=1,nqpar
!
         boundx=bcqx; boundy=bcqy; boundz=bcqz
!
!  Cartesian boundaries: Boundary condition in the x-direction. The physical
!  domain is in the interval
!
!    x \in [x0,x1[
!    y \in [y0,y1[
!    z \in [z0,z1[
!
        if (nxgrid/=1) then
          if (boundx=='p') then
!  xp < x0
            if (fq(k,ixq)< xyz0(1)) then
              fq(k,ixq)=fq(k,ixq)+Lxyz(1)
!
!  Particle position must never need more than one addition of Lx to get back
!  in the box. Often a NaN or Inf in the particle position will show up as a
!  problem here.
              if (fq(k,ixq)< xyz0(1)) then
                print*, 'boundconds_pointmasses: ERROR - particle ', k, &
                     ' was further than Lx outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'xxq=', fq(k,ixq:izq)
                call fatal_error_local('boundconds_pointmasses','')
              endif
            endif
!  xp > x1
            if (fq(k,ixq)>=xyz1(1)) then
              fq(k,ixq)=fq(k,ixq)-Lxyz(1)
!
              if (fq(k,ixq)>=xyz1(1)) then
                print*, 'boundconds_pointmasses: ERROR - particle ', k, &
                     ' was further than Lx outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'xxq=', fq(k,ixq:izq)
                call fatal_error_local('boundconds_pointmasses','')
              endif
            endif
          elseif (boundx=='out') then
!
!  Do nothing. A massive particle, can be out of the box. A star, for example,
!  in a cylindrical simulation
!
          else
            print*, 'boundconds_pointmasses: No such boundary condition =', boundx
            call stop_it('boundconds_pointmasses')
          endif
        endif
!
!  Boundary condition in the y-direction.
!
        if (nygrid/=1) then
          if (boundy=='p') then
!  yp < y0
            if (fq(k,iyq)< xyz0(2)) then
              fq(k,iyq)=fq(k,iyq)+Lxyz(2)
              if (fq(k,iyq)< xyz0(2)) then
                print*, 'boundconds_pointmasses: ERROR - particle ', k, &
                     ' was further than Ly outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'xxq=', fq(k,ixq:izq)
                call fatal_error_local('boundconds_pointmasses','')
              endif
            endif
!  yp > y1
            if (fq(k,iyq)>=xyz1(2)) then
              fq(k,iyq)=fq(k,iyq)-Lxyz(2)
              if (fq(k,iyq)>=xyz1(2)) then
                print*, 'boundconds_pointmasses: ERROR - particle ', k, &
                     ' was further than Ly outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'xxq=', fq(k,ixq:izq)
                call fatal_error_local('boundconds_pointmasses','')
              endif
            endif
          elseif (boundy=='out') then
            ! massive particles can be out of the box
            ! the star, for example, in a cylindrical simulation
          else
            print*, 'boundconds_pointmasses: No such boundary condition =', boundy
            call stop_it('boundconds_pointmasses')
          endif
        endif
!
!  Boundary condition in the z-direction.
!
        if (nzgrid/=1) then
          if (boundz=='p') then
!  zp < z0
            if (fq(k,izq)< xyz0(3)) then
              fq(k,izq)=fq(k,izq)+Lxyz(3)
              if (fq(k,izq)< xyz0(3)) then
                print*, 'boundconds_pointmasses: ERROR - particle ', k, &
                     ' was further than Lz outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'xxq=', fq(k,ixq:izq)
                call fatal_error_local('boundconds_pointmasses','')
              endif
            endif
!  zp > z1
            if (fq(k,izq)>=xyz1(3)) then
              fq(k,izq)=fq(k,izq)-Lxyz(3)
              if (fq(k,izq)>=xyz1(3)) then
                print*, 'boundconds_pointmasses: ERROR - particle ', k, &
                     ' was further than Lz outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'xxq=', fq(k,ixq:izq)
                call fatal_error_local('boundconds_pointmasses','')
              endif
            endif
          elseif (boundz=='out') then
            ! massive particles can be out of the box
            ! the star, for example, in a cylindrical simulation
          else
            print*, 'boundconds_pointmasses: No such boundary condition=', boundz
            call stop_it('boundconds_pointmasses')
          endif
        endif
      enddo
!
    endsubroutine boundconds_pointmasses
!***********************************************************************
    subroutine pointmasses_timestep_first(f)
!    
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (lroot) then 
        if (lfirst) then
          dfq(1:nqpar,:)=0.0
        else
          dfq(1:nqpar,:)=alpha_ts(itsub)*dfq(1:nqpar,:)            
        endif
      endif
!
   endsubroutine pointmasses_timestep_first
!***********************************************************************
    subroutine pointmasses_timestep_second(f)
!
!  Time evolution of particle variables.
!
!  07-jan-05/anders: coded
!
      use Mpicomm, only: mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (lroot) & 
           fq(1:nqpar,1:mqvar) = fq(1:nqpar,1:mqvar) + dt_beta_ts(itsub)*dfq(1:nqpar,1:mqvar)
      call mpibcast_real(fq,(/nqpar,mqarray/))
      !else
      !  call advance_particles_in_cartesian(fp,dfp)
      !endif
!
    endsubroutine pointmasses_timestep_second
!***********************************************************************
    subroutine correct_curvilinear
!
!  Curvilinear corrections to acceleration only.
!  Corrections to velocity were already taken into account
!  in the dxx_dt of particles_dust.f90
!
!  In the case that the N-body code is used, the update in polar grids
!  in done by transforming the variables first to Cartesian, to achieve a
!  better conservation of the Jacobi constant, and this code is not called.
!
!  15-sep-07/wlad: coded
!
      real :: rad,raddot,phidot,thtdot,sintht,costht
      integer :: k
!
      do k=1,nqpar
!
!  Correct acceleration for curvilinear coordinates.
!
        if (lcylindrical_coords) then
          rad=fq(k,ixq);raddot=fq(k,ivxq);phidot=fq(k,ivyq)/max(rad,tini)
          dfq(k,ivxq) = dfq(k,ivxq) + rad*phidot**2
          dfq(k,ivyq) = dfq(k,ivyq) - 2*raddot*phidot
        elseif (lspherical_coords) then
          rad=fq(k,ixq)
          sintht=sin(fq(k,iyq));costht=cos(fq(k,iyq))
          raddot=fq(k,ivxq);thtdot=fq(k,ivyq)/max(rad,tini)
          phidot=fq(k,ivzq)/(max(rad,tini)*sintht)
!
          dfq(k,ivxq) = dfq(k,ivxq) &
               + rad*(thtdot**2 + (sintht*phidot)**2)
          dfq(k,ivyq) = dfq(k,ivyq) &
               - 2*raddot*thtdot + rad*sintht*costht*phidot**2
          dfq(k,ivzq) = dfq(k,ivzq) &
               - 2*phidot*(sintht*raddot + rad*costht*thtdot)
        endif
      enddo
!
    endsubroutine correct_curvilinear
!***********************************************************************
  endmodule PointMasses
