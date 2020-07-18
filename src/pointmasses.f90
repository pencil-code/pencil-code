! $Id: pointmasses.f90,v 1.1 2019/02/02 03:54:41 wlyra Exp $
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
! MPVAR CONTRIBUTION 3
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
  character (len=labellen), dimension(mqarray) :: qvarname
  real, dimension(nqpar,mqarray) :: fq=0.0
  real, dimension(nqpar,mqvar) :: dfq=0.0
  real, dimension(nqpar,3) :: dfq_cart=0.0
  real, dimension(nqpar) :: xq0=0.0, yq0=0.0, zq0=0.0
  real, dimension(nqpar) :: vxq0=0.0, vyq0=0.0, vzq0=0.0
  real, dimension(nqpar) :: pmass=0.0, r_smooth=impossible, pmass1
  real, dimension(nqpar) :: frac_smooth=0.4
  real, dimension(nqpar) :: accrete_hills_frac=0.2, final_ramped_mass=0.0
  real :: eccentricity=0.0, semimajor_axis=1.0
  real :: totmass, totmass1
  real :: GNewton1, GNewton=impossible, density_scale=0.001
  real :: cdtq=0.1
  real :: hills_tempering_fraction=0.8
  real, pointer :: rhs_poisson_const, tstart_selfgrav
  integer :: ramp_orbits=5
  integer :: iprimary=1, isecondary=2
  integer :: ivpx_cart=0,ivpy_cart=0,ivpz_cart=0
!
  integer :: imass=0, ixq=0, iyq=0, izq=0, ivxq=0, ivyq=0, ivzq=0
  integer :: nqvar=0
!
  logical, dimension(nqpar) :: lcylindrical_gravity_nbody=.false.
  logical, dimension(nqpar) :: lfollow_particle=.false., laccretion=.false.
  logical :: llive_secondary=.false., lnorm=.true.
  logical :: lreset_cm=.true., lnogravz_star=.false., lexclude_frozen=.false.
  logical :: lramp=.false.
  logical :: ldt_pointmasses=.true.
  logical :: ldust=.false.
  logical :: ltempering=.false.
  logical :: lretrograde=.false.
  logical :: lnoselfgrav_primary=.true.
  logical :: lgas_gravity=.true.,ldust_gravity=.false.
  logical :: lcorrect_gasgravity_lstart=.false.
  logical :: lexclude_hills=.false.
!
  character (len=labellen) :: initxxq='random', initvvq='nothing'
  character (len=labellen), dimension (nqpar) :: ipotential_pointmass='newton'
  character (len=2*bclen+1) :: bcqx='p', bcqy='p', bcqz='p'
!
  logical :: ladd_dragforce=.false.,lquadratic_drag=.false.,llinear_drag=.true.
  logical :: lcoriolis_force=.false.
  logical :: l2D,l3D
  real :: ugas=0.0,Omega_coriolis=0.0
  real, dimension(nqpar) :: StokesNumber=1.
!
  type IndexDustParticles
    integer :: ixw=0,iyw=0,izw=0
    integer :: ivxw=0,ivyw=0,ivzw=0
  endtype IndexDustParticles
  type (IndexDustParticles), save :: index
!
  namelist /pointmasses_init_pars/ &
      initxxq, initvvq, xq0, yq0, zq0, vxq0, vyq0, vzq0,  &
      pmass, r_smooth, lcylindrical_gravity_nbody, lexclude_frozen, GNewton, &
      bcqx, bcqy, bcqz, ramp_orbits, lramp, final_ramped_mass, &
      laccretion, &
      accrete_hills_frac, iprimary, &
      ldt_pointmasses, cdtq, lretrograde, &
      eccentricity, semimajor_axis, & 
      ipotential_pointmass, density_scale,&
      lgas_gravity,ldust_gravity,lcorrect_gasgravity_lstart,&
      frac_smooth
!
  namelist /pointmasses_run_pars/ &
      lreset_cm, &
      lnogravz_star, lfollow_particle, llive_secondary, lexclude_frozen, &
      GNewton, bcqx, bcqy, bcqz, &
      laccretion, accrete_hills_frac, iprimary, &
      ldt_pointmasses, cdtq, hills_tempering_fraction, &
      ltempering, & 
      ipotential_pointmass, density_scale,&
      lgas_gravity,ldust_gravity,&
      ladd_dragforce,ugas,StokesNumber,&
      lquadratic_drag,llinear_drag,lcoriolis_force,Omega_coriolis,&
      frac_smooth,lexclude_hills
!
  integer, dimension(nqpar,3) :: idiag_xxq=0,idiag_vvq=0
  integer, dimension(nqpar)   :: idiag_torqint=0,idiag_torqext=0
  integer, dimension(nqpar)   :: idiag_torqext_gas=0,idiag_torqext_par=0
  integer, dimension(nqpar)   :: idiag_torqint_gas=0,idiag_torqint_par=0
  integer, dimension(nqpar)   :: idiag_period=0,idiag_torque=0
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
      use Particles_main, only: append_particle_index
!      
      integer :: iqvar
!
      if (lroot) call svn_id( &
          "$Id: pointmasses.f90,v 1.1 2019/02/02 03:54:41 wlyra Exp $")
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
!  Auxiliary variables for polar coordinates
!
      ixq = nqvar+1
      qvarname(ixq)='ixq'
      iyq = nqvar+2
      qvarname(iyq)='iyq'
      izq = nqvar+3
      qvarname(izq)='izq'
      nqvar=nqvar+3
!
      ivxq = nqvar+1
      qvarname(ivxq)='ivxq'
      ivyq = nqvar+2
      qvarname(ivyq)='ivyq'
      ivzq = nqvar+3
      qvarname(ivzq)='ivzq'
      nqvar=nqvar+3
!
!  Set up mass as particle index. Plus seven, since the other 6 are
!  used by positions and velocities.
!
      imass=nqvar+1
      qvarname(imass)='imass'
      nqvar=nqvar+1
!
!  Check that the fq and dfq arrays are big enough.
!
      if (nqvar > mqvar) then
        if (lroot) write(0,*) 'nqvar = ', nqvar, ', mqvar = ', mqvar
        call fatal_error('register_pointmasses','nqvar > mqvar')
      endif
!
      if (lparticles) then
        call append_particle_index('ivpx_cart',ivpx_cart)
        call append_particle_index('ivpy_cart',ivpy_cart)
        call append_particle_index('ivpz_cart',ivpz_cart)
      endif
!
      if (lroot) then
        open(3,file=trim(datadir)//'/qvarname.dat',status='replace')
        do iqvar=1,mqarray
          write(3,"(i4,2x,a)") iqvar, qvarname(iqvar)
        enddo
        close(3)
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
      logical :: lxpresent,lypresent,lzpresent,l2Dcyl,l2Dsph
!
      integer :: ks
!
!  Look for initialized masses.
!
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
      call rprint_pointmasses(.false.,LWRITE=lroot)
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
            if (ks/=iprimary) pmass(ks) = epsi
          enddo
          pmass(iprimary)=1-epsi*(nqpar-1)
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
      if (totmass/=1) &
           call warning("initialize_pointmasses","The masses do not sum up to one!")
      totmass1=1./max(totmass, tini)
!
!  Check for consistency between the cylindrical gravities switches
!  from cdata and the one from point mass.
!
      if (((lcylindrical_gravity).and.&
        (.not.lcylindrical_gravity_nbody(iprimary))).or.&
             (.not.lcylindrical_gravity).and.&
             (lcylindrical_gravity_nbody(iprimary))) then
        call fatal_error('initialize_pointmasses','inconsistency '//&
            'between lcylindrical_gravity from cdata and the '//&
            'one from point mass')
      endif
!
!  Smoothing radius 
!
      if (any(r_smooth == impossible)) then 
        do ks=1,nqpar
          if (ks/=iprimary) then
            r_smooth(ks) = frac_smooth(ks) * xq0(ks) * (pmass(ks)/3.)**(1./3)
          else
            r_smooth(ks) = rsmooth
          endif
        enddo
      endif
!
      if (rsmooth/=r_smooth(iprimary)) then
        print*,'rsmooth from cdata=',rsmooth
        print*,'r_smooth(iprimary)=',r_smooth(iprimary)
        call fatal_error('initialize_pointmasses','inconsistency '//&
            'between rsmooth from cdata and the '//&
            'one from point mass')
      endif
!
!  The presence of dust particles needs to be known.
!
      if (npar > nqpar) ldust=.true.
!
!  Check if the run is 2D or 3D to use cylindrical gravity or not
!
      lxpresent=(nxgrid/=1)
      lypresent=(nygrid/=1)
      lzpresent=(nzgrid/=1)
!
      l3D  =                              lxpresent.and.      lypresent .and.      lzpresent
      l2Dcyl=    lcylindrical_coords.and.(lxpresent.and.      lypresent .and.(.not.lzpresent))
      l2Dsph=    lspherical_coords  .and.(lxpresent.and.(.not.lypresent).and.      lzpresent )
      l2D=l2Dcyl.or.l2Dsph
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
      if (ldust.and.llive_secondary) &
          lpenc_requested(i_rhop)=.true.
!
      if (llive_secondary.and.(.not.lselfgravity)) then
        if (l2D.or.(l3D.and.lcylindrical_gravity)) then
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
      if (l2D.or.(l3D.and.lcylindrical_gravity)) then
        lpenc_diagnos(i_rcyl_mn)=.true.
      else
        lpenc_diagnos(i_r_mn)=.true.
      endif
      if (ldust) lpenc_diagnos(i_rhop)=.true.
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
      real :: absolute_offset_star,baricenter_secondaries
      real :: velocity_baricenter_secondaries,mass_secondaries
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
        mass_secondaries = 0.
        baricenter_secondaries=0.
        do ks=1,nqpar
          if (ks/=iprimary) then
            sma(ks) = abs(positions(ks,1))
            mass_secondaries = mass_secondaries + pmass(ks)
            baricenter_secondaries = baricenter_secondaries + positions(ks,1)*pmass(ks)
          endif
        enddo
        absolute_offset_star = abs(baricenter_secondaries)
!
!  Fixed-cm assumes that the total mass is always one. The mass of the
!  star is adjusted to ensure this.
!
        pmass(iprimary)=1.- mass_secondaries
        pmass1=1./max(pmass,tini)
        totmass=1.;totmass1=1.
!
        absolute_offset_star = absolute_offset_star*totmass1
        if (mass_secondaries >= 1.0) &
            call fatal_error('init_pointmasses', &
            'The mass of one '//&
            '(or more) of the particles is too big! The masses should '//&
            'never be bigger than g0. Please scale your assemble so that '//&
            'the combined mass of the (n-1) particles is less than that. '//&
            'The mass of the last particle in the pmass array will be '//&
            'reassigned to ensure that the total mass is g0')
!
!  Correct the semimajor of the secondaries by the offset they generate. 
!
        do ks=1,nqpar
          ! sign(A,B) returns the value of A with the sign of B
          if (ks/=iprimary) &
              positions(ks,1)=sign(1.,positions(ks,1))* (sma(ks) - absolute_offset_star)
        enddo
!
!  The last one (star) fixes the CM at Rcm=zero
!
        if (lcartesian_coords) then
          !put the star opposite to the baricenter of planets
          positions(iprimary,1)=-sign(1.,baricenter_secondaries)*absolute_offset_star
        elseif (lcylindrical_coords) then
          !put the star in positive coordinates, with pi for azimuth
          positions(iprimary,1)=absolute_offset_star
          positions(iprimary,2)=pi
        elseif (lspherical_coords) then
          positions(iprimary,1)=absolute_offset_star
          positions(iprimary,3)=pi
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
!  Here we substitute the first nqpar dust particles by massive ones,
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
!  Define isecondary. iprimary=1 and isecondary=2 is default
!
        if (iprimary == 2) isecondary=1
!
!  Reassign total mass of the star so that totmass=1
!
        pmass(iprimary)=1-pmass(isecondary)
        pmass1=1./max(pmass,tini)
        totmass=1.;totmass1=1.
!
!  Radial position at barycentric coordinates. Start both at apocenter,
!
!     r_i=(1+e)*a_i, where a_i = sma * m_j /(mi+mj)
!
!  See, i.e., Murray & Dermott, p.45, barycentric orbits.
!
        positions(isecondary,1)=(1+eccentricity) * semimajor_axis * pmass(  iprimary)/totmass
        positions(  iprimary,1)=(1+eccentricity) * semimajor_axis * pmass(isecondary)/totmass
!
!  Azimuthal position. Planet and star phased by pi.
!
        positions(isecondary,2)=0
        positions(  iprimary,2)=pi
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
        velocity_baricenter_secondaries=0.
        do ks=1,nqpar
          if (ks/=iprimary) then
            kep_vel(ks)=sqrt(1./sma(ks)) !circular velocity
            velocity_baricenter_secondaries = velocity_baricenter_secondaries + kep_vel(ks)*pmass(ks)
          endif
        enddo
        velocity_baricenter_secondaries = velocity_baricenter_secondaries*totmass
        do ks=1,nqpar
          if (ks/=iprimary) then
            if (lcartesian_coords) then
              velocity(ks,2) = sign(1.,positions(ks,1))*(kep_vel(ks) - velocity_baricenter_secondaries)
            elseif (lcylindrical_coords) then
              !positive for the planets
              velocity(ks,2) = kep_vel(ks) - velocity_baricenter_secondaries
            elseif (lspherical_coords) then
              velocity(ks,3) = kep_vel(ks) - velocity_baricenter_secondaries
            endif
          endif
       enddo
!       
        if (lcorrect_gasgravity_lstart) call initcond_correct_gasgravity(f,velocity,isecondary)
!
!  The last one (star) fixes the CM also with velocity zero
!
        if (lcartesian_coords) then
          velocity(iprimary,2)= -sign(1.,baricenter_secondaries)*velocity_baricenter_secondaries
        elseif (lcylindrical_coords) then
          velocity(iprimary,2)= velocity_baricenter_secondaries
        elseif (lspherical_coords) then
          velocity(iprimary,3)= velocity_baricenter_secondaries
        endif
!
!  Revert all velocities if retrograde
!
        if (lretrograde) velocity=-velocity
!
!  Allocate the point mass particles
!
        do k=1,nqpar
          fq(k,ivxq:ivzq) = velocity(k,1:3)
        enddo
!
      case ('eccentric')
!
!  Coded only for 2 bodies
!
        if (nqpar /= 2) call fatal_error("init_pointmasses",&
             "This initial condition is currently coded for 2 massive particles only.")
!
!  Define isecondary. iprimary=1 and isecondary=2 is default.
!
        if (iprimary == 2) isecondary=1
        velocity(isecondary,2) = sqrt((1-eccentricity)/(1+eccentricity) * GNewton/semimajor_axis) * pmass(  iprimary)/totmass
!
!  Correct secondary by gas gravity 
!
        if (lcorrect_gasgravity_lstart) call initcond_correct_gasgravity(f,velocity,isecondary)
!
        velocity(  iprimary,2) = velocity(isecondary,2) * pmass(isecondary)/pmass(iprimary)
!
!  Revert all velocities if retrograde.
!
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
    subroutine dvvq_dt_pointmasses_pencil(f,df,p)
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
      use Mpicomm, only: mpireduce_sum
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx,nqpar) :: rp_mn, rpcyl_mn
      real, dimension (mx,3) :: ggt
      real, dimension (3) :: xxq,rpsecondary
      real, dimension (nx) :: pot_energy,torque
      real, dimension (nx,3) :: accg_mn
      real :: accg_local,accg
      integer :: ks,j,ju
      logical :: lintegrate, lparticle_out
      logical :: ldiagnostic_only
!
      intent (in) :: f, p
      intent (inout) :: df
!
!  Get the pre-calculated gravity field.
!
      lhydroif: if (lhydro) then
        call get_total_gravity(ggt)
        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + ggt(l1:l2,:)
!
!  Add the gas gravity to the pointmasses by integration if llive_secondary
!  is true.   
!
!  The integration is needed for these two cases:
!
!   1. We are not solving the Poisson equation
!   2. We are, but a particle is out of the box (a star, for instance)
!      and therefore the potential cannot be interpolated.
!
        ldiagnostic_only=ldiagnos.and.idiag_torque(isecondary)/=0
!
        if (llive_secondary.or.ldiagnostic_only) then
          pointmasses1: do ks=1,nqpar
            lparticle_out=.false.
            if ((fq(ks,ixq)< xyz0(1)).or.(fq(ks,ixq) > xyz1(1)) .or. &
                 (fq(ks,iyq)< xyz0(2)).or.(fq(ks,iyq) > xyz1(2)) .or. &
                 (fq(ks,izq)< xyz0(3)).or.(fq(ks,izq) > xyz1(3))) then
               !particle out of box
               lparticle_out=.true.
             endif
             lintegrate=(.not.lselfgravity).or.lparticle_out.or.ldiagnostic_only
!
!  Sometimes making the star feel the selfgravity of the disk leads to
!  numerical troubles as the star is too close to the origin (in cylindrical
!  coordinates).
!
             if ((ks==iprimary).and.lnoselfgrav_primary) lintegrate=.false.
!
             integrategas: if (lintegrate) then
!
!  Get the acceleration particle ks suffers due to self-gravity.
!
               call get_radial_distance(rp_mn(:,ks),rpcyl_mn(:,ks),&
                    E1_=fq(ks,ixq),E2_=fq(ks,iyq),E3_=fq(ks,izq))
               xxq = fq(ks,ixq:izq)
!
               if (lcylindrical_gravity_nbody(ks)) then
                 call integrate_gasgravity(p,rpcyl_mn(:,ks),&
                      xxq,accg_mn,r_smooth(ks))
               else
                 call integrate_gasgravity(p,rp_mn(:,ks),&
                      xxq,accg_mn,r_smooth(ks))
               endif
!
!  Add it to its dfp
!
               if (llive_secondary) then
                 do j=1,3
                   accg_local =   sum(accg_mn(:,j))
                   call mpireduce_sum(accg_local,accg,1)
                   if (lroot) then
                     ju=j-1+ivxq
                     dfq(ks,ju) = dfq(ks,ju) + accg
                   endif
                 enddo
               endif
!
!  Calculate torques for output, if needed
!
               if (ldiagnos.and.idiag_torque(ks)/=0) then
                  rpsecondary=(/fq(ks,ixq)*cos(fq(ks,iyq)-y(m)),&
                       fq(ks,ixq)*sin(fq(ks,iyq)-y(m)),&
                       fq(ks,izq)                      &
                       /)
!
!  Integrate will add the cell volume, so we first remove it from the
!  calculation.
!
                  torque=(rpsecondary(1)*accg_mn(:,2)-rpsecondary(2)*accg_mn(:,1))*&
                       dVol1_x(l1:l2)*dVol1_y(m)*dVol1_z(n)
                  !call cross(rpsecondary,accg,torque)
                  call integrate_mn_name(pmass(ks)*torque,idiag_torque(ks))
               endif
!
             endif integrategas
           enddo pointmasses1
         endif
!
!  Diagnostic
!
        diagnos: if (ldiagnos) then
          pointmasses2: do ks=1,nqpar
!
            if (idiag_totenergy/=0.or.&
                idiag_torqext(ks)/=0.or.&
                idiag_torqint(ks)/=0) &
                call get_radial_distance(rp_mn(:,ks),rpcyl_mn(:,ks),&
                   E1_=fq(ks,ixq),E2_=fq(ks,iyq),E3_=fq(ks,izq))
!
!  Total energy
!
            if (idiag_totenergy/=0) then 
              pot_energy=0.0
              !potential energy
              pot_energy = pot_energy - &
                   GNewton*pmass(ks)*(rpcyl_mn(:,ks)**2+r_smooth(ks)**2)**(-0.5)
              if (ks==nqpar) then
                if (lcartesian_coords) then 
                  call sum_lim_mn_name(.5*p%rho*p%u2 + pot_energy,idiag_totenergy,p)
                else
                  call integrate_mn_name(.5*p%rho*p%u2 + pot_energy,idiag_totenergy)
                endif
              endif
            endif
!
!  Calculate torques splitting inner and outer, for backward compatibility
!
            if ((idiag_torqext(ks)/=0).or.(idiag_torqint(ks)/=0)) &
                 call calc_torque_split_int_ext(p,rpcyl_mn(:,ks),ks)
!
         enddo pointmasses2
        endif diagnos
      endif lhydroif
!
      call keep_compiler_quiet(f)
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
      logical :: ldiagnostic_only
!
      if (lroot) then
        call calc_hill_radius(hill_radius_square)
!        
        call dxxq_dt_pointmasses
        call dvvq_dt_pointmasses(hill_radius_square)
!
      endif
!
!  Add gas selfgravity if present
!
      if (lselfgravity) then
        ldiagnostic_only=ldiagnos.and.idiag_torque(isecondary)/=0
        if (llive_secondary.or.ldiagnostic_only) &
             call selfgravity_gas_on_pointmass
      endif
!
      if (lparticles) then
        call mpibcast_real(hill_radius_square,nqpar)
        call dvvp_dt_dustparticles(hill_radius_square)
      endif
!
      call keep_compiler_quiet(f,df)
!    
    endsubroutine pointmasses_pde
!***********************************************************************
    subroutine calc_hill_radius(hill_radius_square)

      real, dimension(nqpar), intent(out) :: hill_radius_square
      integer :: ks
      real :: rr,w2,sma2
!
      do ks=1,nqpar
        if (laccretion(ks).and.(ks/=iprimary)) then
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
          hill_radius_square(ks)=sma2*(pmass(ks)*pmass1(iprimary)/3)**(2./3.)
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
      if (lheader) print*, 'dvvq_dt_pointmasses: Calculate dvvq_dt_pointmasses'
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
        call gravity_pointmasses(k,hill_radius_square,.true.)
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
            if (idiag_period(ks)/=0) call point_par_name(2*pi*fq(ks,ixq)/fq(ks,ivyq),idiag_period(ks))
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

      use Particles_main, only: fetch_nparloc,fetch_fp_array,return_fp_array

      real, dimension (mpar_loc,mparray) :: fp_aux
      real, dimension (mpar_loc,mpvar) :: dfp_aux
      integer :: np_aux,k
      logical, dimension(mpar_loc) :: flag
      real, dimension(nqpar) :: hill_radius_square
!
      call fetch_nparloc(np_aux)
      if (np_aux/=0) then
        call fetch_fp_array(fp_aux,dfp_aux,&
             index%ixw,index%iyw,index%izw,&
             index%ivxw,index%ivyw,index%ivzw)
        do k=1,np_aux
          flag(k)=.false.  
          call gravity_pointmasses(k,hill_radius_square,.false.,&
               fp_aux(k,:),dfp_aux(k,:),flag(k))
        enddo
        call return_fp_array(fp_aux,dfp_aux,flag)
      endif
!
    endsubroutine dvvp_dt_dustparticles
!************************************************************
    subroutine gravity_pointmasses(k,hill_radius_square,&
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
      real :: rr2, r2_ij, rs2, v_ij, rsmooth2, Omega2_pm, rhill1, invr3_ij
      integer :: ks
!
      real, dimension (3) :: evr_cart,evr,positions
!
      logical, intent(in) :: lcallpointmass
!
      real, dimension (mparray), optional :: fp_pt
      real, dimension (mpvar), optional :: dfp_pt
      logical, optional :: flag_pt
!
      if (lcallpointmass) then
        positions    =  fq(k,ixq:izq)
      else
        positions    =  fp_pt(index%ixw:index%izw)
      endif
!
      do ks=1,nqpar
        if (.not.(lcallpointmass.and.k==ks)) then !prevent self-acceleration
!
!  r_ij = sqrt(ev1**2 + ev2**2 + ev3**2)
!
          call get_evr(positions,fq(ks,ixq:izq),evr_cart)!,.true.)
          rr2=sum(evr_cart**2)
          rsmooth2=r_smooth(ks)**2
!
!  Particles relative distance from each other:
!
          select case (ipotential_pointmass(ks))
!
          case ('newton-hill','newton','newtonian')
!
!  Newtonian potential. Should be used for a particle outside the grid (that does not need smoothing).
!
            r2_ij=max(rr2,rsmooth2)
            if (r2_ij > 0) then
              invr3_ij = r2_ij**(-1.5)
            else                ! can happen during pencil_check
              invr3_ij = 0.0
            endif
            Omega2_pm = GNewton*pmass(ks)*invr3_ij
!
          case ('plummer')
!
!  Potential of a Plummer sphere. It's better to use Boley's (below) but here for testing purposes.
!
            r2_ij=rr2+rsmooth2
            Omega2_pm = GNewton*pmass(ks)*r2_ij**(-1.5)
!
          case ('boley')
!
!  Potential that Aaron Boley uses (should ask him for reference). Exactly matches Newtonian outside Hill radius. 
!
            r2_ij=rr2
            if (r2_ij .gt. hill_radius_square(ks)) then
              Omega2_pm =  GNewton*pmass(ks)*r2_ij**(-1.5)
            else
              rhill1=1./sqrt(hill_radius_square(ks))
              Omega2_pm = -GNewton*pmass(ks)*(3*sqrt(r2_ij)*rhill1 - 4)*rhill1**3
            endif
!

          case default
            if (lroot) print*, 'gravity_pointmasses: '//&
                 'No such value for ipotential_pointmass: ', trim(ipotential_pointmass(ks))
            call fatal_error("","")
          endselect  
!
!  If there is accretion, remove the accreted particles from the simulation, if any.
!
          if (.not.(lcallpointmass).and.laccretion(ks)) then
            rs2=(accrete_hills_frac(ks)**2)*hill_radius_square(ks)
            if (r2_ij<=rs2) then
              !flag particle for removal 
              flag_pt=.true.
              return
            endif
          endif
!
!  Gravitational acceleration: g=-g0/|r-r0|^3 (r-r0)
!
!  The acceleration is in non-coordinate basis (all have dimension of length).
!  The main dxx_dt of particle_dust takes care of transforming the linear
!  velocities to angular changes in position.
!
          if (lcallpointmass) then 
            dfq_cart(k,:) = dfq_cart(k,:) - Omega2_pm*evr_cart(1:3)
          else
            dfp_pt(ivpx_cart:ivpz_cart) = &
                 dfp_pt(ivpx_cart:ivpz_cart) - Omega2_pm*evr_cart(1:3)
          endif
!
          if (ladd_dragforce) call dragforce_pointmasses(k)
!
          if (lcoriolis_force) then
            dfq(k,ivxq) = dfq(k,ivxq) + 2*Omega_Coriolis*fq(k,ivyq)
            dfq(k,ivyq) = dfq(k,ivyq) - 2*Omega_Coriolis*fq(k,ivxq)
          endif
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
                 sqrt(Omega2_pm)/cdtq)
          endif
!
        endif !if (ipar(k)/=ks)
!
     enddo !nbody loop
!
    endsubroutine gravity_pointmasses
!**********************************************************
    subroutine dragforce_pointmasses(k)
!
!  Adds dragforce on massive particles.
!  06-feb-18/wlad: coded
!
      real, dimension (3) :: uup
      integer, intent(in) :: k
!
!  Supports only Cartesian with ugas=uy so far.
!
      uup=(/0.,ugas,0./)
      if (llinear_drag) then
        dfq(k,ivxq:ivzq) = dfq(k,ivxq:ivzq) - (fq(k,ivxq:ivzq)-uup)/StokesNumber(k)
      else if (lquadratic_drag) then
        dfq(k,ivxq:ivzq) = dfq(k,ivxq:ivzq) - &
             abs(fq(k,ivxq:ivzq)-uup)*(fq(k,ivxq:ivzq)-uup)/StokesNumber(k)
      else
        call fatal_error("drag should be linear or quadratic","")
      endif
!
    endsubroutine dragforce_pointmasses
!**********************************************************
    subroutine selfgravity_gas_on_pointmass()
!
!  Adds selfgravity on massive particles. This subroutine
!  fetchs the (global) acceleration array calculated by the
!  logspirals method and interpolates the acceleration to
!  the position of the pointmass. TODO: other poisson methods
!  should get the potential, take the gradient and interpolate.
!
!  14-may-20/wlad: coded
!
      use Mpicomm
      use Sub, only: get_radial_distance,cross
      use Poisson,     only:inverse_laplacian,get_acceleration
!
      integer :: ks
      real, dimension (nx,ny,nz,3) :: acceleration
      real, dimension (3) :: accg,torque
!
      call get_acceleration(acceleration)
!
      do ks=1,nqpar
        if (ks/=iprimary) then
          call bilinear_interpolate(acceleration,fq(ks,ixq:izq),accg)
          if (lroot.and.llive_secondary) &
               dfq(ks,ivxq:ivzq) = dfq(ks,ivxq:ivzq) + accg
!
          if (ldiagnos) call calc_torque(ks,accg)
        endif
      enddo
!
    endsubroutine selfgravity_gas_on_pointmass
!**********************************************************
    subroutine get_evr(xxp,xxq,evr_output)
!
!  Point-to-point vector distance, in different coordinate systems.
!  Return always in Cartesian.
!
!  14-feb-14/wlad: coded
!
      real, dimension(3), intent(in) :: xxp,xxq
      real, dimension(3), intent(out) :: evr_output
      real :: x1,y1,x2,y2,z1,z2
      real :: e1,e2,e3,e10,e20,e30
!
      e1=xxp(1);e10=xxq(1)
      e2=xxp(2);e20=xxq(2)
      e3=xxp(3);e30=xxq(3)
!
      if (lcartesian_coords) then
        x1=e1 ; x2=e10
        y1=e2 ; y2=e20
        z1=e3 ; z2=e30
      elseif (lcylindrical_coords) then
!
        x1=e1*cos(e2)
        y1=e1*sin(e2)
        z1=e3
!
        x2=e10*cos(e20)
        y2=e10*sin(e20)
        z2=e30
!
      elseif (lspherical_coords) then
        x1=e1*sin(e2)*cos(e3)
        y1=e1*sin(e2)*sin(e3)
        z1=e1*cos(e2)
!
        x2=e10*sin(e20)*cos(e30)
        y2=e10*sin(e20)*sin(e30)
        z2=e10*cos(e20)
      endif
!
      evr_output(1)=x1-x2
      evr_output(2)=y1-y2
      evr_output(3)=z1-z2
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
    subroutine calc_torque(ks,accg)
!
      use Sub, only: cross
!
      integer, intent(in) :: ks
      real, dimension(3), intent(in) :: accg
      real, dimension(3) :: torque
!
      if (lroot) then
        if (idiag_torque(ks)/=0) then
          call cross((/fq(ks,ixq),0.,0./),accg,torque)
          call point_par_name(pmass(ks)*torque(3),idiag_torque(ks))
        endif
      endif
!
    endsubroutine calc_torque
!***********************************************************************             
    subroutine calc_torque_split_int_ext(p,dist,ks)
!
!  Output torque diagnostic for nbody particle ks
!  (maintained for backward compatibility). As of May 2020
!  the torque diagnostic should be computed after integrate_gasgravity,
!  from the same acceleration used on the planet, instead of integrating
!  separately.
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
      if (ks==iprimary) call fatal_error('calc_torque', &
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
      hills = smap*(pmass(ks)*pmass1(iprimary)/3.)**(1./3.)
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
        if (lexclude_hills) then 
          do i=1,nx
            if (dist(i)<hills) then
              torque_gas(i)=0.
              torque_par(i)=0.
            endif
          enddo
        endif
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
    endsubroutine calc_torque_split_int_ext
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
          if (ks/=iprimary) then
            pmass(ks)= max(dble(tini),&
                 final_ramped_mass(ks)*(sin((.5*pi)*(t/ramping_period))**2))
            tmp=tmp+pmass(ks)
          endif
        enddo
        pmass(iprimary)= 1-tmp
      else
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
!
      real, dimension (mx)       :: Omega2_pm,rrp
      real                       :: rr,rp1,rhill,rhill1
      integer                    :: ks,i
      !real, dimension (mx)       :: grav_particle,rrp
      !integer                    :: ks
!
      intent(out) :: ggt
!
      ggt=0.
      do ks=1,nqpar
!
!  Hill radius
!
        if (lcartesian_coords) then 
          rp1 = sqrt(fq(ks,ixq)**2+fq(ks,iyq)**2+fq(ks,izq)**2)
        elseif (lcylindrical_coords) then
          rp1 = sqrt(fq(ks,ixq)**2+              fq(ks,izq)**2)
        elseif (lspherical_coords) then
          rp1 =      fq(ks,ixq)
        endif
!
        rhill  = rp1*(GNewton*pmass(ks)/3.)**(1./3)
        rhill1 = 1./rhill
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
        select case (ipotential_pointmass(ks))
!
        case ('plummer')
!
!  Potential of a Plummer sphere
!
          if (ks==iprimary) call fatal_error("get_total_gravity",&
               "The primary can only be newtonian, please switch ipotential_pointmass")
          Omega2_pm =-GNewton*pmass(ks)*(rrp**2+r_smooth(ks)**2)**(-1.5)
!
        case ('boley')
!
!  Correct potential outside Hill sphere
!
          if (ks==iprimary) call fatal_error("get_total_gravity",&
               "The primary can only be newtonian, please switch ipotential_pointmass")
          do i=1,mx
            if (rrp(i) .gt. rhill) then
              Omega2_pm(i) = -GNewton*pmass(ks)*rrp(i)**(-3)
            else
              Omega2_pm(i) =  GNewton*pmass(ks)*(3*rrp(i)*rhill1 - 4)*rhill1**3
            endif
          enddo
!
        case ('newton-hill','newton','newtonian')
!
!  Newtonian potential; same as boley but constant inside rsmooth
!
          if (ks==iprimary.and.r_smooth(ks)/=0) call fatal_error("get_total_gravity",&
               "Use r_smooth=0 for the primary's potential")
!
          do i=1,mx
            rr=max(rrp(i),r_smooth(ks))
            if (rr > 0) then
              Omega2_pm(i) = -GNewton*pmass(ks)*rr**(-3)
            else                ! can happen during pencil_check
              Omega2_pm(i) = 0.
            endif
          enddo
!                                                                                                                                     
        case default
!                                                                                                                                     
!  Catch unknown values                                                                                                               
!                                                                                                                                     
          if (lroot) print*, 'get_total_gravity: '//&
               'No such value for ipotential_pointmass: ', trim(ipotential_pointmass(ks))
          call fatal_error("","")
        endselect
!
        call get_gravity_field_pointmasses(Omega2_pm,ggp,ks)
!
        if ((ks==iprimary).and.lnogravz_star) &
            ggp(:,3) = 0.
!
!  Sum up the accelerations of the massive particles
!
        ggt=ggt+ggp
!
!  Add indirect term if the star is fixed at the center
!
      enddo
!
    endsubroutine get_total_gravity
!***********************************************************************
    subroutine initcond_correct_gasgravity(f,velocity,k)
!
!  Calculates acceleration on the point (x,y,z)=xxpar
!  due to the gravity of gas. 
!
!  14-may-20/wlad : coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(nqpar,3), intent(inout) :: velocity
      real, dimension(3) :: accg
      real :: vphi,phidot2,OO,rr
      integer, intent(in) :: k
!
      if (lselfgravity) then
        call correct_gasgravity_selfgravity(f,k,accg)
      else
        call correct_gasgravity_integrate(f,k,accg)
      endif
!
!  Correct original velocity by this acceleration
!
!  phidot**2 = OmegaK**2 + 1/r (d/dr PHI_sg)
!  phidot**2 = OmegaK**2 - 1/r gr
!
!  Current azimuthal velocity and azimuthal frequency
!
      if (lcylindrical_coords) then
        vphi = velocity(k,2)
        rr = fq(k,ixq)
      elseif (lspherical_coords) then
        vphi = velocity(k,3)
        rr = fq(k,ixq)*sin(fq(k,iyq))
      endif
      OO = vphi/rr
!
!  Update angular velocity by selfgravitational acceleration
!  This line assumes axisymmetry
!
      phidot2 = OO**2 - accg(1)/rr
!
!  Corrected velocity
!
      vphi = sqrt(phidot2) * rr
!
      if (lcylindrical_coords) then
        velocity(k,2) = vphi
      elseif (lspherical_coords) then
        velocity(k,3) = vphi
      endif      
!
    endsubroutine initcond_correct_gasgravity
!***********************************************************************
    subroutine correct_gasgravity_selfgravity(f,k,accg)
!
!  Calculates acceleration on the point (x,y,z)=xxpar
!  due to the gravity of gas, by interpolating from the
!  gravitational acceleration from the selfgravity module. 
!
!  14-may-20/wlad : coded
!
      use Mpicomm
      use Sub, only: get_radial_distance
      use Poisson,     only:inverse_laplacian,get_acceleration
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny,nz,3) :: acceleration
      real, dimension (nx,ny,nz) :: rho
      real, dimension (3) :: accg
      integer, intent (in) :: k
!
      if (ldensity_nolog) then
         rho=f(l1:l2,m1:m2,n1:n2,irho)
      else
         rho=exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
      endif
      call inverse_laplacian(rho)
      call get_acceleration(acceleration)
      call bilinear_interpolate(acceleration,fq(k,ixq:izq),accg)
!
    endsubroutine correct_gasgravity_selfgravity
!***********************************************************************
    subroutine correct_gasgravity_integrate(f,k,accg)
!
!  Calculates acceleration on the point (x,y,z)=xxpar
!  due to the gravity of gas, by integrating the whole grid. 
!
!  29-aug-18/wlad : coded
!
      use Mpicomm
      use Sub, only: get_radial_distance
!
      implicit none
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(3) :: xxpar,sum_loc,accg
      real, dimension(nx,3) :: dist
      real, dimension(nx) :: rrp,rp_mn,rpcyl_mn,selfgrav,density
      real, dimension(nx) :: dv,jac,dqy,tmp
      real :: vphi,phidot2,OO,rr
      real :: dqx,dqz,rp0,fac
      integer :: j,k
!
!  Sanity check
!
      !if (.not.(lgas_gravity.or.ldust_gravity)) &
      !     call fatal_error("initcond_correct_selfgravity",&
      !     "No gas gravity or dust gravity to add. "//&
      !     "Switch on lgas_gravity or ldust_gravity in n-body parameters")
!
      xxpar = fq(k,ixq:izq)
      rp0=r_smooth(k)
!      
      sum_loc=0.
!      
      mloop: do m=m1,m2
      nloop: do n=n1,n2
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
          call fatal_error('correct_gasgravity_integrate', &
               ' not yet implemented for spherical polars')
           dqx=0.;dqy=0.;dqz=0.
        else
          call fatal_error('correct_gasgravity_integrate','wrong coord_system')
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
!  selfgrav = - G*((rho+rhop)*dv)*mass*r*(r**2 + r0**2)**(-1.5)
!  gx = selfgrav * r\hat dot x\hat
!  -> gx = selfgrav * (x-x0)/r = - G*((rho+rhop)*dv)*mass*(r**2+r0**2)**(-1.5) * (x-x0)
!
        density=0.
        if (ldensity_nolog) then
           density=f(l1:l2,m,n,irho)
        else
           density=exp(f(l1:l2,m,n,ilnrho))
        endif
!
!  Add the particle gravity if npar>mspar (which means dust is being used)
!
        !if (ldust.and.ldust_gravity) density=density+f(l1:l2,m,n,irhop)
!
        call get_radial_distance(rp_mn,rpcyl_mn,E1_=xxpar(1),E2_=xxpar(2),E3_=xxpar(3))
        if (lcylindrical_gravity_nbody(k)) then
           rrp=rpcyl_mn
        else
           rrp=rp_mn
        endif
!
        selfgrav = -GNewton*density_scale*&
             density*jac*dv*(rrp**2 + rp0**2)**(-1.5)
!
!  Exclude the frozen zones
!
        if (lcartesian_coords) call fatal_error("","")
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
            sum_loc(j) = sum_loc(j)+fac*sum(tmp)
          else
            sum_loc(j) = sum_loc(j)+fac*(sum(tmp(2:nx-1))+.5*(tmp(1)+tmp(nx)))
          endif
        enddo
      enddo nloop
      enddo mloop
!
      do j=1,3
        call mpireduce_sum(sum_loc(j),accg(j))
      enddo
!
!  Broadcast particle acceleration
!
      call mpibcast_real(accg,3)
!
    endsubroutine correct_gasgravity_integrate
!***********************************************************************
    subroutine bilinear_interpolate(v,q,vp)
!
!  Interpolate the gravity vector with bilinear interpolation to the position
!  of the pointmass(es). 
!
      use Mpicomm, only: mpisend_real,mpirecv_real,mpibcast_real
!      
      real, dimension (nx,ny,nz,3) :: v
      real, dimension (3) :: q
      real, dimension (3) :: vp
      integer :: j

      real :: rdx1,rdy1,rdz1
      real :: rdx2,rdy2,rdz2
      integer :: ix0_global,iy0_global,iz0_global
      integer :: ix1_global,iy1_global,iz1_global

      integer :: ipx0,ipy0,ipz0
      integer :: ipx1,ipy1,ipz1
      integer :: ix0,iy0,iz0
      integer :: ix1,iy1,iz1
!
      integer :: iproc_q11,iproc_q21,iproc_q12,iproc_q22
      real, dimension(3) :: Q11,Q21,Q12,Q22
!
      intent(in) :: v,q
      intent(out) :: vp
!
      call calc_indices_and_processors(q(1),xgrid,dx1grid,nprocx,nx,&
           lperi(1),ix0,ix1,ipx0,ipx1,rdx1,rdx2)
      call calc_indices_and_processors(q(2),ygrid,dy1grid,nprocy,ny,&
           lperi(2),iy0,iy1,ipy0,ipy1,rdy1,rdy2)
!         
      if (nz==1) then
        iz0=1
      else
        call fatal_error("trilinear_interpolate","logpsirals works only for 2D")
        !call calc_indices_and_processors(q(3),zgrid,dz1grid,nprocz,nz,&
        !lperi(3),iz0,iz1,ipz0,ipz1,rdz1,rdz2)
      endif
!
      iproc_q11 = ipx0 + nprocx*ipy0 + nprocx*nprocy*ipz
      iproc_q21 = ipx1 + nprocx*ipy0 + nprocx*nprocy*ipz
      iproc_q12 = ipx0 + nprocx*ipy1 + nprocx*nprocy*ipz
      iproc_q22 = ipx1 + nprocx*ipy1 + nprocx*nprocy*ipz
!
      if (iproc_q11 /= root) then
        if (iproc==iproc_q11) call mpisend_real(v(ix0,iy0,iz0,:),3,root     ,111)
        if (lroot)            call mpirecv_real(Q11             ,3,iproc_q11,111)
      else
        Q11 = v(ix0,iy0,iz0,:)
      endif
!
      if (iproc_q21 /= root) then
        if (iproc==iproc_q21) call mpisend_real(v(ix1,iy0,iz0,:),3,root     ,121)
        if (lroot)            call mpirecv_real(Q21             ,3,iproc_q21,121)
      else
        Q21 = v(ix1,iy0,iz0,:)
      endif
!
      if (iproc_q12 /= root) then
        if (iproc==iproc_q12) call mpisend_real(v(ix0,iy1,iz0,:),3,root     ,112)
        if (lroot)            call mpirecv_real(Q12             ,3,iproc_q12,112)
      else
        Q12 = v(ix0,iy1,iz0,:)
      endif
!
      if (iproc_q22 /= root) then
        if (iproc==iproc_q22) call mpisend_real(v(ix1,iy1,iz0,:),3,root     ,122)
        if (lroot)            call mpirecv_real(Q22             ,3,iproc_q22,122)
      else
        Q22 = v(ix1,iy1,iz0,:)
      endif
!
      if (lroot) then 
        do j=1,3         
!
! Interpolation formula
!
          vp(j) &
               =   Q11(j) * rdx2* rdy2 &   ! Q11
                 + Q21(j) * rdx1* rdy2 &   ! Q21
                 + Q12(j) * rdx2* rdy1 &   ! Q12
                 + Q22(j) * rdx1* rdy1     ! Q22
        enddo
      endif
!
      call mpibcast_real(vp,3)
!
    endsubroutine bilinear_interpolate
!***********************************************************************
    subroutine calc_indices_and_processors(q,qgrid,dq1grid,nprocq,nq,lperiodic_direction,&
         iq0,iq1,ipq0,ipq1,rdq1,rdq2)
!
      use Sub, only: find_index_by_bisection
!
      real, dimension(:) :: qgrid
      real, dimension(:) :: dq1grid
      real :: q,rdq1,rdq2,dq1,dq
      integer :: iq0_global,iq1_global,ipq0,ipq1,iq0,iq1,nprocq,nq
      logical :: lperiodic_direction
      integer :: nqgrid
!
      intent(in) :: q,qgrid,lperiodic_direction,nprocq,nq
      intent(out) :: iq0,iq1,ipq0,ipq1,rdq1,rdq2

      lperiodic: if (.not.lperiodic_direction) then
!
!  Non-periodic grid          
!
        call find_index_by_bisection(q,qgrid,iq0_global)      
        if (qgrid(iq0_global)>q) iq0_global = iq0_global-1
        iq1_global=iq0_global+1
!      
        dq1=1./(qgrid(iq1_global)-qgrid(iq0_global))
        rdq1 = (q-qgrid(iq0_global))*dq1
        rdq2 = (qgrid(iq1_global)-q)*dq1
!
        if (nprocq/=1) then 
          ipq0 = (iq0_global-1)/nq
          ipq1 = (iq1_global-1)/nq
!
          iq0  = iq0_global-ipq0*nq
          iq1  = iq1_global-ipq1*nq
        else
          ipq0=0
          ipq1=0
!
          iq0=iq0_global
          iq1=iq1_global
        endif
!
      else
!
!  Periodic grid
!
        call find_index_by_bisection(q,qgrid,iq0_global)
        dq1=dq1grid(iq0_global)
        dq=1./dq1
        if (qgrid(iq0_global)>q) iq0_global = iq0_global-1
        iq1_global=iq0_global+1
!
        nqgrid=size(qgrid)
        if (iq0_global < 1) then
          iq0_global = iq0_global + nqgrid
          rdq1 = (q-(qgrid(1)-dq))*dq1
        else
          rdq1 = (q-qgrid(iq0_global))*dq1
        endif
!
        if (iq1_global > nqgrid) then
          iq1_global = iq1_global - nqgrid
          rdq2 = (qgrid(nqgrid)+dq-q)*dq1
        else
          rdq2 = (qgrid(iq1_global)-q)*dq1
        endif
!
        if (nprocq/=1) then 
          ipq0=(iq0_global-1)/nq
          ipq1=(iq1_global-1)/nq
!
          iq0=iq0_global-ipq0*nq
          iq1=iq1_global-ipq1*nq
        else
          ipq0=0
          ipq1=0
!
          iq0=iq0_global
          iq1=iq1_global
        endif
!
      endif lperiodic
!
    endsubroutine calc_indices_and_processors
!***********************************************************************    
    subroutine integrate_gasgravity(p,rrp,xxpar,accg,rp0)
!
!  Calculates acceleration on the point (x,y,z)=xxpar
!  due to the gravity of the gas+dust.
!
!  15-sep-06/wlad : coded
!
      use Mpicomm
!
      real, dimension(nx,3) :: dist,accg
      real, dimension(nx) :: rrp,rr,gasgravity,density,cellmass
      real :: rp0
      real, dimension(3) :: xxpar
      integer :: j
      type (pencil_case) :: p
      logical :: lfirstcall=.true.
!
      intent(out) :: accg
!
      if (lfirstcall.and.lroot) &
          print*,'Adding gas+dust gravity to the massive particles'
!
!  Sanity check
!
      if (.not.(lgas_gravity.or.ldust_gravity)) &
           call fatal_error("lintegrate_gasgravity",&
           "No gas gravity or dust gravity to add. "//&
           "Switch on lgas_gravity or ldust_gravity in n-body parameters")
!
      if (coord_system=='cartesian') then
        dist(:,1)=x(l1:l2)-xxpar(1)
        dist(:,2)=y(  m  )-xxpar(2)
        dist(:,3)=z(  n  )-xxpar(3)
      elseif (coord_system=='cylindric') then
        !re = r-rp
        !rp = rp*cos(phip-phi) \hat{r} + rp*sin(phi-phi) \hat{phi}
        dist(:,1)=x(l1:l2)-xxpar(1)*cos(y(m)-xxpar(2))
        dist(:,2)=         xxpar(1)*sin(y(m)-xxpar(2))
        dist(:,3)=z(  n  )-xxpar(3)
      elseif (coord_system=='spherical') then
        call fatal_error('integrate_gasgravity', &
            ' not yet implemented for spherical polars')
      else
        call fatal_error('integrate_gasgravity','wrong coord_system')
      endif
!
!  The gravity of every single cell - should exclude inner and outer radii...
!
!  selfgrav = G*((rho+rhop)*dv)*mass*r*(r**2 + r0**2)**(-1.5)
!  gx = selfgrav * r\hat dot x\hat
!  -> gx = selfgrav * (x-x0)/r = G*((rho+rhop)*dv)*mass*(r**2+r0**2)**(-1.5) * (x-x0)
!
      density=0.
      if (lgas_gravity) density=p%rho
!
!  Add the particle gravity if npar>mspar (which means dust is being used)
!
      if (ldust.and.ldust_gravity) density=density+p%rhop
!
      cellmass=density*dVol_x(l1:l2)*dVol_y(m)*dVol_z(n)
      gasgravity = GNewton*density_scale*cellmass*(rrp**2 + rp0**2)**(-1.5)
!
!  Exclude the frozen zones
!
      if (lexclude_frozen.and.lcylinder_in_a_box) then
       where ((p%rcyl_mn<=r_int).or.(p%rcyl_mn>=r_ext))
          gasgravity = 0
        endwhere
      else
        if (l2D.or.(l3D.and.lcylindrical_gravity)) then
          rr=p%rcyl_mn
        else
          rr=p%r_mn
        endif
        where ((rr<=r_int).or.(rr>=r_ext))
          gasgravity = 0
        endwhere
      endif
!
!  Integrate the accelerations on this processor
!  And sum over processors with mpireduce
!
      do j=1,3
        accg(:,j)=gasgravity*dist(:,j)
      enddo
!
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine integrate_gasgravity
!***********************************************************************
    subroutine pointmasses_read_snapshot(filename)
!
!  Read nbody particle info
!
!  01-apr-08/wlad: dummy
!
      use IO, only: input_pointmass
!
      character (len=*), intent(in) :: filename
!
      call input_pointmass(filename, qvarname, fq, nqpar, mqarray)
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
    subroutine pointmasses_write_snapshot(file,enum,flist)
!
      use General, only:safe_character_assign
      use Sub, only: update_snaptime, read_snaptime
      use IO, only: log_filename_to_file, output_pointmass
!
!  Input and output of information about the massive particles
!
!  01-apr-08/wlad: coded
!
      character (len=*), intent(in) :: file
      logical, intent(in) :: enum
      character (len=*), optional, intent(in) :: flist
!
      logical, save :: lfirst_call=.true.
      integer, save :: nsnap
      real, save :: tsnap
      character (len=fnlen) :: filename, filename_diag
      logical :: lsnap
      character (len=intlen) :: nsnap_ch
!
!  Input is either file=qvar if enum=F or file=QVAR if enum=T.
!  If the latter, append a number to QVAR.
!      
      filename = trim(file)
      if (enum) then
!
!  Prepare to read tsnap.dat         
!
        call safe_character_assign(filename_diag,trim(datadir)//'/tsnap.dat')
!        
!  Read the tsnap.dat to get N for QVARN (N=nsnap).
!
        if (lfirst_call) then
          call read_snaptime(filename_diag,tsnap,nsnap,dsnap,t)
          lfirst_call=.false.
        endif
!
!  Update snaptime to set nsnap_ch=N and append N to filename QVARN.
!
        call update_snaptime(filename_diag,tsnap,nsnap,dsnap,t,lsnap,nsnap_ch,nowrite=.true.)
        call safe_character_assign(filename,trim(filename)//trim(nsnap_ch))
      endif
!
      if (lsnap .or. .not. enum ) then
        call output_pointmass (filename, qvarname, fq, nqpar, mqarray)
        if (lroot .and. present(flist)) call log_filename_to_file(filename,flist)
      endif
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
      use HDF5_IO, only: pointmass_index_append
      use General, only: itoa
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname,ks,j
      character :: str
      character (len=intlen) :: sks
      logical :: lwr
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset
!
      if (lreset) then
        idiag_xxq=0;idiag_vvq=0
        idiag_torqint=0;idiag_torqext=0
        idiag_totenergy=0
        idiag_torqext_gas=0;idiag_torqext_par=0
        idiag_torqint_gas=0;idiag_torqint_par=0
        idiag_period=0
        idiag_torque=0
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
            call pointmass_index_append('i_'//trim(str)//'q'//trim(sks),idiag_xxq(ks,j))
            call pointmass_index_append('i_v'//trim(str)//'q'//trim(sks),idiag_vvq(ks,j))
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
          call parse_name(iname,cname(iname),cform(iname),&
               'period'//trim(sks),idiag_period(ks))
          call parse_name(iname,cname(iname),cform(iname),&
               'torque_'//trim(sks),idiag_torque(ks))
        enddo
!
        if (lwr) then
          call pointmass_index_append('i_torqint_'//trim(sks),idiag_torqint(ks))
          call pointmass_index_append('i_torqext_'//trim(sks),idiag_torqext(ks))
          call pointmass_index_append('i_torqint_gas'//trim(sks),idiag_torqint(ks))
          call pointmass_index_append('i_torqext_gas'//trim(sks),idiag_torqext(ks))
          call pointmass_index_append('i_torqint_par'//trim(sks),idiag_torqint(ks))
          call pointmass_index_append('i_torqext_par'//trim(sks),idiag_torqext(ks))
          call pointmass_index_append('i_period'//trim(sks),idiag_period(ks))
          call pointmass_index_append('i_torque_'//trim(sks),idiag_torque(ks))
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
         call pointmass_index_append('i_totenergy',idiag_totenergy)
       endif
!
    endsubroutine rprint_pointmasses
!***********************************************************************
    subroutine boundconds_pointmasses
!
!  Global boundary conditions for particles.
!
!  30-dec-04/anders: coded
!
      use Mpicomm
      use General, only: random_number_wrapper
      use SharedVariables, only: get_shared_variable
!
      integer :: k
      character (len=2*bclen+1) :: boundx, boundy, boundz
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
          elseif (boundy=='p2pi') then
!  yp < y0
            if (fq(k,iyq)< -pi) then
              fq(k,iyq)=fq(k,iyq)+2*pi
              if (fq(k,iyq)<-pi) then
                print*, 'boundconds_pointmasses: ERROR - particle ', k, &
                     ' was further than Ly outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'xxq=', fq(k,ixq:izq)
                call fatal_error_local('boundconds_pointmasses','')
              endif
            endif
!  yp > y1
            if (fq(k,iyq)>=pi) then
              fq(k,iyq)=fq(k,iyq)-2*pi
              if (fq(k,iyq)>=pi) then
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
      use Particles_main, only: fetch_nparloc
!    
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: np_aux
!
      if (lroot) then
        if (itsub==1) then
          dfq(1:nqpar,:) = 0.0
          dfq_cart(1:nqpar,:)=0.0
        else
          dfq(1:nqpar,:)=alpha_ts(itsub)*dfq(1:nqpar,:)
          dfq_cart(1:nqpar,:)=alpha_ts(itsub)*dfq_cart(1:nqpar,:)
        endif
      endif
!
      call keep_compiler_quiet(f)
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
      use Particles_main, only: fetch_nparloc,fetch_fp_array,return_fp_array
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(3) :: xx_polar,vv_polar,xxdot_polar,vvdot_polar,aa_nbody_cart
      integer :: ks
!
      real, dimension (mpar_loc,mparray) :: fp_aux
      real, dimension (mpar_loc,mpvar) :: dfp_aux
      integer :: np_aux,k
!
      if (lroot) then 
        do ks=1,nqpar
          xx_polar    = fq(ks,ixq:izq)
          vv_polar    = fq(ks,ivxq:ivzq)
          xxdot_polar = dfq(ks,ixq:izq)
          vvdot_polar = dfq(ks,ivxq:ivzq)
!            
          aa_nbody_cart=dfq_cart(ks,:)
          call advance_particles_in_cartesian(xx_polar,vv_polar,xxdot_polar,vvdot_polar,aa_nbody_cart)
          fq(ks,ixq:izq)=xx_polar
          fq(ks,ivxq:ivzq)=vv_polar
          dfq(ks,ixq:izq)=xxdot_polar
        enddo
      endif
      call mpibcast_real(fq,(/nqpar,mqarray/))
!
      if (lparticles) then
        call fetch_nparloc(np_aux)
        if (np_aux/=0) then
          call fetch_fp_array(fp_aux,dfp_aux,&
               index%ixw,index%iyw,index%izw,&
               index%ivxw,index%ivyw,index%ivzw)
          do k=1,np_aux
            xx_polar    = fp_aux(k,index%ixw:index%izw)
            vv_polar    = fp_aux(k,index%ivxw:index%ivzw)
            xxdot_polar = dfp_aux(k,index%ixw:index%izw)
            vvdot_polar = dfp_aux(k,index%ivxw:index%ivzw)
            aa_nbody_cart=dfp_aux(k,ivpx_cart:ivpz_cart)
!            
            call advance_particles_in_cartesian(xx_polar,vv_polar,xxdot_polar,vvdot_polar,aa_nbody_cart)
!
            fp_aux(k,index%ixw:index%izw)=xx_polar
            fp_aux(k,index%ivxw:index%ivzw)=vv_polar
            dfp_aux(k,index%ixw:index%izw)=xxdot_polar
          enddo
          call return_fp_array(fp_aux,dfp_aux)
        endif
!
      endif
!
      call keep_compiler_quiet(f)
!      
    endsubroutine pointmasses_timestep_second
!***********************************************************************
    subroutine advance_particles_in_cartesian(xx_polar,vv_polar,xxdot_polar,vvdot_polar,aa_nbody_cart)
!
! With N-body gravity, the particles should have their position advanced in
! Cartesian coordinates, for better conservation of the Jacobi constant, even
! if the grid is polar.
!
! 14-feb-14/wlad: coded
!
      real :: rad,phi,tht,zed,xp,yp,zp
      real :: vrad,vtht,vphi,vzed,vx,vy,vz
      real :: raddot,thtdot,phidot,zeddot,xdot,ydot,zdot
      real :: vraddot,vthtdot,vphidot,vzeddot,vxdot,vydot,vzdot
      real :: cosp,sinp,cost,sint
!
      real, dimension(3), intent(inout) :: xx_polar,vv_polar
      real, dimension(3) :: aa_nbody_cart,vv_cart,dvv_cart,xx_cart,dxx_cart,xxdot_polar,vvdot_polar
!
      if (lcartesian_coords) then 
!
        xp = xx_polar(1); yp = xx_polar(2); zp = xx_polar(3)
        vx = vv_polar(1); vy = vv_polar(2); vz = vv_polar(3)

        xdot = xxdot_polar(1) + vx
        ydot = xxdot_polar(2) + vy
        zdot = xxdot_polar(3) + vz
!
        vxdot = vvdot_polar(1) + aa_nbody_cart(1) ! ivqx_cart
        vydot = vvdot_polar(2) + aa_nbody_cart(2) ! ivqy_cart
        vzdot = vvdot_polar(3) + aa_nbody_cart(3)
!
      else if (lcylindrical_coords) then
!
!  The input position, velocities, and accelerations in cylindrical coordinates.
!
        rad  = xx_polar(1) ; phi   = xx_polar(2) ; zed  = xx_polar(3)
        vrad = vv_polar(1) ; vphi  = vv_polar(2) ; vzed = vv_polar(3)

        raddot = xxdot_polar(1)  ; phidot = xxdot_polar(2)  ; zeddot = xxdot_polar(3)
        vraddot = vvdot_polar(1) ; vphidot = vvdot_polar(2) ; vzeddot = vvdot_polar(3)
!
!  Shortcuts.
!
        cosp = cos(phi) ; sinp = sin(phi)
!
!  Transform the positions and velocities to Cartesian coordinates.
!
        xp = rad*cosp ; yp = rad*sinp ; zp = zed
!
        vx = vrad*cosp - vphi*sinp
        vy = vrad*sinp + vphi*cosp
        vz = vzed
!
!  Transform the position and velocity derivatives to Cartesian coordinates.        
!  Add the cartesian velocity to position and cartesian acceleration to velocity.
!
        xdot = raddot*cosp - phidot*sinp + vx
        ydot = raddot*sinp + phidot*cosp + vy
        zdot = zeddot                    + vz
!
        vxdot = vraddot*cosp - vphidot*sinp + aa_nbody_cart(1) ! ivqx_cart
        vydot = vraddot*sinp + vphidot*cosp + aa_nbody_cart(2) ! ivqy_cart
        vzdot =                               aa_nbody_cart(3) ! ivqz_cart
!
      else if (lspherical_coords) then
!
!  The input position, velocities, and accelerations in spherical coordinates.
!
        rad  = xx_polar(1) ; tht  = xx_polar(2) ; phi  = xx_polar(3)
        raddot = xxdot_polar(1) ; thtdot = xxdot_polar(2) ; phidot = xxdot_polar(3)
!
        vrad = vv_polar(1) ; vtht = vv_polar(2) ; vphi = vv_polar(3)
        vraddot = vvdot_polar(1) ; vthtdot = vvdot_polar(2) ; vphidot = vvdot_polar(3)
!
!  Shortcuts.
!
        cosp = cos(phi) ; sinp = sin(phi)
        cost = cos(tht) ; sint = sin(tht)
!
!  Transform the positions and velocities to Cartesian coordinates.
!
        xp = rad*sint*cosp ; yp = rad*sint*sinp ; zp = rad*cost
!
        vx = vrad*sint*cosp + vtht*cost*cosp - vphi*sinp
        vy = vrad*sint*sinp + vtht*cost*sinp + vphi*cosp
        vz = vrad*cost      - vtht*sint
!
!  Transform the position and velocity derivatives to Cartesian coordinates.        
!  Add the cartesian velocity to position and cartesian acceleration to velocity.
!        
        xdot  =  raddot*sint*cosp +  thtdot*cost*cosp -  phidot*sinp + vx
        ydot  =  raddot*sint*sinp +  thtdot*cost*sinp +  phidot*cosp + vy
        zdot  =  raddot*cost      -  thtdot*sint                     + vz
!
        vxdot = vraddot*sint*cosp + vthtdot*cost*cosp - vphidot*sinp + aa_nbody_cart(1) !ivxq_cart)
        vydot = vraddot*sint*sinp + vthtdot*cost*sinp + vphidot*cosp + aa_nbody_cart(2) !ivyq_cart)
        vzdot = vraddot*cost      - vthtdot*sint                     + aa_nbody_cart(3) !ivzq_cart)
!
      endif
!
!  Now the time-stepping in Cartesian coordinates.
!
       xx_cart(1) = xp    ;   xx_cart(2) = yp   ;   xx_cart(3) = zp
      dxx_cart(1) = xdot  ;  dxx_cart(2) = ydot ;  dxx_cart(3) = zdot

       vv_cart(1) = vx    ;  vv_cart(2) = vy    ;  vv_cart(3) = vz
      dvv_cart(1) = vxdot ; dvv_cart(2) = vydot ; dvv_cart(3) = vzdot

      call update_position(xx_cart,dxx_cart,xx_polar,xxdot_polar)
      call update_velocity(vv_cart,dvv_cart,vv_polar,xx_polar)
!
    endsubroutine advance_particles_in_cartesian
!***********************************************************************
    subroutine update_position(xx_cart,dxx_cart,xx_polar,xxdot_polar)
!
!  Update position if N-body is used in polar coordinates.
!
!  14-feb-14:wlad/coded
!
      real, dimension(3), intent(inout) :: xx_polar,xx_cart,dxx_cart
      real, dimension(3), intent(out) :: xxdot_polar
      real :: xp,yp,zp,xdot,ydot,zdot,cosp,sinp,sint,cost,phi,tht
!
!  Update.
!
      xx_cart = xx_cart + dt_beta_ts(itsub)*dxx_cart
!
!  Convert back to polar coordinates.
!
      xp=xx_cart(1); yp=xx_cart(2); zp=xx_cart(3)
      if (lcartesian_coords) then
        xx_polar(1) = xp
        xx_polar(2) = yp
        xx_polar(3) = zp
      else if (lcylindrical_coords) then
        xx_polar(1) = sqrt(xp**2+yp**2+zp**2)
        xx_polar(2) = atan2(yp,xp)
        xx_polar(3) = zp
      else if (lspherical_coords) then
        xx_polar(1) = sqrt(xp**2+yp**2+zp**2)
        xx_polar(2) = acos(zp/xx_polar(1))
        xx_polar(3) = atan2(yp,xp)
      endif
!
      xdot=dxx_cart(1); ydot=dxx_cart(2); zdot=dxx_cart(3)
      if (lcartesian_coords) then
        xxdot_polar(1) = xdot
        xxdot_polar(2) = ydot
        xxdot_polar(3) = zdot
      elseif (lcylindrical_coords) then
        phi=xx_polar(2)
        cosp=cos(phi) ; sinp=sin(phi)
        xxdot_polar(1) =  xdot*cosp + ydot*sinp !=vrad
        xxdot_polar(2) = -xdot*sinp + ydot*cosp !=vphi
        xxdot_polar(3) =  zdot
      else if (lspherical_coords) then
        tht=xx_polar(2)
        phi=xx_polar(3)
        cost=cos(tht) ; sint=sin(tht)
        cosp=cos(phi) ; sinp=sin(phi)
        xxdot_polar(1) =  xdot*sint*cosp + ydot*sint*sinp + zdot*cost !=vrad
        xxdot_polar(2) =  xdot*cost*cosp + ydot*cost*sinp - zdot*sint !=vphi
        xxdot_polar(3) = -xdot*sinp      + ydot*cosp                  !=vz
      endif

    endsubroutine update_position
!***********************************************************************
    subroutine update_velocity(vv_cart,dvv_cart,vv_polar,xx_polar)
!
!  Update velocity if N-body is used in polar coordinates.
!
!  14-feb-14:wlad/coded
!
      !real, intent(in) :: ax,ay,az
      real :: vx,vy,vz
      real, dimension(3), intent(inout) :: vv_polar,vv_cart
      real, dimension(3), intent(in) :: dvv_cart,xx_polar
!
      real :: phi,tht
      real :: cosp,sinp,sint,cost
!
!  dvv_cart contain the accelerations transformed to Cartesian, plus
!  the n-body acceleration, also in Cartesian. Do the RK timestepping.    
!
      vv_cart = vv_cart + dt_beta_ts(itsub)*dvv_cart
!
!  Convert back to polar coordinates.
!
      vx=vv_cart(1); vy=vv_cart(2); vz=vv_cart(3)
      if (lcartesian_coords) then
        vv_polar(1) = vx
        vv_polar(2) = vy
        vv_polar(3) = vz
      else if (lcylindrical_coords) then
        phi=xx_polar(2)
        cosp=cos(phi) ; sinp=sin(phi)
        vv_polar(1) =  vx*cosp + vy*sinp !=vrad
        vv_polar(2) = -vx*sinp + vy*cosp !=vphi
        vv_polar(3) =  vz
      else if (lspherical_coords) then
        tht=xx_polar(2)
        phi=xx_polar(3)
        cost=cos(tht) ; sint=sin(tht)
        cosp=cos(phi) ; sinp=sin(phi)
        vv_polar(1) =  vx*sint*cosp + vy*sint*sinp + vz*cost !=vrad
        vv_polar(2) =  vx*cost*cosp + vy*cost*sinp - vz*sint !=vphi
        vv_polar(3) = -vx*sinp      + vy*cosp                !=vz
      endif
!
    endsubroutine update_velocity
!***********************************************************************
  endmodule PointMasses
