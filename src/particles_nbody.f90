! $Id: particles_nbody.f90,v 1.48 2007-09-07 01:10:17 wlyra Exp $
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

  use Cdata
  use Messages
  use Particles_cdata
  use Particles_sub

  implicit none

  include 'particles_nbody.h'

  real, dimension(nspar,mpvar) :: fsp
  real, dimension(nspar) :: xsp0=0.0, ysp0=0.0, zsp0=0.0
  real, dimension(nspar) :: vspx0=0.0, vspy0=0.0, vspz0=0.0
  real, dimension(nspar) :: pmass=1.,r_smooth,pmass1
  real :: delta_vsp0=1.0,totmass,totmass1
  character (len=labellen) :: initxxsp='origin', initvvsp='nothing'
  logical :: lcalc_orbit=.true.,lmigrate=.false.,lnorm=.true.
  logical :: lreset_cm=.false.,lnogravz_star=.false.,lexclude_frozen=.true.
  logical, dimension(nspar) :: lcylindrical_gravity_nbody=.false.
  logical, dimension(nspar) :: lfollow_particle=.false.
  real :: GNewton=impossible

  namelist /particles_nbody_init_pars/ &
       initxxsp, initvvsp, xsp0, ysp0, zsp0, vspx0, vspy0, vspz0, delta_vsp0, &
       pmass, r_smooth, lcylindrical_gravity_nbody, &
       lexclude_frozen, GNewton


  namelist /particles_nbody_run_pars/ &
       dsnap_par_minor, linterp_reality_check, lcalc_orbit, lreset_cm, &
       lnogravz_star,lfollow_particle, lmigrate, lexclude_frozen, GNewton

  integer, dimension(nspar,3) :: idiag_xxspar=0,idiag_vvspar=0
  integer, dimension(nspar)   :: idiag_torqint=0,idiag_torqext=0
  integer                     :: idiag_totenergy=0,idiag_totangmom=0

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
           "$Id: particles_nbody.f90,v 1.48 2007-09-07 01:10:17 wlyra Exp $")
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
      use Mpicomm,only:stop_it
!
      integer :: k
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
    endsubroutine initialize_particles_nbody
!***********************************************************************
    subroutine pencil_criteria_par_nbody()
!
!  All pencils that the Particles_nbody module depends on are specified here.
!
!  22-sep-06/wlad: adapted
!
      lpenc_requested(i_rho)=.true.
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
      real :: aux,aux2,rr,fac,parc
      integer :: k,ks

      intent (in) :: f
      intent (out) :: fp
!
!  Initial particle position.
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
      case ('fixed-cm')
!
! Ok, I have the masses and the positions of all sinks except the last, 
! which will have a position determined to fix the center of mass on 
! the center of the grid
!
        position(:,1) = xsp0; position(:,2) = ysp0; position(:,3) = zsp0
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
        aux = 0.;parc=0
        do ks=1,nspar-1 
          sma(ks)=abs(position(ks,1))
          aux=aux+pmass(ks)
          parc = parc - sma(ks)*pmass(ks)
        enddo
        pmass(nspar)=1.- aux;pmass1=1./pmass;totmass=1.;totmass1=1.
        parc = parc*totmass1
        if (aux .ge. 1.) &
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
        if (lroot) print*,'pmass =',pmass
        if (lroot) print*,'position (x)=',position(:,1)
        if (lroot) print*,'position (y)=',position(:,2)
        if (lroot) print*,'position (z)=',position(:,3)
!
! Loop through ipar to allocate the sink particles
!
        do k=1,npar_loc
          if (ipar(k) <= nspar) then
!
            print*,&
                 'initparticles_nbody. Slot for sink particle ',ipar(k),&
                 ' was at fp position ',k,&
                 ' at processor ',iproc
!
            fp(k,ixp:izp)=position(ipar(k),1:3)
!
! Correct for non-existing dimensions (not really needed, I think)
!
            if (nxgrid==1) fp(k,ixp)=x(nghost+1)
            if (nygrid==1) fp(k,iyp)=y(nghost+1)
            if (nzgrid==1) fp(k,izp)=z(nghost+1)
!
            print*,&
                 'initparticles_nbody. Sink particle ',ipar(k),&
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
    endsubroutine init_particles_nbody
!***********************************************************************
    subroutine dvvp_dt_nbody_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Add gravity from the particles to the gas
!  and the backreaction from the gas onto the particles
!
!  07-sep-06/wlad: coded
!
      use Messages, only: fatal_error
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      real, dimension (nx,3) :: ggp
      real, dimension (nx,nspar) :: rp_mn,rpcyl_mn
      real, dimension (nx) :: grav_particle,rrp,pot_energy
      real, dimension (3) :: xxpar,accg
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
      integer :: ks,k
!
      intent (in) :: f, p, fp, ineargrid
      intent (inout) :: df,dfp
!
! Output the positions of all particles
!
      if (lhydro) then
        if (headtt) then
          print*,'dvvp_dt_nbody_pencil: Add particles gravity to the gas'
          if (nxgrid/=1) print*,'dvvp_dt_nbody_pencil: Particles '//&
               'located at fsp(x)=',fsp(:,ixp)
          if (nygrid/=1) print*,'dvvp_dt_nbody_pencil: Particles '//&
               'located at fsp(y)=',fsp(:,iyp)
          if (nzgrid/=1) print*,'dvvp_dt_nbody_pencil: Particles '//&
               'located at fsp(z)=',fsp(:,izp)
        endif
!
! Initialize tmp for potential energy
!
        if (idiag_totenergy/=0) pot_energy=0.
!
! Calculate gas-particles distances
!
        do ks=1,nspar
!
! Spherical and cylindrical distances
!
          call get_radial_distance(rp_mn(:,ks),rpcyl_mn(:,ks),&
               E1_=fsp(ks,ixp),E2_=fsp(ks,iyp),E3_=fsp(ks,izp))
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
! Add this acceleration to the gas
!
          df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + ggp
!
! Backreaction of the gas onto the particles
!
          if (lmigrate) then
!
! Get the acceleration particle ks suffers due to gas gravity
!
            xxpar = fsp(ks,ixp:izp)
            call gas_gravity(p,rrp,xxpar,accg,r_smooth(ks))
!
! Add it to its dfp
!
            do k=1,npar_loc
              if (ipar(k)==ks) &
                   dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + accg(1:3)
            enddo
          endif
!
! Calculate torques for output, if needed
!
          if (ldiagnos) then
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
                   call sum_lim_mn_name(p%rho*0.5*p%u2 + pot_energy,&
                   idiag_totenergy,p)
            endif
          endif
        enddo
!
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
!  Evolution of sink particles velocities due to
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
      real, dimension (nspar,3) :: acc
!
      real, dimension (npar_loc,nspar,3) :: xxspar,vvspar
!
      real :: Omega2,invr3_ij
      integer :: i, k, ks, ki, kj, j, jvel, jpos
      logical :: lheader, lfirstcall=.true.
!
      real :: e1,e2,e3,e10,e20,e30,ev1,ev2,ev3
      real :: rad,raddot,sintht,costht,thtdot,phidot
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
!
!  Evolve the position of the sink particles due to their mutual gravity
!
      do ks=1,nspar
!
!  Calculate in the root only, then broadcast
!
        if (lroot) then
          acc(ks,:) = 0.
          do kj=1,nspar
            if (kj/=ks) then
!
              e1=fsp(ks,ixp);e10=fsp(kj,ixp)
              e2=fsp(ks,iyp);e20=fsp(kj,iyp)
              e3=fsp(ks,izp);e30=fsp(kj,izp)
!
!  Get the distances in each ortogonal component. 
!  These are NOT (x,y,z) for all.
!  For cartesian it is (x,y,z), for cylindrical (s,phi,z)
!  for spherical (r,theta,phi)
! 
              if (lcartesian_coords) then
                ev1 = e1 - e10  
                ev2 = e2 - e20  
                ev3 = e3 - e30  
              elseif (lcylindrical_coords) then
                ev1 = e1 - e10*cos(e2-e20)  
                ev2 = e10*sin(e2-e20)       
                ev3 = e3 - e30              
              elseif (lspherical_coords) then
                call stop_it("dvvp_dt_nbody: not yet implemented for "//&
                     " spherical polars")
              endif
!
!  Particles relative distance from each other
!
!  r_ij = sqrt(ev1**2 + ev2**2 + ev3**2)
!  invr3_ij = r_ij**(-3)
!
              invr3_ij = (ev1**2 + ev2**2 + ev3**2)**(-1.5)
!
!  Gravitational acceleration: g=g0/|r-r0|^3 (r-r0)
!
              acc(ks,1) = acc(ks,1) - GNewton*pmass(kj)*invr3_ij*ev1
              acc(ks,2) = acc(ks,2) - GNewton*pmass(kj)*invr3_ij*ev2
              acc(ks,3) = acc(ks,3) - GNewton*pmass(kj)*invr3_ij*ev3
!
            endif
          enddo
        endif
!
!  Broadcast particle acceleration
!
        call mpibcast_real(acc(ks,:),3)
!
!  Put it back on the dfp array on this processor
!
        do k=1,npar_loc
          if (ipar(k)==ks) then
!
!  Non-coordinate basis (all have dimension of length). 
!  The main dxx_dt of particle_dust takes care of 
!  transforming the linear velocities to angular changes 
!  in position.
!
            dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + acc(ks,1:3)
!
!  Curvilinear coordinates corrections
!
            if (lcylindrical_coords) then
!
              rad    =  fsp(ks,ixp)
              raddot =  fsp(ks,ivpx)
              phidot =  fsp(ks,ivpy)/rad
!
              dfp(k,ivpx) = dfp(k,ivpx) + rad*phidot**2
              dfp(k,ivpy) = dfp(k,ivpy) - 2*raddot*phidot
!
            elseif (lspherical_coords) then
!
              rad    = fsp(ks,ixp)
              sintht = sin(fsp(ks,iyp))
              costht = cos(fsp(ks,iyp))
              raddot = fsp(ks,ivpx)
              thtdot = fsp(ks,ivpy)/rad
              phidot = fsp(ks,ivpz)/(rad*sintht)
!             
              dfp(k,ivpx) = dfp(k,ivpx) &
                   + rad*(thtdot**2 + (sintht*phidot)**2)
              dfp(k,ivpy) = dfp(k,ivpy) &
                   - 2*raddot*thtdot + rad*sintht*costht*phidot**2
              dfp(k,ivpz) = dfp(k,ivpz) &
                   - 2*phidot*(sintht*raddot + rad*costht*thtdot)
!
            endif
          endif
        enddo
!
!  Position and velocity diagnostics (per sink particle)
!
        if (ldiagnos) then
          if (lfollow_particle(ks)) then
            do j=1,3
              jpos=j+ixp-1 ; jvel=j+ivpx-1
              xxspar(1:npar_loc,ks,j) = fsp(ks,jpos)
              vvspar(1:npar_loc,ks,j) = fsp(ks,jvel)
!
              if (idiag_xxspar(ks,j)/=0) &
                   call sum_par_name(xxspar(1:npar_loc,ks,j),idiag_xxspar(ks,j))
              if (idiag_vvspar(ks,j)/=0) &
                   call sum_par_name(vvspar(1:npar_loc,ks,j),idiag_vvspar(ks,j))
            enddo
          endif
        endif
!
      enddo
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
    subroutine gas_gravity(p,rrp,xxpar,accg,rp0)
!
! Calculates acceleration on the point (x,y,z)=xxpar
! due to the gravity of the gas.
!
! 15-sep-06/wlad : coded
!
      use Mpicomm
!
      real, dimension(nx,3) :: xxgas
      type (pencil_case) :: p
      real, dimension(nx) :: rrp,grav_gas
      real, dimension(3) :: accg,sum_loc,xxpar
      real :: dv,rp0
      integer :: k,j
!
      intent(out) :: accg
!
      if (headtt) &
           print*,'Adding gas gravity to particles'
!
      if (.not.lcartesian_coords) &
           call stop_it("particles_nbody: gas gravity not yet "//&
           "implemented for curvilinear coordinates.")
!
      if (nzgrid==1) then
        dv=dx*dy
      else
        dv=dx*dy*dz
      endif
!
! The gravity of every single gas cell - should exclude inner and outer radii...
!
! grav_gas = G*(rho*dv)*mass*r*(r**2 + r0**2)**(-1.5)
! gx = grav_gas * r\hat dot x\hat
! -> gx = grav_gas * (x-x0)/r = G*(rho*dv)*mass*(r**2+r0**2)**(-1.5) * (x-x0)
!
      grav_gas = GNewton*p%rho*dv*(rrp**2 + rp0**2)**(-1.5)
!
! Everything inside the accretion radius of the particle shouldn't exert gravity
!
      where (rrp.le.rp0)
        grav_gas = 0
      endwhere
!
      if (lexclude_frozen) then
        if (lcylinder_in_a_box) then
          where ((p%rcyl_mn.le.r_int).or.(p%rcyl_mn.ge.r_ext))
            grav_gas = 0
          endwhere
        else
          where ((p%r_mn.le.r_int).or.(p%r_mn.ge.r_ext))
            grav_gas = 0
          endwhere
        endif
      endif
!
! Sum the accelerations on this processor
! And over processors with mpireduce
!
      xxgas(:,1) = x(l1:l2) ; xxgas(:,2) = y(m) ; xxgas(:,3) = z(n)
!
      do j=1,3
        sum_loc(j) = sum(grav_gas * (xxgas(:,j) - xxpar(j)))
        call mpireduce_sum_scl(sum_loc(j),accg(j))
      enddo
!
!  Broadcast particle acceleration
!
      call mpibcast_real(accg,3)
!
    endsubroutine gas_gravity
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
              call mpisend_real(fsp(ks,:),mpvar,root,ks)
              if (ip<=6) print*,'logical for particle ',ks,&
                   ' set to true on processor ',iproc, &
                   ' with tag=',ks
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
