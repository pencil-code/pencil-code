! $Id$
!
!  This module takes care of everything related to inertial particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 6
! MAUX CONTRIBUTION 2
! CPARAM logical, parameter :: lparticles=.true.
! PENCILS PROVIDED np; rhop; vol; peh
! PENCILS PROVIDED np_rad(5); npvz(5); npuz(5); sherwood
! PENCILS PROVIDED epsp; grhop(3)
! PENCILS PROVIDED tausupersat
!
!
!***************************************************************
module Particles
!
  use Cdata
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages
  use Mpicomm
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
  use Particles_sub
  use Particles_radius
!
  implicit none
!
  include 'particles.h'
  include 'particles_common.h'
!

!
  character (len=labellen) :: interp_pol_uu ='ngp'
  character (len=labellen) :: interp_pol_oo ='ngp'
  character (len=labellen) :: interp_pol_TT ='ngp'
  character (len=labellen) :: interp_pol_rho='ngp'
  character (len=labellen) :: interp_pol_pp ='ngp'
  character (len=labellen) :: interp_pol_species='ngp'
  character (len=labellen) :: interp_pol_gradTT='ngp'
  character (len=labellen) :: interp_pol_nu='ngp'
!
  character (len=labellen), dimension (ninit) :: initxxp='nothing'
  character (len=labellen), dimension (ninit) :: initvvp='nothing'
!
  integer :: init_repeat=0       !repeat particle initialization for distance statistics
!
!



  character (len=labellen) :: ppotential='nothing'
  real :: ppower=19.,skin_factor=2.,fampl=1.
  real :: rescale_diameter=1.
  logical :: lpotential=.true.

! 
!
! ----------- About Potential ----------------------
! This module calculates all quantities related to particle-particle
! interaction. There may or may not actually be a potential, but 
! we shall use this module to calculate diagnostic of particle
! pairs. 
!--------But if we do have a potential then ..----------
! Note: psigma below is psigma = psigma_by_dx * dx ; because we always like to think of the
! range of the potential in units of dx. While calculating neighbours we look around for 
! sigma_in_grid number of grid points. This should be >= 1 .  Also, if
! psigma is larger than dx then sigma_in_grid must be larger too. 
! Note : psigma should be much smaller than the dissipation range, maybe even less than a grid-spacing
! as we are not including many terms in the Maxey-Riley equations. As far as the interaction
! with the fluid is concerned our particles are point particles with inertia. But the
! interaction potential gives them an effective radius. The interaction potential is
! typically of the form
!   V(r) = function of (r/psigma) , \xi = r/psigma
! The default potential is repulsive
!  V(r) = fampl*(1/xi)^(beta)
! with beta = 2*ppowerby2
! This potential is quite steep (almost hard-sphere) hence the effective force on a particle
! due to other particles which are within a distance of skin_factor*psigma. Particles
! within this distance are included in the neighbourlist.
!
  real :: cell_length=0.
  integer :: ncell=0
  integer :: mcellx=0,mcelly=0,mcellz=0
  integer :: arb_factor=10
  logical :: lhead_allocated=.false.
  integer, allocatable, dimension(:,:,:) :: head
  integer, allocatable,dimension(:) :: link_list
  logical :: ldragforce_gas_par, ldragforce_dust_par
  logical :: lup_as_aux=.false.
  real :: xp0, yp0, zp0, vpx0, vpy0, vpz0, delta_vp0, &
      xp1,yp1, zp1, vpx1, vpy1, vpz1, xp2, yp2, zp2, vpx2, vpy2, vpz2, &
      xp3, yp3, zp3, vpx3, vpy3, vpz3
  real :: cdtp=0.2, cdtpgrav=0.1, cdtp_drag=0.2,gas_nu
  real :: Lx0=0.0, Ly0=0.0, Lz0=0.0
  character (len=labellen) :: draglaw='none'
  integer :: idiag_xpm=0, idiag_ypm=0, idiag_zpm=0      ! DIAG_DOC: $x_{part}$
  integer :: idiag_xp2m=0, idiag_yp2m=0, idiag_zp2m=0   ! DIAG_DOC: $x^2_{part}$
  integer :: idiag_vrelpabsm=0                          ! DIAG_DOC: $\rm{Absolute value of mean relative velocity}$
  integer :: idiag_rpm=0, idiag_rp2m=0
  integer :: idiag_vpxm=0, idiag_vpym=0, idiag_vpzm=0   ! DIAG_DOC: $u_{part}$
  integer :: idiag_vpx2m=0, idiag_vpy2m=0, idiag_vpz2m=0 ! DIAG_DOC: $u^2_{part}$
  integer :: idiag_ekinp=0     ! DIAG_DOC: $E_{kin,part}$
  integer :: idiag_vdotupm=0
  integer :: idiag_vpxmax=0, idiag_vpymax=0, idiag_vpzmax=0, idiag_vpmax=0 ! DIAG_DOC: $MAX(u_{part})$
  integer :: idiag_vpxmin=0, idiag_vpymin=0, idiag_vpzmin=0    ! DIAG_DOC: $MIN(u_{part})$
  integer :: idiag_vpxvpym=0, idiag_vpxvpzm=0, idiag_vpyvpzm=0
  integer :: idiag_rhopvpxm=0, idiag_rhopvpym=0, idiag_rhopvpzm=0
  integer :: idiag_rhopvpxt=0, idiag_rhopvpyt=0, idiag_rhopvpzt=0
  integer :: idiag_rhopvpysm=0
  integer :: idiag_lpxm=0, idiag_lpym=0, idiag_lpzm=0
  integer :: idiag_lpx2m=0, idiag_lpy2m=0, idiag_lpz2m=0
  integer :: idiag_npm=0, idiag_np2m=0, idiag_npmax=0, idiag_npmin=0 ! DIAG_DOC: $\rm{mean particle number density}$
  integer :: idiag_dtdragp=0
  integer :: idiag_nparmin=0, idiag_nparmax=0, idiag_npargone=0
  integer :: idiag_nparsum=0
  integer :: idiag_rhopm=0, idiag_rhoprms=0, idiag_rhop2m=0, idiag_rhopmax=0
  integer :: idiag_rhopmin=0, idiag_decollp=0, idiag_rhopmphi=0
  integer :: idiag_epspmin=0, idiag_epspmax=0, idiag_epspm=0
  integer :: idiag_npmx=0, idiag_npmy=0, idiag_npmz=0
  integer :: idiag_rhopmx=0, idiag_rhopmy=0, idiag_rhopmz=0
  integer :: idiag_rhop2mx=0, idiag_rhop2my=0, idiag_rhop2mz=0
  integer :: idiag_epspmx=0, idiag_epspmy=0, idiag_epspmz=0
  integer :: idiag_mpt=0, idiag_dedragp=0, idiag_rhopmxy=0, idiag_rhopmr=0
  integer :: idiag_sigmap=0
  integer :: idiag_dvpx2m=0, idiag_dvpy2m=0, idiag_dvpz2m=0
  integer :: idiag_dvpm=0, idiag_dvpmax=0, idiag_epotpm=0
  integer :: idiag_rhopmxz=0, idiag_nparpmax=0, idiag_npmxy=0
  integer :: idiag_eccpxm=0, idiag_eccpym=0, idiag_eccpzm=0
  integer :: idiag_eccpx2m=0, idiag_eccpy2m=0, idiag_eccpz2m=0
  integer :: idiag_vprms=0, idiag_vpyfull2m=0, idiag_deshearbcsm=0
  integer :: idiag_Shm=0
  integer, dimension(ninit)  :: idiag_npvzmz=0, idiag_nptz=0
  integer, dimension(ninit)  :: idiag_npuzmz=0
!------------------------!
  namelist /particles_init_pars/ &
      initxxp, initvvp, xp0, yp0, zp0, vpx0, vpy0, vpz0, delta_vp0, &
      ldragforce_gas_par, ldragforce_dust_par, bcpx, bcpy, bcpz, &
      rhopmat, xp1, &
      yp1, zp1, vpx1, vpy1, vpz1, xp2, yp2, zp2, vpx2, vpy2, vpz2, &
      xp3, yp3, zp3, vpx3, vpy3, vpz3, &
      lpotential,arb_factor,ppotential,skin_factor,rescale_diameter 
!
  namelist /particles_run_pars/ &
      bcpx, bcpy, bcpz,  &
      ldragforce_gas_par
!-----------------------
  contains
!***********************************************************************
    subroutine register_particles
!
!  Set up indices for access to the fp and dfp arrays
!
!  29-dec-04/anders: coded
!
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  Indices for particle position.
!
      ixp=npvar+1
      pvarname(npvar+1)='ixp'
      iyp=npvar+2
      pvarname(npvar+2)='iyp'
      izp=npvar+3
      pvarname(npvar+3)='izp'
!
!  Indices for particle velocity.
!
      ivpx=npvar+4
      pvarname(npvar+4)='ivpx'
      ivpy=npvar+5
      pvarname(npvar+5)='ivpy'
      ivpz=npvar+6
      pvarname(npvar+6)='ivpz'
!
!  Increase npvar accordingly.
!
      npvar=npvar+6
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_particles','npvar > mpvar')
      endif
!
    endsubroutine register_particles
!***********************************************************************
    subroutine initialize_particles(f,fp)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  29-dec-04/anders: coded
!   5-mar-15/MR: reference state included in calculation of mean density
!
      use EquationOfState, only: rho0, cs0
      use SharedVariables, only: put_shared_variable, get_shared_variable
      use Density, only: mean_density
      use Viscosity, only: getnu
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mpar_loc,mparray), intent (in) :: fp
!
      real :: nu_
      character (len=labellen) :: ivis='nothing'
      integer :: ierr, jspec
      logical, pointer :: lshearadvection_as_shift
!
!  This module is incompatible with particle block domain decomposition.
!
      if (lparticles_blocks) then
        if (lroot) then
          print*, 'initialize_particles: must use PARTICLES =  PARTICLES_DUST_BLOCKS'
          print*, '                      with particle block domain decomposition'
        endif
        call fatal_error('initialize_particles','')
      endif
      !
! prepare array for fluid velocity at particle positions if necessary.
!  Register an extra aux slot for uup if requested. This is needed
!  for calculating the correlation time from <u.intudt>. For this to work
!  you must reserve enough auxiliary workspace by setting, for example,
!     ! MPAUX CONTRIBUTION 3
!  in the beginning of your src/cparam.local file, *before* setting
!  ncpus, nprocy, etc.
      if (lup_as_aux) then
         iuf = mpvar+npaux+1
         iufx=iuf;iufy=iuf+1;iufz=iuf+2
         npaux=npaux+3
      endif
! Check that the fp and dfp arrays are big enough.
      if (npaux > mpaux) then
         if (lroot) write (0,*) 'npaux = ', npaux, ', mpaux = ', mpaux
         call fatal_error('register_particles_potential: npaux > mpaux','')
      endif
!
!  Inverse material density.
!
      if (rhopmat/=0.0) rhopmat1=1/rhopmat
!
!  Gas density (also sometime gas-viscosity) is needed for back-reaction friction force.
!
      if (ldragforce_gas_par) then
         if (.not. ldensity) then
            if (lroot) then
               print*, 'initialize_particles: friction force on gas only works '
               print*, '                      together with gas density module!'
            endif
            call fatal_error('initialize_particles','density not found!')
         endif
      
!
!  Find the kinematic viscosity
!
         call getnu(nu_input=nu_,ivis=ivis)
         select case(ivis)
         case('nu-const')
            gas_nu=nu_
         case default
            call fatal_error("particles_potential","works only with nu_const")
         endselect
      endif
!
!
! The size of a cell is twice the radius of the biggest particle. This assumes that
! the size of the particles are now going to change over time. Otherwise the input      
! parameter cell_length, if not equal to zero, sets the size of the cell. 
!
      if (cell_length.eq.0.) then 
         cell_length=2*maxval(fp(:,iap))
      endif
!
! the following line assumes that the domain is roughly size in all three      
! directions. If not, we need to code some more
!      
      ncell=int(abs(x(l2)-x(l1))/cell_length)+1
      cell_length=(x(l2)-x(l1))/ncell
!
! Assuming uniform distribution we can estimate the number of particles
! in a slab. These number are then multiplied
! an arbitrary factor (arb_factor) for which the default value is 10        
!
      nslab=arb_factor*(npar/ncpus)/ncell
!
! If we are using many processors then our domain effectively includes
! three neighbouring processors in each directions.
!
      mcellx=ncell;mcelly=ncell;mcellz=ncell
!
! Allocate the arrays head and link_list (only if they have
! not been allocated before)      
!
      if(.not.lhead_allocated) then
         if (lmpicomm) then
            lpar_max=max(arb_factor*(npar/ncpus)+6*nslab,mpar_loc)
            allocate(fp_buffer_in(nslab,mparray))
            allocate(fp_buffer_out(nslab,mparray))
         else
           lpar_max=mpar_loc         
        endif
        allocate(head(-1:mcellx,-1:mcelly,-1:mcellz))
        allocate(link_list(lpar_max))
!
! We also need to allocate a larger array in case of parallel communications
!
        allocate(fpwn(lpar_max,mparray))
        lhead_allocated=.true.
        fpwn=0.
        head=0
        link_list=0
      endif


      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles
!***********************************************************************
    subroutine init_particles(f,fp,ineargrid)
!
!  Initial positions and velocities of dust particles.
!
!  29-dec-04/anders: coded
!
      use EquationOfState, only: beta_glnrho_global, cs20
      use General, only: random_number_wrapper, normal_deviate
      use Mpicomm, only: mpireduce_sum, mpibcast_real
      use InitialCondition, only: initial_condition_xxp, initial_condition_vvp
      use Particles_diagnos_dv, only: repeated_init
!
      real, dimension (mx,my,mz,mfarray), intent (out) :: f
      real, dimension (mpar_loc,mparray), intent (out) :: fp
      integer, dimension (mpar_loc,3), intent (out) :: ineargrid
      real, dimension (mpar_loc) :: rr_tmp, az_tmp
!
      real, dimension (3) :: uup, Lxyz_par, xyz0_par, xyz1_par
      real :: vpx_sum, vpy_sum, vpz_sum
      real :: r, p, q, px, py, pz, eps, cs, k2_xxp, rp2
      real :: dim1, npar_loc_x, npar_loc_y, npar_loc_z, dx_par, dy_par, dz_par
      real :: rad,rad_scl,phi,tht,tmp,OO,xx0,yy0,r2
      integer :: l, j, k, ix0, iy0, iz0, n_kill
      logical :: lequidistant=.false.
      real :: rpar_int,rpar_ext
!
        Lxyz_par=Lxyz_loc
        xyz0_par=xyz0_loc
        xyz1_par=xyz1_loc
!
!  Initial particle position.
!
      do j=1,ninit
!
        select case (initxxp(j))
!
        case ('nothing')
          if (lroot .and. j==1) print*, 'init_particles: nothing'
!
        case ('origin')
          if (lroot) print*, 'init_particles: All particles at origin'
          fp(1:npar_loc,ixp:izp)=0.0
!
        case ('zero-z')
          if (lroot) print*, 'init_particles: Zero z coordinate'
          fp(1:npar_loc,izp)=0.0
!
        case ('constant')
          if (lroot) &
              print*, 'init_particles: All particles at x,y,z=', xp0, yp0, zp0
          fp(1:npar_loc,ixp)=xp0
          fp(1:npar_loc,iyp)=yp0
          fp(1:npar_loc,izp)=zp0
          !
        case ('random-constz')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k=1,npar_loc
            if (nxgrid/=1) then
              call random_number_wrapper(r)
              fp(k,ixp)=r
            endif
            if (nygrid/=1) then
              call random_number_wrapper(r)
              fp(k,iyp)=r
            endif
          enddo
          if (nxgrid/=1) &
              fp(1:npar_loc,ixp)=xyz0_par(1)+fp(1:npar_loc,ixp)*Lxyz_par(1)
          if (nygrid/=1) &
              fp(1:npar_loc,iyp)=xyz0_par(2)+fp(1:npar_loc,iyp)*Lxyz_par(2)
          if (nzgrid/=1) &
              fp(1:npar_loc,izp)=zp0
!
        case ('random')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k=1,npar_loc
            if (nxgrid/=1) then
              call random_number_wrapper(r)
              fp(k,ixp)=r
            endif
            if (nygrid/=1) then
              call random_number_wrapper(r)
              fp(k,iyp)=r
            endif
            if (nzgrid/=1) then
              call random_number_wrapper(r)
              fp(k,izp)=r
            endif
          enddo
          if (nxgrid/=1) &
              fp(1:npar_loc,ixp)=xyz0_par(1)+fp(1:npar_loc,ixp)*Lxyz_par(1)
          if (nygrid/=1) &
              fp(1:npar_loc,iyp)=xyz0_par(2)+fp(1:npar_loc,iyp)*Lxyz_par(2)
          if (nzgrid/=1) &
              fp(1:npar_loc,izp)=xyz0_par(3)+fp(1:npar_loc,izp)*Lxyz_par(3)
!
        case default
          call fatal_error('init_particles','Unknown value initxxp="'//trim(initxxp(j))//'"')
!
        endselect
!
      enddo ! do j=1,ninit
!
!  Interface for user's own initial condition for position
!
      if (linitial_condition) call initial_condition_xxp(f,fp)
!
!  Particles are not allowed to be present in non-existing dimensions.
!  This would give huge problems with interpolation later.
!
      if (nxgrid==1) fp(1:npar_loc,ixp)=x(nghost+1)
      if (nygrid==1) fp(1:npar_loc,iyp)=y(nghost+1)
      if (nzgrid==1) fp(1:npar_loc,izp)=z(nghost+1)
!
      if (init_repeat/=0) call repeated_init(fp,init_repeat)
!
!  Redistribute particles among processors (now that positions are determined).
!
      call boundconds_particles(fp,ipar)
!
!  Map particle position on the grid.
!
      call map_nearest_grid(fp,ineargrid)
      call map_xxp_grid(f,fp,ineargrid)
!
!  Initial particle velocity.
!
      do j=1,ninit
!
        select case (initvvp(j))
!
        case ('nothing')
          if (lroot.and.j==1) print*, 'init_particles: No particle velocity set'
!
        case ('zero')
          if (lroot) print*, 'init_particles: Zero particle velocity'
          fp(1:npar_loc,ivpx:ivpz)=0.0
!
        case ('constant')
          if (lroot) print*, 'init_particles: Constant particle velocity'
          if (lroot) &
              print*, 'init_particles: vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
          if (lcylindrical_coords) then
            fp(1:npar_loc,ivpx)=vpx0*cos(fp(k,iyp))+vpy0*sin(fp(k,iyp))
            fp(1:npar_loc,ivpy)=vpy0*cos(fp(k,iyp))-vpx0*sin(fp(k,iyp))
            fp(1:npar_loc,ivpz)=vpz0
          else
            fp(1:npar_loc,ivpx)=vpx0
            fp(1:npar_loc,ivpy)=vpy0
            fp(1:npar_loc,ivpz)=vpz0
          endif
!
        case ('random')
          if (lroot) print*, 'init_particles: Random particle velocities; '// &
              'delta_vp0=', delta_vp0
          do k=1,npar_loc
            call random_number_wrapper(r)
            fp(k,ivpx) = fp(k,ivpx) + delta_vp0*(2*r-1)
            call random_number_wrapper(r)
            fp(k,ivpy) = fp(k,ivpy) + delta_vp0*(2*r-1)
            call random_number_wrapper(r)
            fp(k,ivpz) = fp(k,ivpz) + delta_vp0*(2*r-1)
          enddo
!
        case ('follow-gas')
          if (lroot) &
              print*, 'init_particles: Particle velocity equal to gas velocity'
          do k=1,npar_loc
            call interpolate_linear(f,iux,iuz,fp(k,ixp:izp),uup, &
                ineargrid(k,:),0,0)
            fp(k,ivpx:ivpz) = uup
          enddo
!
        case default
          call fatal_error('init_particles','Unknown value initvvp="'//trim(initvvp(j))//'"')
!
        endselect
!
      enddo ! do j=1,ninit
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_vvp(f,fp)
!
!  Sort particles (must happen at the end of the subroutine so that random
!  positions and velocities are not displaced relative to when there is no
!  sorting).
!
      call sort_particles_imn(fp,ineargrid,ipar)
!
    endsubroutine init_particles
!***********************************************************************
    subroutine insert_lost_particles(f,fp,ineargrid)
!
!  14-oct-12/dhruba: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray), intent (inout) :: fp
      integer, dimension (mpar_loc,3), intent (inout) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine insert_lost_particles
!***********************************************************************
    subroutine pencil_criteria_particles
!
!  All pencils that the Particles module depends on are specified here.
!
!  20-04-06/anders: coded
!
      if (ldragforce_gas_par) then
        lpenc_requested(i_epsp)=.true.
        lpenc_requested(i_np)=.true.
        lpenc_requested(i_rho1) = .true.
        if (draglaw=='epstein') lpenc_requested(i_cs2)=.true.
      endif
!
      if (idiag_npm/=0 .or. idiag_np2m/=0 .or. idiag_npmax/=0 .or. &
          idiag_npmin/=0 .or. idiag_npmx/=0 .or. idiag_npmy/=0 .or. &
          idiag_npmz/=0 .or. idiag_nparpmax/=0) lpenc_diagnos(i_np)=.true.
      if (idiag_rhopm/=0 .or. idiag_rhoprms/=0 .or. idiag_rhop2m/=0 .or. &
          idiag_rhopmax/=0 .or. idiag_rhopmin/=0 .or. idiag_rhopmphi/=0 .or. &
          idiag_rhopmx/=0 .or. idiag_rhopmy/=0 .or. idiag_rhopmz/=0) &
          lpenc_diagnos(i_rhop)=.true.
      if (idiag_rhop2mx /= 0 .or. idiag_rhop2my /= 0 .or. idiag_rhop2mz /= 0) lpenc_diagnos(i_rhop) = .true.
      if (idiag_dedragp/=0 .or. idiag_decollp/=0) then
        lpenc_diagnos(i_TT1)=.true.
        lpenc_diagnos(i_rho1)=.true.
      endif
      if (idiag_epspmx/=0 .or. idiag_epspmy/=0 .or. idiag_epspmz/=0 .or. &
          idiag_epspmin/=0 .or. idiag_epspmax/=0 .or. idiag_epspm/=0) &
          lpenc_diagnos(i_epsp)=.true.
      if (idiag_rhopmxy/=0 .or. idiag_rhopmxz/=0 .or. idiag_rhopmphi/=0) &
          lpenc_diagnos2d(i_rhop)=.true.
      if (idiag_npmxy/=0 ) lpenc_diagnos2d(i_np)=.true.
      if (idiag_sigmap /= 0) lpenc_diagnos2d(i_rhop) = .true.
!
      if (maxval(idiag_npvzmz) > 0) lpenc_requested(i_npvz)=.true.
      if (maxval(idiag_npuzmz) > 0) lpenc_requested(i_npuz)=.true.
      if (maxval(idiag_nptz) > 0)   lpenc_requested(i_np_rad)=.true.
      if (idiag_Shm /= 0) lpenc_requested(i_sherwood)=.true.
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
      if (lpencil_in(i_rhop) .and. irhop==0) then
        lpencil_in(i_np)=.true.
      endif
!
      if (lpencil_in(i_epsp)) then
        lpencil_in(i_rhop)=.true.
        lpencil_in(i_rho1)=.true.
      endif
!
      if (lsupersat) lpencil_in(i_tausupersat)=.true.
!
    endsubroutine pencil_interdep_particles
!***********************************************************************
    subroutine calc_pencils_particles(f,p)
!
      use Sub, only: grad
!
!  Calculate Particles pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  16-feb-06/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      if (lpencil(i_np)) then
        if (inp/=0) then
          p%np=f(l1:l2,m,n,inp)
        else
          p%np=0.0
        endif
      endif
!
    endsubroutine calc_pencils_particles
!***********************************************************************
    subroutine dxxp_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of dust particle position.
!
!  02-jan-05/anders: coded
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mx,my,mz,mvar), intent (inout) :: df
      real, dimension (mpar_loc,mparray), intent (in) :: fp
      real, dimension (mpar_loc,mpvar), intent (inout) :: dfp
      integer, dimension (mpar_loc,3), intent (in) :: ineargrid
!
      logical :: lheader, lfirstcall=.true.
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
!  If pointmasses are used do the evolution in that module, for better
!  conservation of the Jacobi constant.
!
      if (.not.lpointmasses) then
        if (lcartesian_coords) then
!
          if (nxgrid/=1) &
              dfp(1:npar_loc,ixp) = dfp(1:npar_loc,ixp) + fp(1:npar_loc,ivpx)
          if (nygrid/=1) &
              dfp(1:npar_loc,iyp) = dfp(1:npar_loc,iyp) + fp(1:npar_loc,ivpy)
          if (nzgrid/=1) &
              dfp(1:npar_loc,izp) = dfp(1:npar_loc,izp) + fp(1:npar_loc,ivpz)
!
        elseif (lcylindrical_coords) then
          if (nxgrid/=1) &
              dfp(1:npar_loc,ixp) = dfp(1:npar_loc,ixp) + fp(1:npar_loc,ivpx)
          if (nygrid/=1) &
              dfp(1:npar_loc,iyp) = dfp(1:npar_loc,iyp) + &
                fp(1:npar_loc,ivpy)/max(fp(1:npar_loc,ixp),tini)
          if (nzgrid/=1) &
              dfp(1:npar_loc,izp) = dfp(1:npar_loc,izp) + fp(1:npar_loc,ivpz)
!
        elseif (lspherical_coords) then
!
          if (nxgrid/=1) &
              dfp(1:npar_loc,ixp) = dfp(1:npar_loc,ixp) + fp(1:npar_loc,ivpx)
          if (nygrid/=1) &
              dfp(1:npar_loc,iyp) = dfp(1:npar_loc,iyp) + &
                fp(1:npar_loc,ivpy)/max(fp(1:npar_loc,ixp),tini)
          if (nzgrid/=1) &
              dfp(1:npar_loc,izp) = dfp(1:npar_loc,izp) + &
                fp(1:npar_loc,ivpz)/(max(fp(1:npar_loc,ixp),tini)*sin(fp(1:npar_loc,iyp)))
        endif
      endif
!
!  With shear there is an extra term due to the background shear flow.
!
      if (lshear.and.nygrid/=1) dfp(1:npar_loc,iyp) = &
          dfp(1:npar_loc,iyp) - qshear*Omega*fp(1:npar_loc,ixp)
!
      if (lfirstcall) lfirstcall=.false.
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dxxp_dt
!***********************************************************************
    subroutine dvvp_dt(f,df,p,fp,dfp,ineargrid)
!
!  Evolution of dust particle velocity.
!
!  29-dec-04/anders: coded
!
      use Diagnostics
      use EquationOfState, only: cs20
      use Sub, only: periodic_fold_back
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      real, dimension (mx,my,mz,mvar), intent (inout) :: df
      real, dimension (mpar_loc,mparray), intent (in) :: fp
      real, dimension (mpar_loc,mpvar), intent (inout) :: dfp
      integer, dimension (mpar_loc,3), intent (in) :: ineargrid
      type (pencil_case) :: p
      intent (in) :: p
      integer,parameter :: Nnab=13
      integer,dimension(Nnab,3) :: neighbours
      integer :: inab,ix2,ip,jp,kp
      real,dimension(3) :: xxij,vvij
      integer,dimension(3) :: cell_vec
      real,dimension(3) :: xxip,my_origin
!
      real :: Omega2
      integer :: npar_found,i
      logical :: lheader, lfirstcall=.true.
!
! first fill up fpwn with particles in local fp (do this even 
! in the serial case )
! 
      lpar_loc=npar_loc
      my_origin(1)=x(l1);my_origin(2)=y(m1);my_origin(3)=z(n1)
      do ip=1,npar_loc
         fpwn(ip,:) = fp(ip,:)
!
! Now make the linked list for the local particles         
!
         xxip=fpwn(ip,ixp:izp)
         cell_vec=floor((xxip-my_origin)/cell_length)
         link_list(ip)=head(cell_vec(1),cell_vec(2),cell_vec(3))
         head(cell_vec(1),cell_vec(2),cell_vec(3))=ip
         call accn_single_particle(dfp,f,fp,p,ineargrid,ip)
      enddo
      if (ldiagnos) call single_particle_diagnos(fp)
!
! We need to know about particles in neighbouring processors      
! only if we are dealing with mutual interaction potential between
! particles. Otherwise we merely access the particles by their
! link-list, this keeps the memory access in the f arrays somewhat
! local.
! Now modify the link list for particles in neighbouring
! processors. This is different in parallel and serial case.
! In the latter, this merely takes care of periodic boundary       
! conditions.
!         
      call particles_neighbour_proc()
!
!
      ! Now access all the cells and calculate their acceleration for
      ! the particles in them.      
!
      do iz=-1,mcellz;do iy=-1,mcelly;do ix=-1,mcellx
         ip = head(ix,iy,iz)
!
! within the same cell 
!
         do while (ip.ne.0) 
            jp = link_list(ip)
            do while (jp.ne.0)
               xxij= fpwn(jp,ixp:izp)-fpwn(ip,ixp:izp)
               vvij=fpwn(jp,ivpx:ivpz)-fpwn(ip,ivpx:ivpz)
               if (lpotential) &
                    call two_particle_int(dfp,ip,jp,xxij,vvij)
               jp = link_list(jp)
            enddo
! Now for neighbouring cells
               
            call get_cell_neighbours(ix,iy,iz,neighbours) 
            do inab = 1,Nnab						
               ix2 = neighbours(inab,1)      
               iy2 = neighbours(inab,2)
               iz2 = neighbours(inab,3)
               if (( ix2.eq.-1).or.(iy2.eq.-1).or.(iz2.eq.-1)&
                 .or.( ix2.eq.mcellx+1).or.(iy2.eq.mcelly+1).or.(iz2.eq.mcellz+1) ) then
               else
                  kp = head(ix2,iy2,iz2)
                  do while (kp.ne.0) 
                     xxij= fpwn(kp,ixp:izp)-fpwn(ip,ixp:izp)
                     call periodic_fold_back(xxij, Lxyz)
                     vvij=fpwn(kp,ivpx:ivpz)-fpwn(ip,ivpx:ivpz)
                     if (lpotential) &
                          call two_particle_int(dfp,ip,kp,xxij,vvij)
                     kp = link_list(kp)       
                  enddo
               endif
            enddo
            ip = link_list(ip)
         enddo
!
      enddo; enddo; enddo
!
! loop over cells done
!      

      call keep_compiler_quiet(f)
!
    endsubroutine dvvp_dt
!***********************************************************************
    subroutine accn_single_particle(dfp,f,fp,p,ineargrid,kp)

      integer,intent(in) :: kp
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mpar_loc,mparray), intent (in) :: fp
      real, dimension (mpar_loc,mpvar), intent (inout) :: dfp
      integer, dimension (mpar_loc,3), intent (in) :: ineargrid
      type (pencil_case) :: p
      intent (in) :: p
      real :: mp_k,radp,taup1
      real :: rhof
      real,dimension(3) :: xxp,vvp,uup
      integer :: inx0
!
      if (ldragforce_gas_par) then
         radp=fp(kp,iap)
         xxp=fp(kp,ixp:izp)
         vvp=fp(kp,ivpx:ivpz)
         call interpolate_linear(f,ilnrho,xxp,rhof,ineargrid(kp,:),0,0)
         if (.not.ldensity_nolog) rhof=exp(rhof)
         call interpolate_linear(f,iux,iuz,xxp,uup,ineargrid(kp,:),0,0)
         if (luf_as_aux) fp(kp,iufx:iufz)=uup
!   get the inverse stopping time
         select case (draglaw)
!
         case ('stokes')
            taup1=18.0*gas_nu/((rhopmat/rhof)*4.*radp**2)
!
         case('epstein')
            inx0=ineargrid(kp,1)-nghost
            taup1=(sqrt(p%cs2(inx0))*rhof)/(radp*rhopmat)
!
         case default
            call fatal_error("particles_potential:","no draglaw of this name")
         endselect
         dfp(kp,ivpx:ivpz) = dfp(kp,ivpx:ivpz)+taup1*(uup-vvp)
      endif
    endsubroutine accn_single_particle
!***********************************************************************
    subroutine construct_link_list(plist_min,plist_max)
      integer :: ip
      integer, intent(in) :: plist_min,plist_max
      integer,dimension(3) :: cell_vec
      real,dimension(3) :: xxip,my_origin
!
      my_origin(1)=x(l1);my_origin(2)=y(m1);my_origin(3)=z(n1)
      do ip=plist_min,plist_max
        xxip=fpwn(ip,ixp:izp)
        cell_vec=floor((xxip-my_origin)/cell_length)
        link_list(ip)=head(cell_vec(1),cell_vec(2),cell_vec(3))
        head(cell_vec(1),cell_vec(2),cell_vec(3))=ip
      enddo
!
    endsubroutine construct_link_list
!***********************************************************************
    subroutine assimilate_incoming(npbuf)
      integer,intent(in) :: npbuf
      fpwn(lpar_loc+1:lpar_loc+npbuf,:)=fp_buffer_in(1:npbuf,:)
      call construct_link_list(lpar_loc+1,lpar_loc+npbuf)
      lpar_loc=lpar_loc+npbuf
!
    endsubroutine assimilate_incoming
!***********************************************************************
    subroutine get_boundary_particles(idirn,porm,npbuf)
!
! fp_buffer in known globally
!
      integer, intent(in) :: idirn,porm
      integer, intent(out) :: npbuf
!
      select case(idirn)
      case(3)
         if (porm.eq.1) then
            call make_fpbuffer(mcellz-1,mcellz-1,0,mcelly-1,0,mcellx-1,npbuf)
         else
            call make_fpbuffer(0,0,0,mcelly-1,0,mcellx-1,npbuf)
         endif
      case(2)
         if (porm.eq.1) then
            call make_fpbuffer(-1,mcellz,mcelly-1,mcelly-1,0,mcellx-1,npbuf)
         else
            call make_fpbuffer(-1,mcellz,0,0,0,mcellx-1,npbuf)
         endif
      case(1)
         if (porm.eq.1) then
            call make_fpbuffer(-1,mcellz,-1,mcelly,mcellx-1,mcellx-1,npbuf)
         else
            call make_fpbuffer(-1,mcellz,-1,mcelly,0,0,npbuf)
         endif
         case default
          !
          !  Catch unknown values
          !
          call fatal_error("particles_potential", &
              "get_boundary_particles is called with wrong idirn")
!
        endselect
!
      endsubroutine get_boundary_particles
!***********************************************************************
    subroutine make_fpbuffer(izmin,izmax,iymin,iymax,ixmin,ixmax,npbuf)
!
! fp_buffer_out is known globally
!
      integer, intent(in) :: izmin,izmax,iymin,iymax,ixmin,ixmax
      integer, intent(out) :: npbuf
      integer :: ipbuf
      fp_buffer_out=0.
      npbuf=0
      ipbuf=0
      do iz=izmin,izmax
         do iy=iymin,iymax
            do ix=ixmin,ixmax
               ip = head(ix,iy,iz)
               do while (ip.ne.0)
                  ipbuf=ipbuf+1
                  fp_buffer_out(ipbuf,:)=fpwn(ip,:)
                  ip = link_list(ip)
               enddo ! loop all the particles in a cell ends
            enddo
         enddo
      enddo
      npbuf=ipbuf
    endsubroutine make_fpbuffer
!***********************************************************************
    subroutine two_particle_int(dfp,ip,jp,xxij,vvij)
!
      use Diagnostics
      use particles_radius, only: get_mass_from_radius
!
      real, dimension (mpar_loc,mpvar) :: dfp
!
      integer,intent(in) :: ip,jp
      real, dimension(3),intent(in) :: xxij,vvij
      real,dimension(3) :: force_ij
      real :: mp_i,mp_j
!!---------------------------------
!
! The particle may be in a different processor. If that is the
! case then we do nothing. 
!
      
      if (ip .le. npar_loc) then
         call get_interaction_force(force_ij,xxij,ip,jp)
         call get_mass_from_radius(mp_i,fpwn,ip)
         dfp(ip,ivpx:ivpz) = dfp(ip,ivpx:ivpz)+force_ij/mp_i
         if (ldiagnos) call two_particle_diagnos(xxij,vvij)
         if (jp .le. npar_loc) then
            call get_mass_from_radius(mp_j,fpwn,jp)
            !don't forget Newton's 3rd law
            dfp(jp,ivpx:ivpz) = dfp(jp,ivpx:ivpz)-force_ij/mp_j
         endif
      endif
!
    endsubroutine two_particle_int
!***********************************************************************
    subroutine get_interaction_force(force_ij,RR,ip,jp)
      integer, intent(in) :: ip,jp
      real,dimension(3),intent(in) :: RR
      real,dimension(3),intent(out) :: force_ij
      real :: RR_mod,sigma
      real,dimension(3) :: Rcap
      real :: radiusi,radiusj,diameter_ij,force_amps
!
      select case (ppotential)
      case ('rep-power-law-cutoff')
!
! repulsive power law, force = -1/(RR/(radi+radj))^p
!
        RR_mod=sqrt(RR(1)*RR(1)+RR(2)*RR(2)+RR(3)*RR(3))
        Rcap=RR/RR_mod
        radiusi=fpwn(ip,iap)
        radiusj=fpwn(jp,iap)
        diameter_ij=rescale_diameter*(radiusi+radiusj)
        sigma=RR_mod/diameter_ij
        if (sigma .lt. 1.) then
          force_ij=0.
        else
          force_amps=fampl*sigma**(-ppower)
          force_ij=-force_amps*Rcap
        endif
      case default
        call fatal_error('particles_potential: no potential coded ','get_interaction_force')
      endselect
!
    endsubroutine get_interaction_force
!***********************************************************************
    subroutine get_cell_neighbours(ix,iy,iz,neighbours)
      integer,intent(in) :: ix,iy,iz
      integer,parameter :: Nnab=13
      integer,dimension(Nnab,3) :: neighbours
      integer :: il,im,in,inab
      !
! All the 9 neighbours above (in z) are covered.
      !
      inab=0
      !
      ! neighbours at one level up (9 of them) in z
      !
      do il=-1,1
         do im=-1,1
            inab=inab+1
            neighbours(inab,1) = ix+il
            neighbours(inab,2) = iy+im
            neighbours(inab,3) = iz+1
         enddo
      enddo
      !
      !neighbours in same z (selected 4)
      !
      do il=-1,1
         inab=inab+1
         neighbours(inab,1)=ix+il
         neighbours(inab,2)=iy+1
         neighbours(inab,3)=iz
      enddo
      !
      !and the final one
      !
      inab=inab+1
      neighbours(inab,1) = ix + 1
      neighbours(inab,2) = iy
      neighbours(inab,3) = iz
!
    endsubroutine get_cell_neighbours
!***********************************************************************
    subroutine two_particle_diagnos(xxij,vvij)
      real, dimension(3) :: xxij,vvij
    endsubroutine two_particle_diagnos
!***********************************************************************
    subroutine single_particle_diagnos(fp)
!      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray), intent (in) :: fp
!
      if (idiag_vpx2m/=0) &
           call sum_par_name(fp(1:npar_loc,ivpx)**2,idiag_vpx2m)
      if (idiag_vpy2m/=0) &
           call sum_par_name(fp(1:npar_loc,ivpy)**2,idiag_vpy2m)
      if (idiag_vpz2m/=0) &
           call sum_par_name(fp(1:npar_loc,ivpz)**2,idiag_vpz2m)
      if (idiag_vprms/=0) &
           call sum_par_name((fp(1:npar_loc,ivpx)**2 &
                              +fp(1:npar_loc,ivpy)**2 &
                              +fp(1:npar_loc,ivpz)**2),idiag_vprms,lsqrt=.true.)
      if (idiag_ekinp/=0) then
         call sum_par_name(0.5*(4./3.)*pi*(fp(1:npar_loc,iap)**3)*rhopmat*&
              ( fp(1:npar_loc,ivpx)**2 &
              +fp(1:npar_loc,ivpy)**2 &
              +fp(1:npar_loc,ivpz)**2),idiag_ekinp)
      endif
      if (idiag_vdotupm/=0) then
         call sum_par_name((fp(1:npar_loc,ivpx)*fp(1:npar_loc,iufx) +& 
              fp(1:npar_loc,ivpy)*fp(1:npar_loc,iufy) +&
              fp(1:npar_loc,ivpz)*fp(1:npar_loc,iufz) ), idiag_vdotupm) 
      endif
    endsubroutine single_particle_diagnos
!***********************************************************************
    subroutine dxxp_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle position (called from main pencil loop).
!
!  25-apr-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      integer :: k, ix0, iy0, iz0
      real :: dt1_advpx, dt1_advpy, dt1_advpz
!
!  Contribution of dust particles to time step.
!

!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(p)
!
    endsubroutine dxxp_dt_pencil
!***********************************************************************
    subroutine dvvp_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of dust particle velocity (called from main pencil loop).
!
!  Jan-2017, dhruba: dummy
!
      use Diagnostics
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ineargrid)
!
!  Diagnostic output.
!
!      if (ldiagnos) then
!
!      endif
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1
!
!      if (l1davgfirst) 

!
    endsubroutine dvvp_dt_pencil
!***********************************************************************
    subroutine dxxp_dt_blocks(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle position in blocks.
!
!  29-nov-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dxxp_dt_blocks
!***********************************************************************
    subroutine dvvp_dt_blocks(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle velocity in blocks.
!
!  29-nov-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dvvp_dt_blocks

!***********************************************************************
    subroutine read_particles_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_init_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_init_pars
!***********************************************************************
    subroutine write_particles_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_init_pars)
!
    endsubroutine write_particles_init_pars
!***********************************************************************
    subroutine read_particles_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_run_pars, IOSTAT=iostat)
!
!
    endsubroutine read_particles_run_pars
!***********************************************************************
    subroutine write_particles_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_run_pars)
!
    endsubroutine write_particles_run_pars
!***********************************************************************
    subroutine powersnap_particles(f)
!
!  Calculate power spectra of dust particle variables.
!
!  01-jan-06/anders: coded
!
      use Power_spectrum, only: power_1d
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!      if (lpar_spec) call power_1d(f,'p',0,irhop)
!
    endsubroutine powersnap_particles
!***********************************************************************
    subroutine rprint_particles(lreset,lwrite)
!
!  Read and register print parameters relevant for particles.
!
!  29-dec-04/anders: coded
!
      use Diagnostics
      use General,   only: itoa
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname,inamez,inamey,inamex,inamexy,inamexz,inamer,inamerz
      integer :: k
      logical :: lwr
      character (len=intlen) :: srad
!
!  Write information to index.pro.
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
      if (lwr) then
        write(3,*) 'ixp=', ixp
        write(3,*) 'iyp=', iyp
        write(3,*) 'izp=', izp
        write(3,*) 'ivpx=', ivpx
        write(3,*) 'ivpy=', ivpy
        write(3,*) 'ivpz=', ivpz
        write(3,*) 'inp=', inp
        write(3,*) 'irhop=', irhop
        write(3,*) 'iufx=', iufx
        write(3,*) 'iufy=', iufy
        write(3,*) 'iufz=', iufz
      endif
!
!  Reset everything in case of reset.
!
      if (lreset) then
        idiag_xpm=0; idiag_ypm=0; idiag_zpm=0
        idiag_vrelpabsm=0
        idiag_xp2m=0; idiag_yp2m=0; idiag_zp2m=0; idiag_rpm=0; idiag_rp2m=0
        idiag_vpxm=0; idiag_vpym=0; idiag_vpzm=0
        idiag_vpxvpym=0; idiag_vpxvpzm=0; idiag_vpyvpzm=0
        idiag_vpx2m=0; idiag_vpy2m=0; idiag_vpz2m=0; idiag_ekinp=0
        idiag_vpxmax=0; idiag_vpymax=0; idiag_vpzmax=0
        idiag_vpxmin=0; idiag_vpymin=0; idiag_vpzmin=0
        idiag_rhopvpxm=0; idiag_rhopvpym=0; idiag_rhopvpzm=0; idiag_rhopvpysm=0
        idiag_rhopvpxt=0; idiag_rhopvpyt=0; idiag_rhopvpzt=0
        idiag_lpxm=0; idiag_lpym=0; idiag_lpzm=0
        idiag_lpx2m=0; idiag_lpy2m=0; idiag_lpz2m=0
        idiag_npm=0; idiag_np2m=0; idiag_npmax=0; idiag_npmin=0
        idiag_dtdragp=0; idiag_dedragp=0
        idiag_rhopm=0; idiag_rhoprms=0; idiag_rhop2m=0; idiag_rhopmax=0
        idiag_rhopmin=0; idiag_decollp=0; idiag_rhopmphi=0
        idiag_epspmin=0; idiag_epspmax=0; idiag_epspm=0;
        idiag_nparmin=0; idiag_nparmax=0; idiag_nparsum=0
        idiag_nmigmax=0; idiag_nmigmmax=0; idiag_mpt=0
        idiag_npmx=0; idiag_npmy=0; idiag_npmz=0; idiag_epotpm=0
        idiag_rhopmx=0; idiag_rhopmy=0; idiag_rhopmz=0
        idiag_rhop2mx=0; idiag_rhop2my=0; idiag_rhop2mz=0
        idiag_epspmx=0; idiag_epspmy=0; idiag_epspmz=0
        idiag_rhopmxy=0; idiag_rhopmxz=0; idiag_rhopmr=0
        idiag_sigmap = 0
        idiag_dvpx2m=0; idiag_dvpy2m=0; idiag_dvpz2m=0
        idiag_dvpmax=0; idiag_dvpm=0; idiag_nparpmax=0
        idiag_eccpxm=0; idiag_eccpym=0; idiag_eccpzm=0
        idiag_eccpx2m=0; idiag_eccpy2m=0; idiag_eccpz2m=0
        idiag_npargone=0; idiag_vpyfull2m=0; idiag_deshearbcsm=0
        idiag_npmxy=0; idiag_vprms=0
        idiag_npvzmz=0; idiag_nptz=0; idiag_Shm=0
        idiag_npuzmz=0
        idiag_vdotupm=0
      endif
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot .and. ip<14) print*,'rprint_particles: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'nparsum',idiag_nparsum)
        call parse_name(iname,cname(iname),cform(iname),'nparmin',idiag_nparmin)
        call parse_name(iname,cname(iname),cform(iname),'nparmax',idiag_nparmax)
        call parse_name(iname,cname(iname),cform(iname),'nparpmax',idiag_nparpmax)
        call parse_name(iname,cname(iname),cform(iname),'xpm',idiag_xpm)
        call parse_name(iname,cname(iname),cform(iname),'vrelpabsm',idiag_vrelpabsm)
        call parse_name(iname,cname(iname),cform(iname),'ypm',idiag_ypm)
        call parse_name(iname,cname(iname),cform(iname),'zpm',idiag_zpm)
        call parse_name(iname,cname(iname),cform(iname),'xp2m',idiag_xp2m)
        call parse_name(iname,cname(iname),cform(iname),'yp2m',idiag_yp2m)
        call parse_name(iname,cname(iname),cform(iname),'zp2m',idiag_zp2m)
        call parse_name(iname,cname(iname),cform(iname),'rpm',idiag_rpm)
        call parse_name(iname,cname(iname),cform(iname),'rp2m',idiag_rp2m)
        call parse_name(iname,cname(iname),cform(iname),'vpxm',idiag_vpxm)
        call parse_name(iname,cname(iname),cform(iname),'vpym',idiag_vpym)
        call parse_name(iname,cname(iname),cform(iname),'vpzm',idiag_vpzm)
        call parse_name(iname,cname(iname),cform(iname),'vpxvpym',idiag_vpxvpym)
        call parse_name(iname,cname(iname),cform(iname),'vpxvpzm',idiag_vpxvpzm)
        call parse_name(iname,cname(iname),cform(iname),'vpyvpzm',idiag_vpyvpzm)
        call parse_name(iname,cname(iname),cform(iname),'vpx2m',idiag_vpx2m)
        call parse_name(iname,cname(iname),cform(iname),'vpy2m',idiag_vpy2m)
        call parse_name(iname,cname(iname),cform(iname),'vpz2m',idiag_vpz2m)
        call parse_name(iname,cname(iname),cform(iname),'ekinp',idiag_ekinp)
        call parse_name(iname,cname(iname),cform(iname),'vpxmax',idiag_vpxmax)
        call parse_name(iname,cname(iname),cform(iname),'vpymax',idiag_vpymax)
        call parse_name(iname,cname(iname),cform(iname),'vpzmax',idiag_vpzmax)
        call parse_name(iname,cname(iname),cform(iname),'vpxmin',idiag_vpxmin)
        call parse_name(iname,cname(iname),cform(iname),'vpymin',idiag_vpymin)
        call parse_name(iname,cname(iname),cform(iname),'vpzmin',idiag_vpzmin)
        call parse_name(iname,cname(iname),cform(iname),'vpmax',idiag_vpmax)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpxm',idiag_rhopvpxm)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpym',idiag_rhopvpym)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpzm',idiag_rhopvpzm)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpysm',idiag_rhopvpysm)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpxt',idiag_rhopvpxt)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpyt',idiag_rhopvpyt)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpzt',idiag_rhopvpzt)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpysm',idiag_rhopvpysm)
        call parse_name(iname,cname(iname),cform(iname),'lpxm',idiag_lpxm)
        call parse_name(iname,cname(iname),cform(iname),'lpym',idiag_lpym)
        call parse_name(iname,cname(iname),cform(iname),'lpzm',idiag_lpzm)
        call parse_name(iname,cname(iname),cform(iname),'lpx2m',idiag_lpx2m)
        call parse_name(iname,cname(iname),cform(iname),'lpy2m',idiag_lpy2m)
        call parse_name(iname,cname(iname),cform(iname),'lpz2m',idiag_lpz2m)
        call parse_name(iname,cname(iname),cform(iname),'eccpxm',idiag_eccpxm)
        call parse_name(iname,cname(iname),cform(iname),'eccpym',idiag_eccpym)
        call parse_name(iname,cname(iname),cform(iname),'eccpzm',idiag_eccpzm)
        call parse_name(iname,cname(iname),cform(iname),'eccpx2m',idiag_eccpx2m)
        call parse_name(iname,cname(iname),cform(iname),'eccpy2m',idiag_eccpy2m)
        call parse_name(iname,cname(iname),cform(iname),'eccpz2m',idiag_eccpz2m)
        call parse_name(iname,cname(iname),cform(iname),'dtdragp',idiag_dtdragp)
        call parse_name(iname,cname(iname),cform(iname),'npm',idiag_npm)
        call parse_name(iname,cname(iname),cform(iname),'np2m',idiag_np2m)
        call parse_name(iname,cname(iname),cform(iname),'npmax',idiag_npmax)
        call parse_name(iname,cname(iname),cform(iname),'npmin',idiag_npmin)
        call parse_name(iname,cname(iname),cform(iname),'rhopm',idiag_rhopm)
        call parse_name(iname,cname(iname),cform(iname),'rhoprms',idiag_rhoprms)
        call parse_name(iname,cname(iname),cform(iname),'rhop2m',idiag_rhop2m)
        call parse_name(iname,cname(iname),cform(iname),'rhopmin',idiag_rhopmin)
        call parse_name(iname,cname(iname),cform(iname),'rhopmax',idiag_rhopmax)
        call parse_name(iname,cname(iname),cform(iname),'epspm',idiag_epspm)
        call parse_name(iname,cname(iname),cform(iname),'epspmin',idiag_epspmin)
        call parse_name(iname,cname(iname),cform(iname),'epspmax',idiag_epspmax)
        call parse_name(iname,cname(iname),cform(iname),'rhopmphi',idiag_rhopmphi)
        call parse_name(iname,cname(iname),cform(iname),'nmigmax',idiag_nmigmax)
        call parse_name(iname,cname(iname),cform(iname),'nmigmmax',idiag_nmigmmax)
        call parse_name(iname,cname(iname),cform(iname),'mpt',idiag_mpt)
        call parse_name(iname,cname(iname),cform(iname),'dvpx2m',idiag_dvpx2m)
        call parse_name(iname,cname(iname),cform(iname),'dvpy2m',idiag_dvpy2m)
        call parse_name(iname,cname(iname),cform(iname),'dvpz2m',idiag_dvpz2m)
        call parse_name(iname,cname(iname),cform(iname),'dvpm',idiag_dvpm)
        call parse_name(iname,cname(iname),cform(iname),'dvpmax',idiag_dvpmax)
        call parse_name(iname,cname(iname),cform(iname),'dedragp',idiag_dedragp)
        call parse_name(iname,cname(iname),cform(iname),'decollp',idiag_decollp)
        call parse_name(iname,cname(iname),cform(iname),'epotpm',idiag_epotpm)
        call parse_name(iname,cname(iname),cform(iname),'npargone',idiag_npargone)
        call parse_name(iname,cname(iname),cform(iname),'vpyfull2m',idiag_vpyfull2m)
        call parse_name(iname,cname(iname),cform(iname),'vprms',idiag_vprms)
        call parse_name(iname,cname(iname),cform(iname),'Shm',idiag_Shm)
        call parse_name(iname,cname(iname),cform(iname),'deshearbcsm',idiag_deshearbcsm)
        call parse_name(iname,cname(iname),cform(iname),'vdotupm',idiag_vdotupm)
      enddo
!
!  Check for those quantities for which we want x-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'npmx',idiag_npmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'rhopmx',idiag_rhopmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'rhop2mx',idiag_rhop2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'epspmx',idiag_epspmx)
      enddo
!
!  Check for those quantities for which we want y-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'npmy',idiag_npmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'rhopmy',idiag_rhopmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'rhop2my',idiag_rhop2my)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'epspmy',idiag_epspmy)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'npmz',idiag_npmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rhopmz',idiag_rhopmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rhop2mz',idiag_rhop2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'epspmz',idiag_epspmz)
        do k=1,ninit
          srad=itoa(k)
          call parse_name(inamez,cnamez(inamez),cformz(inamez),'npvzmz'//trim(srad),idiag_npvzmz(k))
          call parse_name(inamez,cnamez(inamez),cformz(inamez),'npuzmz'//trim(srad),idiag_npuzmz(k))
          call parse_name(inamez,cnamez(inamez),cformz(inamez),'nptz'//trim(srad),idiag_nptz(k))
        enddo

      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'npmxy',idiag_npmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'rhopmxy',idiag_rhopmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'sigmap',idiag_sigmap)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamexz=1,nnamexz
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'rhopmxz',idiag_rhopmxz)
      enddo
!
!  Check for those quantities for which we want phiz-averages.
!
      do inamer=1,nnamer
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'rhopmr',idiag_rhopmr)
      enddo
!
!  Check for those quantities for which we want phi-averages.
!
      do inamerz=1,nnamerz
        call parse_name(inamerz,cnamerz(inamerz),cformrz(inamerz),'rhopmphi',idiag_rhopmphi)
      enddo
!
    endsubroutine rprint_particles
!***********************************************************************
    subroutine particles_final_clean_up
!
!  cleanup (dummy)
!


    endsubroutine particles_final_clean_up
!***********************************************************************
    subroutine periodic_boundcond_on_aux(f)
!
!
! Impose periodic boundary condition on gradu as auxiliary variable
!
      use Boundcond, only: set_periodic_boundcond_on_aux
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f

!
      if (lparticles_grad) then
        if (igradu .ne. 0) then
          call set_periodic_boundcond_on_aux(f,igradu11)
          call set_periodic_boundcond_on_aux(f,igradu12)
          call set_periodic_boundcond_on_aux(f,igradu13)
          call set_periodic_boundcond_on_aux(f,igradu21)
          call set_periodic_boundcond_on_aux(f,igradu22)
          call set_periodic_boundcond_on_aux(f,igradu23)
          call set_periodic_boundcond_on_aux(f,igradu31)
          call set_periodic_boundcond_on_aux(f,igradu32)
          call set_periodic_boundcond_on_aux(f,igradu33)
        else
          call fatal_error('periodic_boundcond_on_aux','particles_grad demands igradu ne 0')
        endif
      endif
!
    endsubroutine periodic_boundcond_on_aux
!***********************************************************************
    subroutine calc_relative_velocity(f,fp,ineargrid)
!
      use Diagnostics
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      real, dimension(3) :: uup,rel_vel_sing
      real, dimension(:), allocatable :: rel_vel
      integer :: k,ix0,iy0,iz0
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
!
!  Calculate particle relative velocity
!
      allocate(rel_vel(npar_loc))
!
      rel_vel = 0.0
!
      do k = 1,npar_loc
        call interpolate_linear(f,iux,iuz,fp(k,ixp:izp),uup,ineargrid(k,:),0,0)
        rel_vel_sing = (fp(k,ivpx:ivpz)-uup)**2
        rel_vel(k) = sqrt(sum(rel_vel_sing))
      enddo
!
      call sum_par_name(rel_vel(1:npar_loc),idiag_vrelpabsm)
      if (allocated(rel_vel)) deallocate(rel_vel)
!
    endsubroutine calc_relative_velocity
!***********************************************************************
    subroutine remove_particles_sink_simple(f,fp,dfp,ineargrid)
!
      real,    dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3)       :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine remove_particles_sink_simple
!***********************************************************************
    subroutine create_particles_sink_simple(f,fp,dfp,ineargrid)
!
      real,    dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3)       :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine create_particles_sink_simple
!***********************************************************************
    subroutine insert_particles(f,fp,ineargrid)
!
! Insert particles continuously (when linsert_particles_continuously == T),
! i.e. in each timestep. If number of particles to be inserted are less
! than unity, accumulate number over several timesteps until the integer value
! is larger than one. Keep the remainder and accumulate this to the next insert.
!
! Works only for particles_dust - add neccessary variable
! declarations in particles_tracers to make it work here.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3)    :: ineargrid
!
      intent (inout) :: fp,ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine insert_particles
!***********************************************************************
!***********************************************************************
    subroutine particles_dragforce_stiff(f,fp,ineargrid)
!
!  10-june-11/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine particles_dragforce_stiff
!***********************************************************************

    
endmodule Particles
