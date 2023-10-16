! $Id$
!  This module takes care of everything related to tracer particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles=.true.
! CPARAM character (len=20), parameter :: particles_module="tracer"
!
! MAUX CONTRIBUTION 2
! MPVAR CONTRIBUTION 3
! MPAUX CONTRIBUTION 3
!
! PENCILS PROVIDED np; rhop; epsp; rhop_swarm; grhop(3); peh
!
!***************************************************************
module Particles
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_sub
  use Particles_mpicomm
!
  implicit none
!
  include 'particles.h'
  include 'particles_common.h'
!
  real :: xp0=0.0, yp0=0.0, zp0=0.0, tausp=0.0
  real :: nu_epicycle=0.0, nu_epicycle2=0.0
  logical :: ldragforce_equi_global_eps=.false.
  logical :: lquadratic_interpolation=.false.
  logical :: ltrace_dust=.false.
  real :: particles_insert_rate=0.,max_particle_insert_time=huge1
  character (len=labellen), dimension (ninit) :: initxxp='nothing'
  character (len=labellen), dimension (ninit) :: insertxxp='nothing'
  character (len=labellen) :: gravz_profile='zero'
  logical :: lglobalrandom=.false.
!
  namelist /particles_init_pars/ &
      initxxp, xp0, yp0, zp0, bcpx, bcpy, bcpz, eps_dtog, tausp, &
      ldragforce_equi_global_eps, lquadratic_interpolation, &
      lparticlemesh_cic, lparticlemesh_tsc, ltrace_dust, &
      gravz_profile, nu_epicycle, lglobalrandom
!
  namelist /particles_run_pars/ &
      bcpx, bcpy, bcpz, lquadratic_interpolation, &
      lparticlemesh_cic, lparticlemesh_tsc, ltrace_dust, &
      lcheck_exact_frontier,particles_insert_rate, &
      linsert_particles_continuously,dsnap_par, &
      max_particle_insert_time,insertxxp
!
  integer :: idiag_xpm=0, idiag_ypm=0, idiag_zpm=0
  integer :: idiag_xp2m=0, idiag_yp2m=0, idiag_zp2m=0
  integer :: idiag_nparmax=0, idiag_npmax=0, idiag_npmin=0
  integer :: idiag_npmx=0, idiag_rhopmx=0, idiag_epspmx=0
  integer :: idiag_npmz=0, idiag_rhopmz=0, idiag_epspmz=0
  integer :: idiag_epsmin=0,idiag_epsmax=0
!
  contains
!***********************************************************************
    subroutine register_particles()
!
!  Set up indices for access to the fp and dfp arrays
!
!  29-dec-04/anders: coded
!
      use FArrayManager
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  Indices for particle position.
!
      call append_npvar('ixp',ixp)
      call append_npvar('iyp',iyp)
      call append_npvar('izp',izp)
!
!  Indices for particle velocity.
!
      call append_npaux('ivpx',ivpx)
      call append_npaux('ivpy',ivpy)
      call append_npaux('ivpz',ivpz)
!
!  Set indices for auxiliary variables.
!
      call farray_register_auxiliary('np',inp)
      call farray_register_auxiliary('rhop',irhop)
!
    endsubroutine register_particles
!***********************************************************************
    subroutine initialize_particles(f,fp)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  29-dec-04/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray), intent (in) :: fp
!
      real :: rhom
!
      if (rhop_swarm==0.0 .or. mp_swarm==0.0) then
!
! For stratification, take into account gas present outside the simulation box.
!
        if ((lgravz .and. lgravz_gas) .or. gravz_profile=='linear' ) then
          rhom=sqrt(2*pi)*1.0*1.0/Lz  ! rhom = Sigma/Lz, Sigma=sqrt(2*pi)*H*rho1
        else
          rhom=1.0
        endif
        if (rhop_swarm==0.0) rhop_swarm = eps_dtog*rhom/(real(npar)/(nxgrid*nygrid*nzgrid))
        if (mp_swarm==0.0) mp_swarm = eps_dtog*rhom*box_volume/(real(npar))
        if (lroot) print*, 'initialize_particles: dust-to-gas ratio eps_dtog=', eps_dtog
      endif
!
      if (lroot) then
        print*, 'initialize_particles: mass per constituent particle mpmat=', mpmat
        print*, 'initialize_particles: mass per superparticle mp_swarm =', mp_swarm
        print*, 'initialize_particles: number density per superparticle np_swarm=', np_swarm
        print*, 'initialize_particles: mass density per superparticle rhop_swarm=',rhop_swarm
      endif
!
!  Calculate nu_epicycle**2 for gravity.
!
      if (nu_epicycle/=0.0) then
        gravz_profile='linear'
        nu_epicycle2=nu_epicycle**2
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles
!***********************************************************************
    subroutine init_particles(f,fp,ineargrid)
!
!  Initial positions and velocities of tracer particles.
!
!  29-dec-04/anders: coded
!
      use EquationOfState, only: cs20
      use General, only: random_number_wrapper
      use SharedVariables, only: get_shared_variable
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (inout) :: f
      intent (out) :: fp
!
      real, dimension (nx) :: eps
      real, dimension (3) :: Lxyz_par, xyz0_par, xyz1_par
      real :: r, p, cs
      integer :: j, k
      logical :: lnothing
!
      real :: rad,rad_scl,phi,tmp
!
      real, dimension(:), pointer :: beta_glnrho_global
!
      call get_shared_variable('beta_glnrho_global',beta_glnrho_global,caller='init_particles')

!  Use either a local random position or a global random position for certain
!  initial conditions. The default is a local random position, but the equal
!  number of particles per processors means that this is not completely random.
!
      if (lglobalrandom) then
        Lxyz_par=Lxyz
        xyz0_par=xyz0
        xyz1_par=xyz1
      else
        Lxyz_par=Lxyz_loc
        xyz0_par=xyz0_loc
        xyz1_par=xyz1_loc
      endif
!
!  Initial particle position.
!
      lnothing=.false.
      do j=1,ninit
!
        select case (initxxp(j))
!
        case ('nothing')
          if (lroot .and. .not. lnothing) print*, 'init_particles: nothing'
          lnothing=.true.
!
        case ('origin')
          if (lroot) print*, 'init_particles: All particles at origin'
          fp(1:npar_loc,ixp:izp)=0.
!
        case ('constant')
          if (lroot) print*, 'init_particles: All particles at x,y,z=', xp0, yp0, zp0
          fp(1:npar_loc,ixp)=xp0
          fp(1:npar_loc,iyp)=yp0
          fp(1:npar_loc,izp)=zp0
!
        case ('random')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k=1,npar_loc
            if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
            if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
            if (nzgrid/=1) call random_number_wrapper(fp(k,izp))
          enddo
          if (nxgrid/=1) fp(1:npar_loc,ixp)=xyz0_loc(1)+fp(1:npar_loc,ixp)*Lxyz_loc(1)
          if (nygrid/=1) fp(1:npar_loc,iyp)=xyz0_loc(2)+fp(1:npar_loc,iyp)*Lxyz_loc(2)
          if (nzgrid/=1) fp(1:npar_loc,izp)=xyz0_loc(3)+fp(1:npar_loc,izp)*Lxyz_loc(3)
!
        case ('random-xy')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k=1,npar_loc
            if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
            if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
          enddo
          if (nxgrid/=1) fp(1:npar_loc,ixp)=xyz0_loc(1)+fp(1:npar_loc,ixp)*Lxyz_loc(1)
          if (nygrid/=1) fp(1:npar_loc,iyp)=xyz0_loc(2)+fp(1:npar_loc,iyp)*Lxyz_loc(2)
          fp(1:npar_loc,izp)=zp0
!
        case ('random-xz')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k=1,npar_loc
            if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
            if (nzgrid/=1) call random_number_wrapper(fp(k,izp))
          enddo
          if (nxgrid/=1) fp(1:npar_loc,ixp)=xyz0_loc(1)+fp(1:npar_loc,ixp)*Lxyz_loc(1)
          if (nzgrid/=1) fp(1:npar_loc,izp)=xyz0_loc(3)+fp(1:npar_loc,izp)*Lxyz_loc(3)
          fp(1:npar_loc,iyp)=yp0
!
       case ('random-cylindrical','random-cyl')
          if (lroot) print*, 'init_particles: Random particle '//&
               'cylindrical positions with power-law =',dustdensity_powerlaw
!
          do k=1,npar_loc
!
! Start the particles obeying a power law
!
            tmp=2-dustdensity_powerlaw
            call random_number_wrapper(rad_scl)
            rad_scl = rp_int**tmp + rad_scl*(rp_ext**tmp-rp_int**tmp)
            rad = rad_scl**(1./tmp)
!
! Random in azimuth
!
            call random_number_wrapper(phi)
!
             if (lcartesian_coords) then
               phi = 2*pi*phi
               if (nxgrid/=1) fp(k,ixp)=rad*cos(phi)
               if (nygrid/=1) fp(k,iyp)=rad*sin(phi)
             elseif (lcylindrical_coords) then
               phi = xyz0_par(2)+phi*Lxyz_par(2)
               if (nxgrid/=1) fp(k,ixp)=rad
               if (nygrid/=1) fp(k,iyp)=phi
             elseif (lspherical_coords) then
               call not_implemented("init_particles','random-cylindrical for spherical coordinates")
             endif
!
             if (nzgrid/=1) call random_number_wrapper(fp(k,izp))
             if (nzgrid/=1) fp(k,izp)=xyz0_par(3)+fp(k,izp)*Lxyz_par(3)
!
          enddo
!
        case ('gaussian-z')
          if (lroot) print*, 'init_particles: Gaussian particle positions'
          do k=1,npar_loc
            if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
            if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
            call random_number_wrapper(r)
            call random_number_wrapper(p)
            if (nprocz==2) then
              if (lfirst_proc_z) fp(k,izp)=-abs(zp0*sqrt(-2*alog(r))*cos(2*pi*p))
              if (llast_proc_z) fp(k,izp)=abs(zp0*sqrt(-2*alog(r))*cos(2*pi*p))
            else
              fp(k,izp)=zp0*sqrt(-2*alog(r))*cos(2*pi*p)
            endif
          enddo
          if (nxgrid/=1) fp(1:npar_loc,ixp)=xyz0_loc(1)+fp(1:npar_loc,ixp)*Lxyz_loc(1)
          if (nygrid/=1) fp(1:npar_loc,iyp)=xyz0_loc(2)+fp(1:npar_loc,iyp)*Lxyz_loc(2)
!
        case ('dragforce_equilibrium')
!
!  Equilibrium between drag forces on dust and gas and other forces
!  (from Nakagawa, Sekiya, & Hayashi 1986).
!
          if (lroot) then
            print*, 'init_particles: drag equilibrium'
            print*, 'init_particles: beta_glnrho_global=', beta_glnrho_global
          endif
!
          if (.not.(lcartesian_coords.and.(all(lequidist)))) &
               call not_implemented("init_particles","dragforce_equilibrium " // &
                                    "initial condition for polar or non-equidistant grids")
!  Calculate average dust-to-gas ratio in box.

          if (ldensity_nolog) then
            eps = rhop_swarm*sum(f(l1:l2,m1:m2,n1:n2,inp))/sum(f(l1:l2,m1:m2,n1:n2,irho))
          else
            eps = rhop_swarm*sum(f(l1:l2,m1:m2,n1:n2,inp))/sum(exp(f(l1:l2,m1:m2,n1:n2,ilnrho)))
          endif
!
          if (lroot) print*, 'init_particles: average dust-to-gas ratio=', eps(1)
!  Set gas velocity field.
          do m=m1,m2; do n=n1,n2
            cs=sqrt(cs20)
!  Take either global or local dust-to-gas ratio.
            if (.not. ldragforce_equi_global_eps) then
              if (ldensity_nolog) then
                eps = rhop_swarm*f(l1:l2,m,n,inp)/f(l1:l2,m,n,irho)
              else
                eps = rhop_swarm*f(l1:l2,m,n,inp)/exp(f(l1:l2,m,n,ilnrho))
              endif
            endif
!
            f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + beta_glnrho_global(1)/(2*(1.0+eps))*cs
!
          enddo; enddo
!
        case default
          call fatal_error("init_particles",'no such initxxp: '//trim(initxxp(j)))
!
        endselect
!
      enddo
!
!  Particles are not allowed to be present in non-existing dimensions.
!  This would give huge problems with interpolation later.
!
      if (nxgrid==1) fp(1:npar_loc,ixp)=x(l1)
      if (nygrid==1) fp(1:npar_loc,iyp)=y(m1)
      if (nzgrid==1) fp(1:npar_loc,izp)=z(n1)
!
!  Redistribute particles among processors (now that positions are determined).
!
      call boundconds_particles(fp,ipar)
!
!  Map particle positions on the grid.
!
      call map_nearest_grid(fp,ineargrid)
      call map_xxp_grid(f,fp,ineargrid)
!
!  Sort particles so that they can be accessed contiguously in the memory.
!
      call sort_particles_imn(fp,ineargrid,ipar)
!
    endsubroutine init_particles
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
    subroutine pencil_criteria_particles()
!
!  All pencils that the Particles module depends on are specified here.
!
!  20-04-06/anders: coded
!
      lpenc_diagnos(i_np)=.true.
!
      if (idiag_epsmin/=0.or.idiag_epsmax/=0) lpenc_diagnos(i_epsp)=.true.
!
    endsubroutine pencil_criteria_particles
!***********************************************************************
    subroutine pencil_interdep_particles(lpencil_in)
!
!  Interdependency among pencils provided by the Particles module
!  is specified here.
!
!  15-feb-06/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_rhop) .and. irhop==0) then
        lpencil_in(i_np)=.true.
      endif
!
      if (lpencil_in(i_epsp)) then
        lpencil_in(i_rho)=.true.
        lpencil_in(i_rhop)=.true.
      endif
!
    endsubroutine pencil_interdep_particles
!***********************************************************************
    subroutine calc_pencils_particles(f,p)
!
!  Calculate particle pencils.
!
!  15-feb-06/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
! np
      if (lpencil(i_np)) p%np=f(l1:l2,m,n,inp)
! rhop
      if (lpencil(i_rhop)) then
        if (irhop/=0) then
          p%rhop=f(l1:l2,m,n,irhop)
        else
          p%rhop=rhop_swarm*p%np
        endif
      endif
! epsp
      if (lpencil(i_epsp)) p%epsp=p%rhop/p%rho
!
    endsubroutine calc_pencils_particles
!***********************************************************************
    subroutine dxxp_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of tracer particle position.
!
!  02-jan-05/anders: coded
!
      use Diagnostics, only: max_name
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (in) :: f, fp, ineargrid
      intent (inout) :: dfp, df
!
      logical :: lheader, lfirstcall=.true.
!
      lheader=lfirstcall .and. lroot
!
!  Vertical gravity in the short friction time approximation.
!
      select case (gravz_profile)
!
        case ('zero')
          if (lheader) print*, 'dxxp_dt: No gravity in z-direction.'
!
        case ('linear')
          if (lheader) print*, 'dxxp_dt: Linear gravity field in z-direction.'
          dfp(1:npar_loc,izp)=dfp(1:npar_loc,izp) - tausp*nu_epicycle2*fp(1:npar_loc,izp)
!
      endselect
!
!  Diagnostic output
!
      if (ldiagnos) then
        call sum_par_name(fp(1:npar_loc,ixp),idiag_xpm)
        call sum_par_name(fp(1:npar_loc,iyp),idiag_ypm)
        call sum_par_name(fp(1:npar_loc,izp),idiag_zpm)
        if (idiag_xp2m/=0) call sum_par_name(fp(1:npar_loc,ixp)**2,idiag_xp2m)
        if (idiag_yp2m/=0) call sum_par_name(fp(1:npar_loc,iyp)**2,idiag_yp2m)
        if (idiag_zp2m/=0) call sum_par_name(fp(1:npar_loc,izp)**2,idiag_zp2m)
        call max_name(npar_loc,idiag_nparmax)
      endif
!
      if (lfirstcall) lfirstcall=.false.
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dxxp_dt
!***********************************************************************
    subroutine dvvp_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of dust particle velocity.
!
!  22-aug-05/anders: dummy
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
    endsubroutine dvvp_dt
!***********************************************************************
    subroutine dxxp_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of tracer particle position (called from main pencil loop).
!
!  25-apr-06/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: uu
      integer :: k
!
      intent (in) :: f, ineargrid
      intent (inout) :: df, dfp, fp
!
!  Identify module and boundary conditions.
!
      if (headtt) then
        print*, 'dxxp_dt: Calculate dxxp_dt'
        print*, 'dxxp_dt: Particles boundary condition bcpx=', bcpx
        print*, 'dxxp_dt: Particles boundary condition bcpy=', bcpy
        print*, 'dxxp_dt: Particles boundary condition bcpz=', bcpz
        print*, 'dxxp_dt: Set rate of change of particle position equal to gas velocity.'
      endif
!
!  Interpolate gas velocity to position of particles.
!  Then set particle velocity equal to the local gas velocity.
!
      if (npar_imn(imn)/=0) then
        do k=k1_imn(imn),k2_imn(imn)
          if (lparticlemesh_tsc) then
            if (ltrace_dust) then
              call interpolate_quadratic_spline(f,iudx(1),iudz(1),fp(k,ixp:izp),uu,ineargrid(k,:),0,ipar(k))
            else
              call interpolate_quadratic_spline(f,iux,iuz,fp(k,ixp:izp),uu,ineargrid(k,:),0,ipar(k))
            endif
          else
            if (ltrace_dust) then
              call interpolate_linear(f,iudx(1),iudz(1),fp(k,ixp:izp),uu,ineargrid(k,:),0,ipar(k))
            else
              call interpolate_linear(f,iux,iuz,fp(k,ixp:izp),uu,ineargrid(k,:),0,ipar(k))
            endif
          endif
!
!  Advance of particle position
!
          if (lcartesian_coords) then
            if (nxgrid/=1) dfp(k,ixp) = dfp(k,ixp) + uu(1) ; fp(k,ivpx) = uu(1)
            if (nygrid/=1) dfp(k,iyp) = dfp(k,iyp) + uu(2) ; fp(k,ivpy) = uu(2)
            if (nzgrid/=1) dfp(k,izp) = dfp(k,izp) + uu(3) ; fp(k,ivpz) = uu(3)
                          
          elseif (lcylindrical_coords) then
            if (nxgrid/=1) dfp(k,ixp) = dfp(k,ixp) + uu(1)
            if (nygrid/=1) dfp(k,iyp) = dfp(k,iyp) + uu(2)/max(fp(k,ixp),tini)
            if (nzgrid/=1) dfp(k,izp) = dfp(k,izp) + uu(3)
          elseif (lspherical_coords) then
            if (nxgrid/=1) dfp(k,ixp) = dfp(k,ixp) + uu(1)
            if (nygrid/=1) dfp(k,iyp) = dfp(k,iyp) + uu(2)/max(fp(k,ixp),tini)
            if (nzgrid/=1) dfp(k,izp) = dfp(k,izp) + uu(3)/(max(fp(k,ixp),tini)*sin(fp(k,iyp)))
          endif
!
!  With shear there is an extra term due to the background shear flow.
!
          if (lshear.and.nygrid/=1) dfp(k,iyp) = dfp(k,iyp) - qshear*Omega*fp(k,ixp)
!
        enddo
      endif
!
      call calc_diagnostics_particles(fp,p,ineargrid)
!
    endsubroutine dxxp_dt_pencil
!***********************************************************************
    subroutine calc_diagnostics_particles(fp,p,ineargrid)

      use Diagnostics

      real, dimension (mpar_loc,mparray) :: fp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      if (ldiagnos) then
        call max_mn_name(p%np,idiag_npmax)
        if (idiag_npmin/=0) call max_mn_name(-p%np,idiag_npmin,lneg=.true.)
        call max_mn_name(p%epsp,idiag_epsmax)
        if (idiag_epsmin/=0) call max_mn_name(-p%epsp,idiag_epsmin,lneg=.true.)
      endif
!
      if (l1davgfirst) then
        call yzsum_mn_name_x(p%np,idiag_npmx)
        call yzsum_mn_name_x(p%rhop,idiag_rhopmx)
        call yzsum_mn_name_x(p%epsp,idiag_epspmx)
        call xysum_mn_name_z(p%np,idiag_npmz)
        call xysum_mn_name_z(p%rhop,idiag_rhopmz)
        call xysum_mn_name_z(p%epsp,idiag_epspmz)
      endif
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine calc_diagnostics_particles
!***********************************************************************
    subroutine dvvp_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of dust particle velocity (called from main pencil loop).
!
!  20-apr-06/anders: dummy
!
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
      real, dimension (mx,my,mz,mvar)    :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar)   :: dfp
      integer, dimension (mpar_loc,3)    :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dvvp_dt_blocks
!***********************************************************************
    subroutine remove_particles_sink_simple(f,fp,dfp,ineargrid)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar)   :: dfp
      integer, dimension (mpar_loc,3)    :: ineargrid
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
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar)   :: dfp
      integer, dimension (mpar_loc,3)    :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine create_particles_sink_simple
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
!  Calculate power spectra of particle variables.
!
!  01-jan-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine powersnap_particles
!***********************************************************************
    subroutine insert_lost_particles(f,fp,ineargrid)
!
!  14-oct-12/dhruba: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (inout) :: fp,ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine insert_lost_particles
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
!
      use General, only: random_number_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3)    :: ineargrid
      logical :: linsertmore=.true.
      real :: avg_n_insert,remaining_particles
!
      integer :: j, k, n_insert, npar_loc_old, iii
!
      intent (inout) :: fp,ineargrid
!
! Stop call to this routine when maximum number of particles is reached!
! Since root inserts all new particles, make sure
! npar_total + n_insert < mpar
! so that a processor can not exceed its maximum number of particles.
!
      if (lroot) then
        avg_n_insert=particles_insert_rate*dt
        n_insert=int(avg_n_insert + remaining_particles)
!
! Remaining particles saved for subsequent timestep:
!
        remaining_particles=avg_n_insert + remaining_particles - n_insert
        if ((n_insert+npar_total <= mpar_loc) .and. (t<max_particle_insert_time)) then
          linsertmore=.true.
        else
          linsertmore=.false.
        endif
!
        if (linsertmore) then
! Actual (integer) number of particles to be inserted at this timestep:
          do iii=npar_loc+1,npar_loc+n_insert
            ipar(iii)=npar_total+iii-npar_loc
          enddo
          npar_total=npar_total+n_insert
          npar_loc_old=npar_loc
          npar_loc=npar_loc + n_insert
!
! Insert particles in chosen position (as in init_particles).
!
          do j=1,ninit
            select case (insertxxp(j))
            case ('nothing')
              if (lroot.and.ip<10) print*, 'init_particles: nothing'
            case ('random-xy')
              if (lroot.and.ip<10) print*, 'init_particles: Random particle positions'
              do k=1,npar_loc
                if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
                if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
              enddo
              if (nxgrid/=1) fp(1:npar_loc,ixp)=xyz0_loc(1)+fp(1:npar_loc,ixp)*Lxyz_loc(1)
              if (nygrid/=1) fp(1:npar_loc,iyp)=xyz0_loc(2)+fp(1:npar_loc,iyp)*Lxyz_loc(2)
              fp(1:npar_loc,izp)=zp0
              !
            case default
              call fatal_error("particles_tracers","no such insertxxp: "//trim(insertxxp(j)))
              !
            endselect
!
          enddo ! do j=1,ninit
          if (nxgrid==1) fp(npar_loc_old+1:npar_loc,ixp)=x(nghost+1)
          if (nygrid==1) fp(npar_loc_old+1:npar_loc,iyp)=y(nghost+1)
          if (nzgrid==1) fp(npar_loc_old+1:npar_loc,izp)=z(nghost+1)
        endif
      endif
!
!  Redistribute particles only when t < max_particle_insert_time.
!  Could have included some other tests here aswell......
!
      if (t<max_particle_insert_time) then
!
!  Redistribute particles among processors.
!
        call boundconds_particles(fp,ipar,linsert=.true.)
!
!  Map particle positions on the grid.
!
        call map_nearest_grid(fp,ineargrid)
        call map_xxp_grid(f,fp,ineargrid)
!
!  Sort particles so that they can be accessed contiguously in the memory.
!
        call sort_particles_imn(fp,ineargrid,ipar)
      endif
!
    endsubroutine insert_particles
!***********************************************************************
    subroutine rprint_particles(lreset,lwrite)
!
!  Read and register print parameters relevant for particles
!
!  29-dec-04/anders: coded
!
      use Diagnostics, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname, inamex, inamez
!
!  Reset everything in case of reset
!
      if (lreset) then
        idiag_xpm=0; idiag_ypm=0; idiag_zpm=0
        idiag_xp2m=0; idiag_yp2m=0; idiag_zp2m=0
        idiag_nparmax=0; idiag_nmigmax=0; idiag_npmax=0; idiag_npmin=0
        idiag_npmx=0; idiag_rhopmx=0; idiag_epspmx=0
        idiag_npmz=0; idiag_rhopmz=0; idiag_epspmz=0
        idiag_epsmin=0; idiag_epsmax=0
      endif
!
!  Run through all possible names that may be listed in print.in
!
      if (lroot .and. ip<14) print*,'rprint_particles: run through parse list'
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'xpm',idiag_xpm)
        call parse_name(iname,cname(iname),cform(iname),'ypm',idiag_ypm)
        call parse_name(iname,cname(iname),cform(iname),'zpm',idiag_zpm)
        call parse_name(iname,cname(iname),cform(iname),'xp2m',idiag_xp2m)
        call parse_name(iname,cname(iname),cform(iname),'yp2m',idiag_yp2m)
        call parse_name(iname,cname(iname),cform(iname),'zp2m',idiag_zp2m)
        call parse_name(iname,cname(iname),cform(iname),'nparmax',idiag_nparmax)
        call parse_name(iname,cname(iname),cform(iname),'npmax',idiag_npmax)
        call parse_name(iname,cname(iname),cform(iname),'npmin',idiag_npmin)
        call parse_name(iname,cname(iname),cform(iname),'nmigmax',idiag_nmigmax)
        call parse_name(iname,cname(iname),cform(iname),'epsmax',idiag_epsmax)
        call parse_name(iname,cname(iname),cform(iname),'epsmin',idiag_epsmin)
      enddo
!
!  check for those quantities for which we want x-averages
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'npmx',idiag_npmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'rhopmx',idiag_rhopmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'epspmx',idiag_epspmx)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'npmz',idiag_npmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rhopmz',idiag_rhopmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'epspmz',idiag_epspmz)
      enddo
!
    endsubroutine rprint_particles
!***********************************************************************
    subroutine particles_final_clean_up
!
    endsubroutine particles_final_clean_up
!***********************************************************************
    subroutine periodic_boundcond_on_aux(f)
!
! dummy
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine periodic_boundcond_on_aux
!***********************************************************************
endmodule Particles
