! $Id: particles_dust.f90,v 1.120 2006-07-19 20:46:16 ajohan Exp $
!
!  This module takes care of everything related to dust particles
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 6
! MAUX CONTRIBUTION 2
! CPARAM logical, parameter :: lparticles=.true.
! CPARAM logical, parameter :: lparticles_nbody=.false.
!
! PENCILS PROVIDED np, rhop
!
!***************************************************************
module Particles

  use Cdata
  use Particles_cdata
  use Particles_sub
  use Messages

  implicit none

  include 'particles.h'

  real :: xp0=0.0, yp0=0.0, zp0=0.0, vpx0=0.0, vpy0=0.0, vpz0=0.0
  real :: delta_vp0=1.0, tausp=0.0, tausp1=0.0, eps_dtog=0.01
  real :: nu_epicycle=0.0, nu_epicycle2=0.0
  real :: beta_dPdr_dust=0.0, beta_dPdr_dust_scaled=0.0
  real :: taus1max=0.0, cdtp=0.2
  real :: gravx=0.0, gravz=0.0, kx_gg=1.0, kz_gg=1.0
  real :: Ri0=0.25, eps1=0.5
  real :: kx_xxp=0.0, ky_xxp=0.0, kz_xxp=0.0, amplxxp=0.0
  real :: kx_vvp=0.0, ky_vvp=0.0, kz_vvp=0.0, amplvvp=0.0
  real :: kx_vpx=0.0, kx_vpy=0.0, kx_vpz=0.0
  real :: ky_vpx=0.0, ky_vpy=0.0, ky_vpz=0.0
  real :: kz_vpx=0.0, kz_vpy=0.0, kz_vpz=0.0
  real :: phase_vpx=0.0, phase_vpy=0.0, phase_vpz=0.0
  real :: tstart_dragforce_par=0.0
  complex, dimension (7) :: coeff=(0.0,0.0)
  logical :: ldragforce_gas_par=.false., ldragforce_dust_par=.true.
  logical :: lpar_spec=.false.
  logical :: lsmooth_dragforce_dust=.false., lsmooth_dragforce_gas=.false.
  logical :: ldragforce_equi_global_eps=.false.
  logical, parameter :: ldraglaw_epstein=.true.
  logical :: lcoldstart_amplitude_correction=.false.
  logical :: linterpolate_spline=.true.
  character (len=labellen), dimension (ninit) :: initxxp='nothing'
  character (len=labellen), dimension (ninit) :: initvvp='nothing'
  character (len=labellen) :: gravx_profile='zero',  gravz_profile='zero'

  namelist /particles_init_pars/ &
      initxxp, initvvp, xp0, yp0, zp0, vpx0, vpy0, vpz0, delta_vp0, &
      bcpx, bcpy, bcpz, tausp, beta_dPdr_dust, rhop_tilde, &
      eps_dtog, nu_epicycle, lsmooth_dragforce_dust, &
      gravx_profile, gravz_profile, gravx, gravz, kx_gg, kz_gg, Ri0, eps1, &
      lmigration_redo, ldragforce_equi_global_eps, coeff, &
      kx_vvp, ky_vvp, kz_vvp, amplvvp, kx_xxp, ky_xxp, kz_xxp, amplxxp, &
      kx_vpx, kx_vpy, kx_vpz, ky_vpx, ky_vpy, ky_vpz, kz_vpx, kz_vpy, kz_vpz, &
      phase_vpx, phase_vpy, phase_vpz, lcoldstart_amplitude_correction, &
      lparticlemesh_cic, lparticlemesh_tsc, linterpolate_spline, &
      tstart_dragforce_par

  namelist /particles_run_pars/ &
      bcpx, bcpy, bcpz, tausp, dsnap_par_minor, beta_dPdr_dust, &
      ldragforce_gas_par, ldragforce_dust_par, lsmooth_dragforce_dust, &
      rhop_tilde, eps_dtog, cdtp, lpar_spec, &
      linterp_reality_check, nu_epicycle, &
      gravx_profile, gravz_profile, gravx, gravz, kx_gg, kz_gg, &
      lmigration_redo, tstart_dragforce_par, &
      lparticlemesh_cic, lparticlemesh_tsc

  integer :: idiag_xpm=0, idiag_ypm=0, idiag_zpm=0
  integer :: idiag_xp2m=0, idiag_yp2m=0, idiag_zp2m=0
  integer :: idiag_vpxm=0, idiag_vpym=0, idiag_vpzm=0
  integer :: idiag_vpx2m=0, idiag_vpy2m=0, idiag_vpz2m=0
  integer :: idiag_vpxmax=0, idiag_vpymax=0, idiag_vpzmax=0
  integer :: idiag_npm=0, idiag_np2m=0, idiag_npmax=0, idiag_npmin=0
  integer :: idiag_rhoptilm=0, idiag_dtdragp=0, idiag_nparmax=0
  integer :: idiag_rhopm=0, idiag_rhoprms=0, idiag_rhop2m=0, idiag_rhopmax=0
  integer :: idiag_rhopmin=0
  integer :: idiag_npmx=0, idiag_npmy=0, idiag_npmz=0
  integer :: idiag_rhopmx=0, idiag_rhopmy=0, idiag_rhopmz=0
  integer :: idiag_epspmx=0, idiag_epspmy=0, idiag_epspmz=0
  integer :: idiag_mpt=0, idiag_dedragp=0, idiag_rhopmxy=0

  contains

!***********************************************************************
    subroutine register_particles()
!
!  Set up indices for access to the fp and dfp arrays
!
!  29-dec-04/anders: coded
!
      use Mpicomm, only: stop_it
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_particles: called twice')
      first = .false.
!
      if (lroot) call cvs_id( &
           "$Id: particles_dust.f90,v 1.120 2006-07-19 20:46:16 ajohan Exp $")
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
!  Set indices for auxiliary variables
!
      inp   = mvar + naux + 1 + (maux_com - naux_com); naux = naux + 1
      irhop = mvar + naux + 1 + (maux_com - naux_com); naux = naux + 1
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
    endsubroutine register_particles
!***********************************************************************
    subroutine initialize_particles(lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  29-dec-04/anders: coded
!
      use EquationOfState, only: cs0
!
      logical :: lstarting
!
      real :: rhom
      integer, dimension (0:ncpus-1) :: ipar1, ipar2
!
!  Distribute particles evenly among processors to begin with.
!
      if (lstarting) call dist_particles_evenly_procs(npar_loc,ipar)
!
!  Size of box at local processor is needed for particle boundary conditions.
!
      Lxyz_loc(1)=Lxyz(1)/nprocx
      Lxyz_loc(2)=Lxyz(2)/nprocy
      Lxyz_loc(3)=Lxyz(3)/nprocz
      xyz0_loc(1)=xyz0(1)
      xyz0_loc(2)=xyz0(2)+ipy*Lxyz_loc(2)
      xyz0_loc(3)=xyz0(3)+ipz*Lxyz_loc(3)
      xyz1_loc(1)=xyz1(1)
      xyz1_loc(2)=xyz0(2)+(ipy+1)*Lxyz_loc(2)
      xyz1_loc(3)=xyz0(3)+(ipz+1)*Lxyz_loc(3)
!
!  The inverse stopping time is needed for drag force.
!
      tausp1=0.0
      if (tausp/=0.) tausp1=1/tausp
!
      if (beta_dPdr_dust/=0.0) then
        beta_dPdr_dust_scaled=beta_dPdr_dust*Omega/cs0
        if (lroot) print*, 'initialize_particles: Global pressure '// &
            'gradient with beta_dPdr_dust=', beta_dPdr_dust
      endif
!
!  Calculate mass density per particle (for back-reaction drag force on gas)
!  following the formula
!    rhop_tilde*N_cell = eps*rhom
!  where rhop_tilde is the mass density per particle, N_cell is the number of
!  particles per grid cell and rhom is the mean gas density in the box. 
!
      if (rhop_tilde==0.0) then
! For stratification, take into account gas present outside the simulation box.
        if (lgravz .and. lgravz_gas) then
          rhom=sqrt(2*pi)*1.0*1.0/Lz  ! rhom = Sigma/Lz, Sigma=sqrt(2*pi)*H*rho1
        else
          rhom=1.0
        endif
        rhop_tilde=eps_dtog*rhom/(real(npar)/(nxgrid*nygrid*nzgrid))
        if (lroot) then
          print*, 'initialize_particles: '// &
            'dust-to-gas ratio eps_dtog=', eps_dtog
          print*, 'initialize_particles: '// &
            'mass density per particle rhop_tilde=', rhop_tilde
        endif
      else
        if (lroot) print*, 'initialize_particles: '// &
            'mass density per particle rhop_tilde=', rhop_tilde
      endif
!
!  Calculate nu_epicycle**2 for gravity.
!
      if (nu_epicycle/=0.0) then
        gravz_profile='linear'
        nu_epicycle2=nu_epicycle**2
      endif
!
!  Gas density is needed for back-reaction friction force.
!
      if (ldragforce_gas_par .and. .not. ldensity) then
        if (lroot) then
          print*, 'initialize_particles: friction force on gas only works '
          print*, '                      together with gas density module!'
        endif
        call fatal_error('initialize_particles','')
      endif      
!
!  Need to map particles on the grid for dragforce on gas.
!      
      if (ldragforce_gas_par) then
        lcalc_np=.true.
!
!  When drag force is smoothed, df is also set in the first ghost zone. This 
!  region needs to be folded back into the df array after pde is finished,
!
        if (lparticlemesh_cic .or. lparticlemesh_tsc) lfold_df=.true.
      endif
!
!  Write constants to disc.
!      
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro')
          write (1,*) 'rhop_tilde=', rhop_tilde
        close (1)
      endif
!
    endsubroutine initialize_particles
!***********************************************************************
    subroutine init_particles(f,fp,ineargrid)
!
!  Initial positions and velocities of dust particles.
!
!  29-dec-04/anders: coded
!
      use EquationOfState, only: gamma, beta_glnrho_global, cs20
      use General, only: random_number_wrapper
      use Mpicomm, only: stop_it
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: uup
      real :: r, p, px, py, pz, eps, cs, k2_xxp
      real :: fac, npar_loc_x, npar_loc_y, npar_loc_z, dx_par, dy_par, dz_par
      integer :: l, j, k, ix0, iy0, iz0
      logical :: lequidistant=.false.
!
      intent (out) :: f, fp, ineargrid
!
!  Initial particle position.
!
      do j=1,ninit

        select case(initxxp(j))

        case ('nothing')
          if (lroot .and. j==1) print*, 'init_particles: nothing'

        case ('origin')
          if (lroot) print*, 'init_particles: All particles at origin'
          fp(1:npar_loc,ixp:izp)=0.
 
        case ('constant')
          if (lroot) &
              print*, 'init_particles: All particles at x,y,z=', xp0, yp0, zp0
          fp(1:npar_loc,ixp)=xp0
          fp(1:npar_loc,iyp)=yp0
          fp(1:npar_loc,izp)=zp0
 
        case ('random')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k=1,npar_loc
            if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
            if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
            if (nzgrid/=1) call random_number_wrapper(fp(k,izp))
          enddo
          if (nxgrid/=1) &
              fp(1:npar_loc,ixp)=xyz0_loc(1)+fp(1:npar_loc,ixp)*Lxyz_loc(1)
          if (nygrid/=1) &
              fp(1:npar_loc,iyp)=xyz0_loc(2)+fp(1:npar_loc,iyp)*Lxyz_loc(2)
          if (nzgrid/=1) &
              fp(1:npar_loc,izp)=xyz0_loc(3)+fp(1:npar_loc,izp)*Lxyz_loc(3)
 
        case ('np-constant')
          if (lroot) print*, 'init_particles: Constant number density'
          k=1
k_loop:   do while (.not. (k>npar_loc))
            do l=l1,l2; do m=m1,m2; do n=n1,n2
              if (nxgrid/=1) call random_number_wrapper(px)
              if (nygrid/=1) call random_number_wrapper(py)
              if (nzgrid/=1) call random_number_wrapper(pz)
              fp(k,ixp)=x(l)+(px-0.5)*dx
              fp(k,iyp)=y(m)+(py-0.5)*dy
              fp(k,izp)=z(n)+(pz-0.5)*dz
              k=k+1
              if (k>npar_loc) exit k_loop
            enddo; enddo; enddo
          enddo k_loop

        case ('equidistant')
          if (lroot) print*, 'init_particles: Particles placed equidistantly'
          fac=1.0/dimensionality
          dx_par=0.0; dy_par=0.0; dz_par=0.0
          fp(1,ixp)=x(l1); fp(1,iyp)=y(m1); fp(1,izp)=z(n1)
!  Number of particles in x-direction.          
          if (nxgrid/=1) then
            if (nygrid/=1) then
              npar_loc_x=(npar_loc*Lxyz_loc(1)/Lxyz_loc(2))**fac
            elseif (nzgrid/=1) then
              npar_loc_x=(npar_loc*Lxyz_loc(1)/Lxyz_loc(3))**fac
            else
              npar_loc_x=(npar_loc*Lxyz_loc(1)**2/(Lxyz_loc(2)*Lxyz_loc(3)))**fac
            endif
            dx_par=Lxyz_loc(1)/npar_loc_x
            fp(1,ixp) = xyz0_loc(1)+dx_par/2
          endif
!  Number of particles in y-direction.          
          if (nygrid/=1) then
            if (nxgrid/=1) then
              npar_loc_y=(npar_loc*Lxyz_loc(2)/Lxyz_loc(1))**fac
            elseif (nzgrid/=1) then
              npar_loc_y=(npar_loc*Lxyz_loc(2)/Lxyz_loc(3))**fac
            else
              npar_loc_y=(npar_loc*Lxyz_loc(2)**2/(Lxyz_loc(1)*Lxyz_loc(3)))**fac
            endif
            npar_loc_y=(npar_loc*Lxyz_loc(2)**2/(Lxyz_loc(1)*Lxyz_loc(3)))**fac
            dy_par=Lxyz_loc(2)/npar_loc_y
            fp(1,iyp) = xyz0_loc(2)+dy_par/2
          endif
!  Number of particles in z-direction.          
          if (nzgrid/=1) then
            if (nxgrid/=1) then
              npar_loc_z=(npar_loc*Lxyz_loc(3)/Lxyz_loc(1))**fac
            elseif (nygrid/=1) then
              npar_loc_z=(npar_loc*Lxyz_loc(3)/Lxyz_loc(2))**fac
            else
              npar_loc_z=(npar_loc*Lxyz_loc(3)**2/(Lxyz_loc(1)*Lxyz_loc(2)))**fac
            endif
            dz_par=Lxyz_loc(3)/npar_loc_z
            fp(1,izp) = xyz0_loc(3)+dz_par/2
          endif
!  Place particles iteratively, making sure that they are always in the box.
          do k=2,npar_loc
            fp(k,ixp)=fp(k-1,ixp)+dx_par
            fp(k,iyp)=fp(k-1,iyp)
            fp(k,izp)=fp(k-1,izp)
            if (fp(k,ixp)>xyz1_loc(1) .or. nxgrid==1) then
              fp(k,ixp)=fp(1,ixp)
              if (nygrid/=1) then
                fp(k,iyp)=fp(k,iyp)+dy_par
              else
                fp(k,izp)=fp(k,izp)+dz_par
              endif
            endif
            if (fp(k,iyp)>xyz1_loc(2)) then
              fp(k,iyp)=fp(1,iyp)
              if (nzgrid/=1) then
                fp(k,izp)=fp(k,izp)+dz_par
              endif
            endif
          enddo
          lequidistant=.true.
!
!  Shift particle locations slightly so that a mode appears.
!
        case ('shift')
          if (lroot) print*, 'init_particles: shift particle positions'
          if (.not. lequidistant) then
            if (lroot) print*, 'init_particles: must place particles equidistantly before shifting!'
            call fatal_error('init_particles','')
          endif
          k2_xxp=kx_xxp**2+ky_xxp**2+kz_xxp**2
          if (k2_xxp==0.0) then
            if (lroot) print*, &
                'init_particles: kx_xxp=ky_xxp=kz_xxp=0.0 is not allowed!'
            call fatal_error('init_particles','')
          endif
          do k=1,npar_loc
            fp(k,ixp) = fp(k,ixp) - kx_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp)+ky_xxp*fp(k,iyp)+kz_xxp*fp(k,izp))
            fp(k,iyp) = fp(k,iyp) - ky_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp)+ky_xxp*fp(k,iyp)+kz_xxp*fp(k,izp))
            fp(k,izp) = fp(k,izp) - kz_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp)+ky_xxp*fp(k,iyp)+kz_xxp*fp(k,izp))
          enddo
 
        case ('gaussian-z')
          if (lroot) print*, 'init_particles: Gaussian particle positions'
          do k=1,npar_loc
            if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
            if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
            call random_number_wrapper(r)
            call random_number_wrapper(p)
            if (nprocz==2) then
              if (ipz==0) fp(k,izp)=-abs(zp0*sqrt(-2*alog(r))*cos(2*pi*p))
              if (ipz==1) fp(k,izp)= abs(zp0*sqrt(-2*alog(r))*cos(2*pi*p))
            else
              fp(k,izp)= zp0*sqrt(-2*alog(r))*cos(2*pi*p)
            endif
          enddo
          if (nxgrid/=1) &
              fp(1:npar_loc,ixp)=xyz0_loc(1)+fp(1:npar_loc,ixp)*Lxyz_loc(1)
          if (nygrid/=1) &
              fp(1:npar_loc,iyp)=xyz0_loc(2)+fp(1:npar_loc,iyp)*Lxyz_loc(2)
 
        case ('streaming')
          call streaming(fp,f)
 
        case ('streaming_coldstart')
          call streaming_coldstart(fp,f)
 
        case ('constant-Ri')
          call constant_richardson(fp,f)
 
        case default
          if (lroot) &
              print*, 'init_particles: No such such value for initxxp: ', &
              trim(initxxp(j))
          call stop_it("")
 
        endselect

      enddo ! do j=1,ninit
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
!  Map particle position on the grid.
!
      call map_nearest_grid(f,fp,ineargrid)
      call map_xxp_grid(f,fp,ineargrid)
!
!  Initial particle velocity.
!
      do j=1,ninit

        select case(initvvp(j))

        case ('nothing')
          if (lroot.and.j==1) print*, 'init_particles: No particle velocity set'
        case ('zero')
          if (lroot) print*, 'init_particles: Zero particle velocity'
          fp(1:npar_loc,ivpx:ivpz)=0.
 
        case ('constant')
          if (lroot) print*, 'init_particles: Constant particle velocity'
          if (lroot) &
              print*, 'init_particles: vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
          fp(1:npar_loc,ivpx)=vpx0
          fp(1:npar_loc,ivpy)=vpy0
          fp(1:npar_loc,ivpz)=vpz0
 
        case ('sinwave-phase')
          if (lroot) print*, 'init_particles: sinwave-phase'
          if (lroot) &
              print*, 'init_particles: vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
          do k=1,npar_loc
            fp(k,ivpx)=fp(k,ivpx)+vpx0*sin(kx_vpx*fp(k,ixp)+ky_vpx*fp(k,iyp)+kz_vpx*fp(k,izp)+phase_vpx)
            fp(k,ivpy)=fp(k,ivpy)+vpy0*sin(kx_vpy*fp(k,ixp)+ky_vpy*fp(k,iyp)+kz_vpy*fp(k,izp)+phase_vpy)
            fp(k,ivpz)=fp(k,ivpz)+vpz0*sin(kx_vpz*fp(k,ixp)+ky_vpz*fp(k,iyp)+kz_vpz*fp(k,izp)+phase_vpz)
          enddo
 
        case ('coswave-phase')
          if (lroot) print*, 'init_particles: coswave-phase'
          if (lroot) &
              print*, 'init_particles: vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
          do k=1,npar_loc
            fp(k,ivpx)=fp(k,ivpx)+vpx0*cos(kx_vpx*fp(k,ixp)+ky_vpx*fp(k,iyp)+kz_vpx*fp(k,izp)+phase_vpx)
            fp(k,ivpy)=fp(k,ivpy)+vpy0*cos(kx_vpy*fp(k,ixp)+ky_vpy*fp(k,iyp)+kz_vpy*fp(k,izp)+phase_vpy)
            fp(k,ivpz)=fp(k,ivpz)+vpz0*cos(kx_vpz*fp(k,ixp)+ky_vpz*fp(k,iyp)+kz_vpz*fp(k,izp)+phase_vpz)
          enddo
 
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
            call interpolate_linear(f,iux,iuz,fp(k,ixp:izp),uup,ineargrid(k,:))
            fp(k,ivpx:ivpz) = uup
          enddo

        case('jeans-wave-dustpar-x')
        ! assumes rhs_poisson_const=1 !
          do k=1,npar_loc
            fp(k,ivpx) = fp(k,ivpx) - amplxxp* &
                (sqrt(1+4*1.0*1.0*tausp**2)-1)/ &
                (2*kx_xxp*1.0*tausp)*sin(kx_xxp*(fp(k,ixp)))
          enddo
 
        case('dragforce_equilibrium')
!
!  Equilibrium between drag forces on dust and gas and other forces
!  (from Nakagawa, Sekiya, & Hayashi 1986).
!
          if (lroot) then
            print*, 'init_particles: drag equilibrium'
            print*, 'init_particles: beta_glnrho_global=', beta_glnrho_global
          endif
!  Calculate average dust-to-gas ratio in box.
          if (ldensity_nolog) then
            eps = sum(f(l1:l2,m1:m2,n1:n2,irhop))/ &
                sum(f(l1:l2,m1:m2,n1:n2,ilnrho))
          else
            eps = sum(f(l1:l2,m1:m2,n1:n2,irhop))/ &
                sum(exp(f(l1:l2,m1:m2,n1:n2, ilnrho)))
          endif

          if (lroot) print*, 'init_particles: average dust-to-gas ratio=', eps
!  Set gas velocity field.
          do l=l1,l2; do m=m1,m2; do n=n1,n2
            cs=sqrt(cs20)
!  Take either global or local dust-to-gas ratio.
            if (.not. ldragforce_equi_global_eps) then
              if (ldensity_nolog) then
                eps = f(l,m,n,irhop)/f(l,m,n,ilnrho)
              else
                eps = f(l,m,n,irhop)/exp(f(l,m,n,ilnrho))
              endif
            endif
 
            f(l,m,n,iux) = f(l,m,n,iux) - &
                1/gamma*beta_glnrho_global(1)*eps*Omega*tausp/ &
                ((1.0+eps)**2+(Omega*tausp)**2)*cs
            f(l,m,n,iuy) = f(l,m,n,iuy) + &
                1/gamma*beta_glnrho_global(1)*(1+eps+(Omega*tausp)**2)/ &
                (2*((1.0+eps)**2+(Omega*tausp)**2))*cs
 
          enddo; enddo; enddo
!  Set particle velocity field.
          do k=1,npar_loc
!  Take either global or local dust-to-gas ratio.
            if (.not. ldragforce_equi_global_eps) then
              ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
              if (ldensity_nolog) then
                eps = f(ix0,iy0,iz0,irhop)/f(ix0,iy0,iz0,ilnrho)
              else
                eps = f(ix0,iy0,iz0,irhop)/exp(f(ix0,iy0,iz0,ilnrho))
              endif
            endif
            
            fp(k,ivpx) = fp(k,ivpx) + &
                1/gamma*beta_glnrho_global(1)*Omega*tausp/ &
                ((1.0+eps)**2+(Omega*tausp)**2)*cs
            fp(k,ivpy) = fp(k,ivpy) + &
                1/gamma*beta_glnrho_global(1)*(1+eps)/ &
                (2*((1.0+eps)**2+(Omega*tausp)**2))*cs
 
          enddo
 
        case('dragforce_equi_dust')
!
!  Equilibrium between drag force and Coriolis force on the dust.
!
          if (lroot) then
            print*, 'init_particles: drag equilibrium dust'
            print*, 'init_particles: beta_dPdr_dust=', beta_dPdr_dust
          endif
!  Set particle velocity field.
          cs=sqrt(cs20)
          do k=1,npar_loc
            fp(k,ivpx) = fp(k,ivpx) + &
                1/gamma*beta_dPdr_dust/ &
                (Omega*tausp+1/(Omega*tausp))*cs
            fp(k,ivpy) = fp(k,ivpy) - &
                1/gamma*beta_dPdr_dust*Omega*tausp*0.5/ &
                (Omega*tausp+1/(Omega*tausp))*cs
          enddo
 
        case default
          if (lroot) &
              print*, 'init_particles: No such such value for initvvp: ', &
              trim(initvvp(j))
          call stop_it("")

        endselect
!
      enddo ! do j=1,ninit
!
!  Sort particles (must happen at the end of the subroutine so that random
!  positions and velocities are not displaced relative to when there is no
!  sorting).
!
      call sort_particles_imn(fp,ineargrid,ipar)
!
    endsubroutine init_particles
!***********************************************************************
    subroutine streaming_coldstart(fp,f)
!
!  Mode that is unstable to the streaming instability of Youdin & Goodman (2005)
!
!  14-apr-06/anders: coded
!
      use EquationOfState, only: gamma, beta_glnrho_global
      use General, only: random_number_wrapper
!      
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      real :: eta_glnrho, v_Kepler, ampluug, dxp, dzp
      integer :: i, i1, i2, j, k, npar_loc_x, npar_loc_z
!
!  The number of particles per grid cell must be a quadratic number.
!
      if ( sqrt(npar/real(nwgrid))/=int(sqrt(npar/real(nwgrid))) .or. &
           sqrt(npar_loc/real(nw))/=int(sqrt(npar_loc/real(nw))) ) then
        if (lroot) then
          print*, 'streaming_coldstart: the number of particles per grid must'
          print*, '                     be a quadratic number!'
        endif
        print*, '                     iproc, npar/nw, npar_loc/nwgrid=', &
            iproc, npar/real(nwgrid), npar_loc/real(nw)
        call fatal_error('streaming_coldstart','')
      endif
!
!  Define a few disc parameters.
!
      eta_glnrho = -0.5*1/gamma*abs(beta_glnrho_global(1))*beta_glnrho_global(1)
      v_Kepler   =  1.0/abs(beta_glnrho_global(1))      
      if (lroot) print*, 'streaming: eta, vK=', eta_glnrho, v_Kepler
!
!  Place particles equidistantly.
!
      npar_loc_x=sqrt(npar_loc/(Lxyz_loc(3)/Lxyz_loc(1)))
      npar_loc_z=npar_loc/npar_loc_x
      dxp=Lxyz_loc(1)/npar_loc_x
      dzp=Lxyz_loc(3)/npar_loc_z
      do i=1,npar_loc_x
        i1=(i-1)*npar_loc_z+1; i2=i*npar_loc_z
        fp(i1:i2,ixp)=mod(i*dxp,Lxyz_loc(1))+dxp/2
        do j=i1,i2
          fp(j,izp)=xyz0_loc(3)+dzp/2+(j-i1)*dzp
        enddo
      enddo
!
!  Shift particle locations slightly so that wanted mode appears.
!
      do k=1,npar_loc
        fp(k,ixp) = fp(k,ixp) - &
            amplxxp/(2*(kx_xxp**2+kz_xxp**2))* &
            (kx_xxp*sin(kx_xxp*fp(k,ixp)+kz_xxp*fp(k,izp))+ &
             kx_xxp*sin(kx_xxp*fp(k,ixp)-kz_xxp*fp(k,izp)))
        fp(k,izp) = fp(k,izp) - &
            amplxxp/(2*(kx_xxp**2+kz_xxp**2))* &
            (kz_xxp*sin(kx_xxp*fp(k,ixp)+kz_xxp*fp(k,izp))- &
             kz_xxp*sin(kx_xxp*fp(k,ixp)-kz_xxp*fp(k,izp)))
        fp(k,ixp) = fp(k,ixp) + &
            kx_xxp/(2*(kx_xxp**2+kz_xxp**2))*amplxxp**2* &
            sin(2*(kx_xxp*fp(k,ixp)+kz_xxp*fp(k,izp)))
        fp(k,izp) = fp(k,izp) + &
            kz_xxp/(2*(kx_xxp**2+kz_xxp**2))*amplxxp**2* &
            sin(2*(kx_xxp*fp(k,ixp)+kz_xxp*fp(k,izp)))
      enddo
!  Set particle velocity.
      do k=1,npar_loc
        fp(k,ivpx) = fp(k,ivpx) + eta_glnrho*v_Kepler*amplxxp* &
            ( real(coeff(1))*cos(kx_xxp*fp(k,ixp)) - &
             aimag(coeff(1))*sin(kx_xxp*fp(k,ixp)))*cos(kz_xxp*fp(k,izp))
        fp(k,ivpy) = fp(k,ivpy) + eta_glnrho*v_Kepler*amplxxp* &
            ( real(coeff(2))*cos(kx_xxp*fp(k,ixp)) - &
             aimag(coeff(2))*sin(kx_xxp*fp(k,ixp)))*cos(kz_xxp*fp(k,izp))
        fp(k,ivpz) = fp(k,ivpz) + eta_glnrho*v_Kepler*(-amplxxp)* &
            (aimag(coeff(3))*cos(kx_xxp*fp(k,ixp)) + &
              real(coeff(3))*sin(kx_xxp*fp(k,ixp)))*sin(kz_xxp*fp(k,izp))
      enddo
!
!  Change the gas velocity amplitude so that the numerical error on the drag
!  force is corrected (the error is due to the interpolation of the gas
!  velocity field to the positions of the particles). A better way to correct
!  this is to go to a quadratic interpolation scheme.
!
      ampluug=amplxxp
      if (lcoldstart_amplitude_correction) &
          ampluug=amplxxp/(1-dx**2/8*(kx_xxp**2+kz_xxp**2))
!
!  Set fluid fields.
!
      do m=m1,m2; do n=n1,n2
        f(l1:l2,m,n,ilnrho) = f(l1:l2,m,n,ilnrho) + &
            amplxxp* &
            ( real(coeff(7))*cos(kx_xxp*x(l1:l2)) - &
             aimag(coeff(7))*sin(kx_xxp*x(l1:l2)))*cos(kz_xxp*z(n))
!
        f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) + &
            eta_glnrho*v_Kepler*ampluug* &
            ( real(coeff(4))*cos(kx_xxp*x(l1:l2)) - &
             aimag(coeff(4))*sin(kx_xxp*x(l1:l2)))*cos(kz_xxp*z(n))
!                
        f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + &
            eta_glnrho*v_Kepler*ampluug* &
            ( real(coeff(5))*cos(kx_xxp*x(l1:l2)) - &
             aimag(coeff(5))*sin(kx_xxp*x(l1:l2)))*cos(kz_xxp*z(n))
!
        f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) + &
            eta_glnrho*v_Kepler*(-ampluug)* &
            (aimag(coeff(6))*cos(kx_xxp*x(l1:l2)) + &
              real(coeff(6))*sin(kx_xxp*x(l1:l2)))*sin(kz_xxp*z(n))
      enddo; enddo
!
    endsubroutine streaming_coldstart
!***********************************************************************
    subroutine streaming(fp,f)
!
!  Mode that is unstable to the streaming instability of Youdin & Goodman (2005)
!
!  30-jan-06/anders: coded
!
      use EquationOfState, only: gamma, beta_glnrho_global
      use General, only: random_number_wrapper
!      
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      real :: eta_glnrho, v_Kepler, kx, kz
      real :: ampl, r, p, xprob, zprob, dxprob, dzprob, fprob, dfprob
      integer :: j, k
      logical :: lmigration_redo_org
!
!  Define a few disc parameters.
!
      eta_glnrho = -0.5*1/gamma*abs(beta_glnrho_global(1))*beta_glnrho_global(1)
      v_Kepler   =  1.0/abs(beta_glnrho_global(1))      
      if (lroot) print*, 'streaming: eta, vK=', eta_glnrho, v_Kepler
!
!  Place particles according to probability function.
!
!  Invert
!    r = x
!    p = int_0^z f(x,z') dz' = z + A/kz*cos(kx*x)*sin(kz*z)
!  where r and p are random numbers between 0 and 1.
      kx=kx_xxp*Lxyz(1); kz=kz_xxp*Lxyz(3)
      do k=1,npar_loc

        call random_number_wrapper(r)
        call random_number_wrapper(p)

        fprob = 1.0
        zprob = 0.0

        j=0
!  Use Newton-Raphson iteration to invert function.
        do while ( abs(fprob)>0.0001 )

          xprob = r
          fprob = zprob + amplxxp/kz*cos(kx*xprob)*sin(kz*zprob) - p
          dfprob= 1.0 + amplxxp*cos(kx*xprob)*cos(kz*zprob) 
          dzprob= -fprob/dfprob
          zprob = zprob+0.2*dzprob

          j=j+1

        enddo

        if ( mod(k,npar_loc/100)==0) then
          print '(i7,i3,4f11.7)', k, j, r, p, xprob, zprob
        endif

        fp(k,ixp)=xprob*Lxyz(1)+xyz0(1)
        fp(k,izp)=zprob*Lxyz(3)+xyz0(3)
!  Set particle velocity.
        fp(k,ivpx) = fp(k,ivpx) + eta_glnrho*v_Kepler*amplxxp* &
            ( real(coeff(1))*cos(kx_xxp*fp(k,ixp)) - &
             aimag(coeff(1))*sin(kx_xxp*fp(k,ixp)))*cos(kz_xxp*fp(k,izp))
        fp(k,ivpy) = fp(k,ivpy) + eta_glnrho*v_Kepler*amplxxp* &
            ( real(coeff(2))*cos(kx_xxp*fp(k,ixp)) - &
             aimag(coeff(2))*sin(kx_xxp*fp(k,ixp)))*cos(kz_xxp*fp(k,izp))
        fp(k,ivpz) = fp(k,ivpz) + eta_glnrho*v_Kepler*(-amplxxp)* &
            (aimag(coeff(3))*cos(kx_xxp*fp(k,ixp)) + &
              real(coeff(3))*sin(kx_xxp*fp(k,ixp)))*sin(kz_xxp*fp(k,izp))

      enddo
!
!  Particles were placed randomly in the entire simulation space, so they need
!  to be send to the correct processors now.
!
      if (lmpicomm) then
        lmigration_redo_org=lmigration_redo
        lmigration_redo=.true.
        call redist_particles_procs(fp,npar_loc,ipar)
        lmigration_redo=lmigration_redo_org
      endif
!
!  Set fluid fields.
!
      do m=m1,m2; do n=n1,n2
        f(l1:l2,m,n,ilnrho) = f(l1:l2,m,n,ilnrho) + &
            (eta_glnrho*v_Kepler)**2*amplxxp* &
            ( real(coeff(7))*cos(kx_xxp*x(l1:l2)) - &
             aimag(coeff(7))*sin(kx_xxp*x(l1:l2)))*cos(kz_xxp*z(n))
!                
        f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) + &
            eta_glnrho*v_Kepler*amplxxp* &
            ( real(coeff(4))*cos(kx_xxp*x(l1:l2)) - &
             aimag(coeff(4))*sin(kx_xxp*x(l1:l2)))*cos(kz_xxp*z(n))
!                
        f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + &
            eta_glnrho*v_Kepler*amplxxp* &
            ( real(coeff(5))*cos(kx_xxp*x(l1:l2)) - &
             aimag(coeff(5))*sin(kx_xxp*x(l1:l2)))*cos(kz_xxp*z(n))
!
        f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) + &
            eta_glnrho*v_Kepler*(-amplxxp)* &
            (aimag(coeff(6))*cos(kx_xxp*x(l1:l2)) + &
              real(coeff(6))*sin(kx_xxp*x(l1:l2)))*sin(kz_xxp*z(n))
      enddo; enddo
!
    endsubroutine streaming
!***********************************************************************
    subroutine constant_richardson(fp,f)
!
!  Setup dust density with a constant Richardson number (Sekiya, 1998).
!    eps=1/sqrt(z^2/Hd^2+1/(1+eps1)^2)-1
!
!  14-sep-05/anders: coded
!
      use EquationOfState, only: beta_glnrho_scaled, gamma, cs20
      use General, only: random_number_wrapper
!      
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      integer, parameter :: nz_inc=10
      real, dimension (nz_inc*nz) :: z_dense, eps
      real, dimension (nx) :: np
      real :: r, p, Hg, Hd, frac, rho1, Sigmad, Sigmad_num, Xi, fXi, dfdXi
      real :: dz_dense, eps_point, z00_dense, rho, lnrho
      integer :: nz_dense=nz_inc*nz, npar_bin
      integer :: i, i0, k
!
!  Calculate dust "scale height".
!
      rho1=1.0
      Hg=1.0
      Sigmad=eps_dtog*rho1*Hg*sqrt(2*pi)
      Hd = sqrt(Ri0)*abs(beta_glnrho_scaled(1))/(2*gamma)*1.0
!
!  Need to find eps1 that results in given dust column density.
!
      Xi = sqrt(eps1*(2+eps1))/(1+eps1)
      fXi=-2*Xi + alog((1+Xi)/(1-Xi))-Sigmad/(Hd*rho1)
      i=0
!
!  Newton-Raphson on equation Sigmad/(Hd*rho1)=-2*Xi + alog((1+Xi)/(1-Xi)).
!  Here Xi = sqrt(eps1*(2+eps1))/(1+eps1).
!
      do while (abs(fXi)>=0.00001)
        
        dfdXi=2*Xi**2/(1-Xi**2)
        Xi=Xi-0.1*fXi/dfdXi
         
        fXi=-2*Xi + alog((1+Xi)/(1-Xi))-Sigmad/(Hd*rho1)
             
        i=i+1
        if (i>=1000) stop
                 
      enddo
!
!  Calculate eps1 from Xi.
!      
      eps1=-1+1/sqrt(-(Xi**2)+1)
      if (lroot) print*, 'constant_richardson: Hd, eps1=', Hd, eps1
!
!  Make z denser for higher resolution in density.
!
      dz_dense=Lxyz_loc(3)/nz_dense
      z00_dense=xyz0_loc(3)+0.5*dz_dense
      do n=1,nz_dense
        z_dense(n)=z00_dense+(n-1)*dz_dense
      enddo
!
!  Dust-to-gas ratio as a function of z (with cutoff).
!
      eps=1/sqrt(z_dense**2/Hd**2+1/(1+eps1)**2)-1
      where (eps<=0.0) eps=0.0
!
!  Calculate the dust column density numerically.
!
      Sigmad_num=sum(rho1*eps*dz_dense)
      if (lroot) print*, 'constant_richardson: Sigmad, Sigmad (numerical) = ', &
          Sigmad, Sigmad_num
!
!  Place particles according to probability function.
!
      i0=0
      do n=1,nz_dense
        frac=eps(n)/Sigmad_num*dz_dense
        npar_bin=int(frac*npar_loc)
        if (npar_bin>=2.and.mod(n,2)==0) npar_bin=npar_bin+1
        do i=i0+1,i0+npar_bin
          if (i<=npar_loc) then
            call random_number_wrapper(r)
            fp(i,izp)=z_dense(n)+(2*r-1.0)*dz_dense/2
          endif
        enddo
        i0=i0+npar_bin
      enddo
      if (lroot) print '(A,i7,A)', 'constant_richardson: placed ', &
          i0, ' particles according to Ri=const'
!
!  Particles left out by round off are just placed randomly.
!      
      if (i0+1<=npar_loc) then
        do k=i0+1,npar_loc
          call random_number_wrapper(fp(k,izp))
          fp(k,izp)=xyz0(3)+fp(k,izp)*Lxyz(3)
        enddo
        if (lroot) print '(A,i7,A)', 'constant_richardson: placed ', &
            npar_loc-i0, ' particles randomly.'
      endif
!
!  Random positions in x and y.
!      
      do k=1,npar_loc
        if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
        if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
      enddo
      if (nxgrid/=1) &
          fp(1:npar_loc,ixp)=xyz0_loc(1)+fp(1:npar_loc,ixp)*Lxyz_loc(1)
      if (nygrid/=1) &
          fp(1:npar_loc,iyp)=xyz0_loc(2)+fp(1:npar_loc,iyp)*Lxyz_loc(2)
!       
!  Set gas velocity according to dust-to-gas ratio and global pressure gradient.
!          
      do imn=1,ny*nz

        n=nn(imn); m=mm(imn)

        if (abs(z(n))<=Hd*sqrt(1-1/(1+eps1)**2)) then
          lnrho = -sqrt(z(n)**2/Hd**2+1/(1+eps1)**2)* &
              gamma*Omega**2*Hd**2/cs20 + gamma*Omega**2*Hd**2/(cs20*(1+eps1))
        else
          lnrho = -0.5*gamma*Omega**2/cs20*z(n)**2 + &
              gamma*Omega**2*Hd**2/cs20*(1/(1+eps1)-1/(2*(1+eps1)**2) - 0.5)
        endif
!
!  Isothermal stratification.
!        
        if (lentropy) f(l1:l2,m,n,iss) = (1/gamma-1.0)*lnrho

        rho=exp(lnrho)

        if (ldensity_nolog) then
          f(l1:l2,m,n,ilnrho)=rho
        else
          f(l1:l2,m,n,ilnrho)=lnrho
        endif

        eps_point=1/sqrt(z(n)**2/Hd**2+1/(1+eps1)**2)-1
        if (eps_point<=0.0) eps_point=0.0

        f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) - &
            1/gamma*cs20*beta_glnrho_scaled(1)*eps_point*tausp/ &
            (1.0+2*eps_point+eps_point**2+(Omega*tausp)**2)
        f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + &
            1/gamma*cs20*beta_glnrho_scaled(1)*(1+eps_point+(Omega*tausp)**2)/ &
            (2*Omega*(1.0+2*eps_point+eps_point**2+(Omega*tausp)**2))
        f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) + 0.0
      enddo
!
!  Set particle velocity.
!      
      do k=1,npar_loc

        eps_point=1/sqrt(fp(k,izp)**2/Hd**2+1/(1+eps1)**2)-1
        if (eps_point<=0.0) eps_point=0.0

        fp(k,ivpx) = fp(k,ivpx) + &
            1/gamma*cs20*beta_glnrho_scaled(1)*tausp/ &
            (1.0+2*eps_point+eps_point**2+(Omega*tausp)**2)
        fp(k,ivpy) = fp(k,ivpy) + &
            1/gamma*cs20*beta_glnrho_scaled(1)*(1+eps_point)/ &
            (2*Omega*(1.0+2*eps_point+eps_point**2+(Omega*tausp)**2))
        fp(k,ivpz) = fp(k,ivpz) - tausp*Omega**2*fp(k,izp)

      enddo
!
    endsubroutine constant_richardson
!***********************************************************************
    subroutine pencil_criteria_particles()
! 
!  All pencils that the Particles module depends on are specified here.
! 
!  20-04-06/anders: coded
!
      use Cdata
!
      if (ldragforce_gas_par) lpenc_requested(i_rho1)=.true.
!
      lpenc_diagnos(i_np)=.true.
      lpenc_diagnos(i_rhop)=.true.
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p
!
      if (lpencil(i_rhop)) then
        if (irhop/=0) then
          p%rhop=f(l1:l2,m,n,irhop)
        else
          p%rhop=rhop_tilde*f(l1:l2,m,n,inp)
        endif
      endif
!
      if (lpencil(i_np)) p%np=f(l1:l2,m,n,inp)
!
    endsubroutine calc_pencils_particles
!***********************************************************************
    subroutine dxxp_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle position (called from main pencil loop).
!
!  25-apr-06/anders: dummy
!
      real, dimension (mx,my,mz,mvar+maux) :: f
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
!  Evolution of dust particle velocity (called from main pencil loop).
!
!  25-apr-06/anders: coded
!
      use Cdata
      use EquationOfState, only: cs20, gamma
      use Mpicomm, only: stop_it
      use Particles_number, only: get_nptilde
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (nx) :: tausg1, dt1_drag
      real, dimension (3) :: uup, dragforce
      real :: np_point, eps_point, rho_point, rho1_point, tausp1_point, up2
      real :: weight, weight_x, weight_y, weight_z
      real :: dt1_advpx, dt1_advpy, dt1_advpz
      integer :: k, l, ix0, iy0, iz0
      integer :: ixx, iyy, izz, ixx0, iyy0, izz0, ixx1, iyy1, izz1
!
      intent (in) :: f, fp, ineargrid
      intent (inout) :: df, dfp
! 
!  Add drag force if stopping time is not infinite.
! 
      if (ldragforce_dust_par .and. t>=tstart_dragforce_par) then
        if (headtt) print*,'dvvp_dt: Add drag force; tausp=', tausp
        if (npar_imn(imn)/=0) then
          do k=k1_imn(imn),k2_imn(imn)
            ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
            call get_frictiontime(f,fp,ineargrid,k,tausp1_point)
!  Use interpolation to calculate gas velocity at position of particles.
            if (lhydro) then
              if (lparticlemesh_cic) then
                call interpolate_linear( &
                    f,iux,iuz,fp(k,ixp:izp),uup,ineargrid(k,:),ipar(k) )
              elseif (lparticlemesh_tsc) then
                if (linterpolate_spline) then
                  call interpolate_quadratic_spline( &
                      f,iux,iuz,fp(k,ixp:izp),uup,ineargrid(k,:),ipar(k) )
                else
                  call interpolate_quadratic( &
                      f,iux,iuz,fp(k,ixp:izp),uup,ineargrid(k,:),ipar(k) )
                endif
              else
                uup=f(ix0,iy0,iz0,iux:iuz)
              endif
            else
              uup=0.0
            endif
            dragforce = -tausp1_point*(fp(k,ivpx:ivpz)-uup)
            dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + dragforce
!
!  Drag force diagnostics
!
            if (ldiagnos) then
              if (idiag_dedragp/=0) then
                if (ldragforce_gas_par) then
                  up2=sum((fp(k,ivpx:ivpz)-uup)**2)
                  call sum_par_name((/-rhop_tilde*tausp1_point*up2/),idiag_dedragp)
                else
                  up2=sum(fp(k,ivpx:ivpz)*(fp(k,ivpx:ivpz)-uup))
                  call sum_par_name((/-rhop_tilde*tausp1_point*up2/),idiag_dedragp)
                endif
              endif
            endif
!            
!  Back-reaction friction force from particles on gas. Three methods are
!  implemented for assigning a particle to the mesh (see Hockney & Eastwood):
!
!    0. NGP (Nearest Grid Point)
!       The entire effect of the particle goes to the nearest grid point. 
!    1. CIC (Cloud In Cell)
!       The particle has a region of influence with the size of a grid cell.
!       This is equivalent to a first order (spline) interpolation scheme.
!    2. TSC (Triangular Shaped Cloud)
!       The particle is spread over a length of two grid cells, but with
!       a density that falls linearly outwards.
!       This is equivalent to a second order spline interpolation scheme.
!
            if (ldragforce_gas_par) then
!
!  Cloud In Cell (CIC) scheme.
!              
              if (lparticlemesh_cic) then
                ixx0=ix0; iyy0=iy0; izz0=iz0
                ixx1=ix0; iyy1=iy0; izz1=iz0
!
!  Particle influences the 8 surrounding grid points. The reference point is
!  the grid point at the lower left corner.
!
                if ( (x(ix0)>fp(k,ixp)) .and. nxgrid/=1) ixx0=ixx0-1
                if ( (y(iy0)>fp(k,iyp)) .and. nygrid/=1) iyy0=iyy0-1
                if ( (z(iz0)>fp(k,izp)) .and. nzgrid/=1) izz0=izz0-1
                if (nxgrid/=1) ixx1=ixx0+1
                if (nygrid/=1) iyy1=iyy0+1
                if (nzgrid/=1) izz1=izz0+1
                do izz=izz0,izz1; do iyy=iyy0,iyy1; do ixx=ixx0,ixx1
                  weight=1.0
                  if (nxgrid/=1) &
                      weight=weight*( 1.0-abs(fp(k,ixp)-x(ixx))*dx_1(ixx) )
                  if (nygrid/=1) &
                      weight=weight*( 1.0-abs(fp(k,iyp)-y(iyy))*dy_1(iyy) )
                  if (nzgrid/=1) &
                      weight=weight*( 1.0-abs(fp(k,izp)-z(izz))*dz_1(izz) )
!  Save the calculation of rho1 when inside pencil.
                  if ( (iyy/=m) .or. (izz/=n) .or. (ixx<l1) .or. (ixx>l2) ) then
                    rho_point=f(ixx,iyy,izz,ilnrho)
                    if (.not. ldensity_nolog) rho_point=exp(rho_point)
                    rho1_point=1/rho_point
                  else
                    rho1_point=p%rho1(ixx-nghost)
                  endif
!  Add friction force to grid point.                  
                  df(ixx,iyy,izz,iux:iuz)=df(ixx,iyy,izz,iux:iuz) - &
                      rhop_tilde*rho1_point*dragforce*weight
                enddo; enddo; enddo
!
!  Triangular Shaped Cloud (TSC) scheme.
!                
              elseif (lparticlemesh_tsc) then
!
!  Particle influences the 27 surrounding grid points, but has a density that
!  decreases with the distance from the particle centre.
!
                if (nxgrid/=1) then
                  ixx0=ix0-1; ixx1=ix0+1
                else
                  ixx0=ix0  ; ixx1=ix0
                endif
                if (nygrid/=1) then
                  iyy0=iy0-1; iyy1=iy0+1
                else
                  iyy0=iy0  ; iyy1=iy0
                endif
                if (nzgrid/=1) then
                  izz0=iz0-1; izz1=iz0+1
                else
                  izz0=iz0  ; izz1=iz0
                endif
!
!  The nearest grid point is influenced differently than the left and right
!  neighbours are. A particle that is situated exactly on a grid point gives
!  3/4 contribution to that grid point and 1/8 to each of the neighbours.
!
                do izz=izz0,izz1; do iyy=iyy0,iyy1; do ixx=ixx0,ixx1
                  if ( ((ixx-ix0)==-1) .or. ((ixx-ix0)==+1) ) then
                    weight_x = 1.125 - 1.5* abs(fp(k,ixp)-x(ixx))*dx_1(ixx) + &
                                       0.5*(abs(fp(k,ixp)-x(ixx))*dx_1(ixx))**2
                  else
                    if (nxgrid/=1) &
                    weight_x = 0.75  -       ((fp(k,ixp)-x(ixx))*dx_1(ixx))**2
                  endif
                  if ( ((iyy-iy0)==-1) .or. ((iyy-iy0)==+1) ) then
                    weight_y = 1.125 - 1.5* abs(fp(k,iyp)-y(iyy))*dy_1(iyy) + &
                                       0.5*(abs(fp(k,iyp)-y(iyy))*dy_1(iyy))**2
                  else
                    if (nygrid/=1) &
                    weight_y = 0.75  -       ((fp(k,iyp)-y(iyy))*dy_1(iyy))**2
                  endif
                  if ( ((izz-iz0)==-1) .or. ((izz-iz0)==+1) ) then
                    weight_z = 1.125 - 1.5* abs(fp(k,izp)-z(izz))*dz_1(izz) + &
                                       0.5*(abs(fp(k,izp)-z(izz))*dz_1(izz))**2
                  else
                    if (nzgrid/=1) &
                    weight_z = 0.75  -       ((fp(k,izp)-z(izz))*dz_1(izz))**2
                  endif
                  
                  weight=1.0

                  if (nxgrid/=1) weight=weight*weight_x
                  if (nygrid/=1) weight=weight*weight_y
                  if (nzgrid/=1) weight=weight*weight_z
!  Save the calculation of rho1 when inside pencil.
                  if ( (iyy/=m) .or. (izz/=n) .or. (ixx<l1) .or. (ixx>l2) ) then
                    rho_point=f(ixx,iyy,izz,ilnrho)
                    if (.not. ldensity_nolog) rho_point=exp(rho_point)
                    rho1_point=1/rho_point
                  else
                    rho1_point=p%rho1(ixx-nghost)
                  endif
!  Add friction force to grid point.
                  df(ixx,iyy,izz,iux:iuz)=df(ixx,iyy,izz,iux:iuz) - &
                      rhop_tilde*rho1_point*dragforce*weight
                enddo; enddo; enddo
              else
!
!  Nearest Grid Point (NGP) scheme.
!                
                l=ineargrid(k,1)
                df(l,m,n,iux:iuz) = df(l,m,n,iux:iuz) - &
                    rhop_tilde*p%rho1(l-nghost)*dragforce
              endif
            endif
          enddo
!          
!  Contribution of friction force to time-step.
!
          if (lfirst.and.ldt) then
            if (ldragforce_gas_par) then
              tausg1  =f(l1:l2,m,n,irhop)*p%rho1*tausp1_point
              dt1_drag=(tausp1_point+tausg1)/cdtp
            else
              dt1_drag=tausp1_point/cdtp
            endif
            dt1_max=max(dt1_max,dt1_drag)
            if (ldiagnos.and.idiag_dtdragp/=0) &
                call max_mn_name(dt1_drag,idiag_dtdragp,l_dt=.true.)
          endif
        endif
      endif
!
!  Contribution of dust particles to drag force.
!
      if (npar_imn(imn)/=0) then
        do k=k1_imn(imn),k2_imn(imn)
          if (lfirst.and.ldt) then
            ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
            dt1_advpx=fp(k,ivpx)*dx_1(ix0)/cdtp
            dt1_advpy=fp(k,ivpy)*dy_1(iy0)/cdtp
            dt1_advpz=fp(k,ivpz)*dz_1(iz0)/cdtp

            dt1_max(ix0-nghost)=max(dt1_max(ix0-nghost),dt1_advpx)
            dt1_max(ix0-nghost)=max(dt1_max(ix0-nghost),dt1_advpy)
            dt1_max(ix0-nghost)=max(dt1_max(ix0-nghost),dt1_advpz)
          endif 
        enddo
      endif
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_npm/=0)     call sum_mn_name(p%np,idiag_npm)
        if (idiag_np2m/=0)    call sum_mn_name(p%np**2,idiag_np2m)
        if (idiag_npmax/=0)   call max_mn_name(p%np,idiag_npmax)
        if (idiag_npmin/=0)   call max_mn_name(-p%np,idiag_npmin,lneg=.true.)
        if (idiag_rhopm/=0)   call sum_mn_name(p%rhop,idiag_rhopm)
        if (idiag_rhop2m/=0)  call sum_mn_name(p%rhop**2,idiag_rhop2m)
        if (idiag_rhoprms/=0) call sum_mn_name(p%rhop**2,idiag_rhoprms,lsqrt=.true.)
        if (idiag_rhopmax/=0) call max_mn_name(p%rhop,idiag_rhopmax)
        if (idiag_rhopmin/=0) call max_mn_name(-p%rhop,idiag_rhopmin,lneg=.true.)
        if (idiag_npmx/=0)    call yzsum_mn_name_x(p%np,idiag_npmx)
        if (idiag_npmy/=0)    call xzsum_mn_name_y(p%np,idiag_npmy)
        if (idiag_npmz/=0)    call xysum_mn_name_z(p%np,idiag_npmz)
        if (idiag_rhopmx/=0)  call yzsum_mn_name_x(p%rhop,idiag_rhopmx)
        if (idiag_rhopmy/=0)  call xzsum_mn_name_y(p%rhop,idiag_rhopmy)
        if (idiag_rhopmz/=0)  call xysum_mn_name_z(p%rhop,idiag_rhopmz)
        if (idiag_epspmx/=0)  call yzsum_mn_name_x(p%rhop*p%rho1,idiag_epspmx)
        if (idiag_epspmy/=0)  call xzsum_mn_name_y(p%rhop*p%rho1,idiag_epspmy)
        if (idiag_epspmz/=0)  call xysum_mn_name_z(p%rhop*p%rho1,idiag_epspmz)
        if (idiag_rhopmxy/=0) call zsum_mn_name_xy(p%rhop,idiag_rhopmxy)
      endif
!
    endsubroutine dvvp_dt_pencil
!***********************************************************************
    subroutine dxxp_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of dust particle position.
!
!  02-jan-05/anders: coded
!
      use General, only: random_number_wrapper, random_seed_wrapper
!      
      real, dimension (mx,my,mz,mvar+maux) :: f
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
      if (lshear.and.nygrid/=1) dfp(1:npar_loc,iyp) = &
          dfp(1:npar_loc,iyp) - qshear*Omega*fp(1:npar_loc,ixp)
!
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine dxxp_dt
!***********************************************************************
    subroutine dvvp_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of dust particle velocity.
!
!  29-dec-04/anders: coded
!
      use Cdata
      use EquationOfState, only: cs20, gamma
      use Mpicomm, only: stop_it
      use Particles_number, only: get_nptilde
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real :: Omega2, np_tilde
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
!  Add constant background pressure gradient beta=alpha*H0/r0, where alpha
!  comes from a global pressure gradient P = P0*(r/r0)^alpha.
!  (the term must be added to the dust equation of motion when measuring
!  velocities relative to the shear flow modified by the global pressure grad.)
!
      if (beta_dPdr_dust/=0.0 .and. t>=tstart_dragforce_par) then
        dfp(1:npar_loc,ivpx) = &
            dfp(1:npar_loc,ivpx) + 1/gamma*cs20*beta_dPdr_dust_scaled
      endif
!
!  Gravity on the particles.
!
      select case (gravx_profile)

        case ('zero')
          if (lheader) print*, 'dvvp_dt: No gravity in x-direction.'
 
        case ('sinusoidal')
          if (lheader) &
              print*, 'dvvp_dt: Sinusoidal gravity field in x-direction.'
          dfp(1:npar_loc,ivpx)=dfp(1:npar_loc,ivpx) - &
              gravx*sin(kx_gg*fp(1:npar_loc,ixp))
 
        case ('default')
          call fatal_error('dvvp_dt','chosen gravx_profile is not valid!')

      endselect
!
      select case (gravz_profile)

        case ('zero')
          if (lheader) print*, 'dvvp_dt: No gravity in z-direction.'
 
        case ('linear')
          if (lheader) print*, 'dvvp_dt: Linear gravity field in z-direction.'
          dfp(1:npar_loc,ivpz)=dfp(1:npar_loc,ivpz) - &
              nu_epicycle2*fp(1:npar_loc,izp)
 
        case ('sinusoidal')
          if (lheader) &
              print*, 'dvvp_dt: Sinusoidal gravity field in z-direction.'
          dfp(1:npar_loc,ivpz)=dfp(1:npar_loc,ivpz) - &
              gravz*sin(kz_gg*fp(1:npar_loc,izp))
 
        case ('default')
          call fatal_error('dvvp_dt','chosen gravz_profile is not valid!')

      endselect
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_nparmax/=0) call max_name(npar_loc,idiag_nparmax)
        if (idiag_xpm/=0)  call sum_par_name(fp(1:npar_loc,ixp),idiag_xpm)
        if (idiag_ypm/=0)  call sum_par_name(fp(1:npar_loc,iyp),idiag_ypm)
        if (idiag_zpm/=0)  call sum_par_name(fp(1:npar_loc,izp),idiag_zpm)
        if (idiag_xp2m/=0) call sum_par_name(fp(1:npar_loc,ixp)**2,idiag_xp2m)
        if (idiag_yp2m/=0) call sum_par_name(fp(1:npar_loc,iyp)**2,idiag_yp2m)
        if (idiag_zp2m/=0) call sum_par_name(fp(1:npar_loc,izp)**2,idiag_zp2m)
        if (idiag_vpxm/=0) call sum_par_name(fp(1:npar_loc,ivpx),idiag_vpxm)
        if (idiag_vpym/=0) call sum_par_name(fp(1:npar_loc,ivpy),idiag_vpym)
        if (idiag_vpzm/=0) call sum_par_name(fp(1:npar_loc,ivpz),idiag_vpzm)
        if (idiag_vpx2m/=0) &
            call sum_par_name(fp(1:npar_loc,ivpx)**2,idiag_vpx2m)
        if (idiag_vpy2m/=0) &
            call sum_par_name(fp(1:npar_loc,ivpy)**2,idiag_vpy2m)
        if (idiag_vpz2m/=0) &
            call sum_par_name(fp(1:npar_loc,ivpz)**2,idiag_vpz2m)
        if (idiag_vpxmax/=0) call max_par_name(fp(1:npar_loc,ivpx),idiag_vpxmax)
        if (idiag_vpymax/=0) call max_par_name(fp(1:npar_loc,ivpy),idiag_vpymax)
        if (idiag_vpzmax/=0) call max_par_name(fp(1:npar_loc,ivpz),idiag_vpzmax)
        if (idiag_rhoptilm/=0) then
          do k=1,npar_loc
            call get_nptilde(fp,k,np_tilde)
            call sum_par_name( &
                (/4/3.*pi*rhops*fp(k,iap)**3*np_tilde/),idiag_rhoptilm)
          enddo
        endif
        if (idiag_mpt/=0) then
          do k=1,npar_loc
            call get_nptilde(fp,k,np_tilde)
            call integrate_par_name( &
                (/4/3.*pi*rhops*fp(k,iap)**3*np_tilde/),idiag_mpt)
          enddo
        endif
      endif
!
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine dvvp_dt
!***********************************************************************
    subroutine get_frictiontime(f,fp,ineargrid,k,tausp1_point)
!
!  Calculate the friction time.
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      real :: tausp1_point
      integer, dimension (mpar_loc,3) :: ineargrid
      integer :: k
!
      if (ldraglaw_epstein) then
        if (iap/=0) then
          tausp1_point=1/(fp(k,iap)*rhops)
        else
          tausp1_point=tausp1
        endif
      endif
!
    endsubroutine get_frictiontime
!***********************************************************************
    subroutine get_distances(f,fp,rp_mn,rpcyl_mn)
!
! 18-jul-06/wlad: dummy routine
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension (nx,mpar_loc) :: rp_mn,rpcyl_mn
      integer :: i
!
      intent(out) :: rp_mn,rpcyl_mn
       do i=1,mpar_loc
          rp_mn(:,i)    = 0.
          rpcyl_mn(:,i) = 0.
       enddo
!
    endsubroutine get_distances
!***********************************************************************
    subroutine read_particles_init_pars(unit,iostat)
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
      integer, intent (in) :: unit
!
      write(unit,NML=particles_init_pars)
!
    endsubroutine write_particles_init_pars
!***********************************************************************
    subroutine read_particles_run_pars(unit,iostat)
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
!  01-jan-06/anders: coded
!
      use Power_spectrum, only: power_1d
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if (lpar_spec) call power_1d(f,'p',0,irhop)
!
    endsubroutine powersnap_particles
!***********************************************************************
    subroutine rprint_particles(lreset,lwrite)
!   
!  Read and register print parameters relevant for particles
!
!  29-dec-04/anders: coded
!
      use Cdata
      use Sub, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname,inamez,inamey,inamex, inamexy
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
        write(3,*) 'inp=', inp
        write(3,*) 'irhop=', irhop
      endif
!
!  Reset everything in case of reset
!
      if (lreset) then
        idiag_xpm=0; idiag_ypm=0; idiag_zpm=0
        idiag_xp2m=0; idiag_yp2m=0; idiag_zp2m=0
        idiag_vpxm=0; idiag_vpym=0; idiag_vpzm=0
        idiag_vpx2m=0; idiag_vpy2m=0; idiag_vpz2m=0
        idiag_vpxmax=0; idiag_vpymax=0; idiag_vpzmax=0
        idiag_npm=0; idiag_np2m=0; idiag_npmax=0; idiag_npmin=0
        idiag_rhoptilm=0; idiag_dtdragp=0; idiag_dedragp=0
        idiag_rhopm=0; idiag_rhoprms=0; idiag_rhop2m=0; idiag_rhopmax=0
        idiag_rhopmin=0
        idiag_nparmax=0; idiag_nmigmax=0; idiag_mpt=0
        idiag_npmx=0; idiag_npmy=0; idiag_npmz=0
        idiag_rhopmx=0; idiag_rhopmy=0; idiag_rhopmz=0
        idiag_epspmx=0; idiag_epspmy=0; idiag_epspmz=0
        idiag_rhopmxy=0
      endif
!
!  Run through all possible names that may be listed in print.in
!
      if (lroot .and. ip<14) print*,'rprint_particles: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'nparmax',idiag_nparmax)
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
        call parse_name(iname,cname(iname),cform(iname),'vpxmax',idiag_vpxmax)
        call parse_name(iname,cname(iname),cform(iname),'vpymax',idiag_vpymax)
        call parse_name(iname,cname(iname),cform(iname),'vpzmax',idiag_vpzmax)
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
        call parse_name(iname,cname(iname),cform(iname),'nmigmax',idiag_nmigmax)
        call parse_name(iname,cname(iname),cform(iname),'mpt',idiag_mpt)
        call parse_name(iname,cname(iname),cform(iname), &
            'rhoptilm',idiag_rhoptilm)
        call parse_name(iname,cname(iname),cform(iname), &
            'dedragp',idiag_dedragp)
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
!  check for those quantities for which we want y-averages
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'npmy',idiag_npmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'rhopmy',idiag_npmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'epspmy',idiag_epspmy)
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
!  check for those quantities for which we want xy-averages
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'rhopmxy',idiag_rhopmxy)
      enddo
!
    endsubroutine rprint_particles
!***********************************************************************

endmodule Particles
