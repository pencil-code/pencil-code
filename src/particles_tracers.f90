! $Id: particles_tracers.f90,v 1.20 2006-06-14 15:01:42 ajohan Exp $
!  This module takes care of everything related to tracer particles
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 1
! CPARAM logical, parameter :: lparticles=.true.
! CPARAM logical, parameter :: lparticles_planet=.false.
!
! PENCILS PROVIDED rhop,epsd
!
!***************************************************************
module Particles

  use Cdata
  use Particles_cdata
  use Particles_sub
  use Messages

  implicit none

  include 'particles.h'

  real :: xp0=0.0, yp0=0.0, zp0=0.0, eps_dtog=0.01
  logical :: ldragforce_equi_global_eps=.false.
  logical :: lquadratic_interpolation=.false.
  character (len=labellen), dimension (ninit) :: initxxp='nothing'

  namelist /particles_init_pars/ &
      initxxp, xp0, yp0, zp0, bcpx, bcpy, bcpz, eps_dtog, &
      ldragforce_equi_global_eps, lquadratic_interpolation

  namelist /particles_run_pars/ &
      bcpx, bcpy, bcpz, lquadratic_interpolation

  integer :: idiag_xpm=0, idiag_ypm=0, idiag_zpm=0
  integer :: idiag_xp2m=0, idiag_yp2m=0, idiag_zp2m=0
  integer :: idiag_nparmax=0

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
           "$Id: particles_tracers.f90,v 1.20 2006-06-14 15:01:42 ajohan Exp $")
!
!  Indices for particle position.
!
      ixp=npvar+1
      iyp=npvar+2
      izp=npvar+3
!
!  Increase npvar accordingly.
!
      npvar=npvar+3
!
!  Set indices for auxiliary variables
!
      inp = mvar + naux + 1 + (maux_com - naux_com); naux = naux + 1
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
      if (rhop_tilde==0.0) then
! For stratification, take into account gas present outside the simulation box.
        if (lgrav) then
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
    endsubroutine initialize_particles
!***********************************************************************
    subroutine init_particles(f,fp,ineargrid)
!
!  Initial positions and velocities of tracer particles.
!
!  29-dec-04/anders: coded
!
      use EquationOfState, only: gamma, cs20, beta_glnrho_global
      use General, only: random_number_wrapper
      use Mpicomm, only: stop_it
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (nx) :: eps
      real :: r, p, cs
      integer :: j, k
      logical :: lnothing
!
      intent (inout) :: f
      intent (out) :: fp
!
!  Initial particle position.
!
      lnothing=.false.
      do j=1,ninit

        select case(initxxp(j))

        case ('nothing')
          if (lroot .and. .not. lnothing) print*, 'init_particles: nothing'
          lnothing=.true.

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
              fp(k,izp)=zp0*sqrt(-2*alog(r))*cos(2*pi*p)
            endif
          enddo
          if (nxgrid/=1) &
              fp(1:npar_loc,ixp)=xyz0_loc(1)+fp(1:npar_loc,ixp)*Lxyz_loc(1)
          if (nygrid/=1) &
              fp(1:npar_loc,iyp)=xyz0_loc(2)+fp(1:npar_loc,iyp)*Lxyz_loc(2)
 
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
            eps = rhop_tilde*sum(f(l1:l2,m1:m2,n1:n2,inp))/ &
                sum(f(l1:l2,m1:m2,n1:n2,ilnrho))
          else
            eps = rhop_tilde*sum(f(l1:l2,m1:m2,n1:n2,inp))/ &
                sum(exp(f(l1:l2,m1:m2,n1:n2,ilnrho)))
          endif
!
          if (lroot) &
              print*, 'init_particles: average dust-to-gas ratio=', eps(1)
!  Set gas velocity field.
          do m=m1,m2; do n=n1,n2
            cs=sqrt(cs20)
!  Take either global or local dust-to-gas ratio.
            if (.not. ldragforce_equi_global_eps) then
              if (ldensity_nolog) then
                eps = rhop_tilde*f(l1:l2,m,n,inp)/f(l1:l2,m,n,ilnrho)
              else
                eps = rhop_tilde*f(l1:l2,m,n,inp)/exp(f(l1:l2,m,n,ilnrho))
              endif
            endif

            f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + &
                1/gamma*beta_glnrho_global(1)/(2*(1.0+eps))*cs

          enddo; enddo

        case default
          if (lroot) &
              print*, 'init_particles: No such such value for initxxp: ', &
              trim(initxxp(j))
          call stop_it("")

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
      call boundconds_particles(fp,npar_loc,ipar)
!
!  Map particle positions on the grid.
!
      call map_nearest_grid(f,fp,ineargrid)
      call map_xxp_grid(f,fp,ineargrid)
!
!  Sort particles so that they can be accessed contiguously in the memory.
!
      call sort_particles_imn(fp,ineargrid,ipar)
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
!  15-feb-06/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_epsd)) then
        lpencil_in(i_rho)=.true.
        lpencil_in(i_rhop)=.true.
        lcalc_np=.true.
      endif
      if (lpencil_in(i_rhop)) lcalc_np=.true.
!
    endsubroutine pencil_interdep_particles
!***********************************************************************
    subroutine calc_pencils_particles(f,p)
!
!  Calculate particle pencils.
!
!  15-feb-06/anders: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p
! rhop
      if (lpencil(i_rhop)) p%rhop=rhop_tilde*f(l1:l2,m,n,inp)
! epsd
      if (lpencil(i_epsd)) p%epsd(:,1)=p%rhop/p%rho
!
    endsubroutine calc_pencils_particles
!***********************************************************************
    subroutine dxxp_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of tracer particle position (called from main pencil loop).
!
!  25-apr-06/anders: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: uu
      integer :: k
!
      intent (in) :: f, fp, ineargrid
      intent (inout) :: df, dfp
!
!  Identify module and boundary conditions.
!
      if (headtt) print*,'dxxp_dt: Calculate dxxp_dt'
      if (headtt) then
        print*, 'dxxp_dt: Particles boundary condition bcpx=', bcpx
        print*, 'dxxp_dt: Particles boundary condition bcpy=', bcpy
        print*, 'dxxp_dt: Particles boundary condition bcpz=', bcpz
      endif
!
      if (headtt) print*, 'dxxp_dt: Set rate of change of particle '// &
          'position equal to gas velocity.'
!
!  Interpolate gas velocity to position of particles. 
!  Then set particle velocity equal to the local gas velocity.
!
      if (npar_imn(imn)/=0) then
        do k=k1_imn(imn),k2_imn(imn)
          if (lquadratic_interpolation) then
            call interpolate_quadratic(f,iux,iuz,fp(k,ixp:izp),uu,ineargrid(k,:))
          else
            call interpolate_linear(f,iux,iuz,fp(k,ixp:izp),uu,ineargrid(k,:))
          endif
          if (nxgrid/=1) dfp(k,ixp) = dfp(k,ixp) + uu(1)
          if (nygrid/=1) dfp(k,iyp) = dfp(k,iyp) + uu(2)
          if (nzgrid/=1) dfp(k,izp) = dfp(k,izp) + uu(3)
        enddo
      endif
!
    endsubroutine dxxp_dt_pencil
!***********************************************************************
    subroutine dvvp_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of dust particle velocity (called from main pencil loop).
!
!  20-apr-06/anders: dummy
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      if (NO_WARN) print*, f, df, fp, dfp, p, ineargrid
!
    endsubroutine dvvp_dt_pencil
!***********************************************************************
    subroutine dxxp_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of tracer particle position.
!
!  02-jan-05/anders: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (in) :: f, fp, ineargrid
      intent (inout) :: dfp, df
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_xpm/=0) call sum_par_name(fp(1:npar_loc,ixp),idiag_xpm)
        if (idiag_ypm/=0) call sum_par_name(fp(1:npar_loc,iyp),idiag_ypm)
        if (idiag_zpm/=0) call sum_par_name(fp(1:npar_loc,izp),idiag_zpm)
        if (idiag_xp2m/=0) call sum_par_name(fp(1:npar_loc,ixp)**2,idiag_xp2m)
        if (idiag_yp2m/=0) call sum_par_name(fp(1:npar_loc,iyp)**2,idiag_yp2m)
        if (idiag_zp2m/=0) call sum_par_name(fp(1:npar_loc,izp)**2,idiag_zp2m)
        if (idiag_nparmax/=0) call max_name(npar_loc,idiag_nparmax)
      endif
!
    endsubroutine dxxp_dt
!***********************************************************************
    subroutine dvvp_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of dust particle velocity.
!
!  22-aug-05/anders: dummy
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      if (NO_WARN) print*, f, df, fp, dfp, ineargrid
!
    endsubroutine dvvp_dt
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
!  Calculate power spectra of particle variables.
!
!  01-jan-06/anders: dummy
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if (NO_WARN) print*, f
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
        write(3,*) 'inp=', inp
      endif
!
!  Reset everything in case of reset
!
      if (lreset) then
        idiag_xpm=0; idiag_ypm=0; idiag_zpm=0
        idiag_xp2m=0; idiag_yp2m=0; idiag_zp2m=0
        idiag_nparmax=0; idiag_nmigmax=0
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
        call parse_name(iname,cname(iname),cform(iname),'nparmax',idiag_nparmax)
        call parse_name(iname,cname(iname),cform(iname),'nmigmax',idiag_nmigmax)
      enddo
!
    endsubroutine rprint_particles
!***********************************************************************

endmodule Particles
