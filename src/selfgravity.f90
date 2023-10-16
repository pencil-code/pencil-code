! $Id$
!
!  This module takes care of self gravity by solving the Poisson equation
!    (d^2/dx^2 + d^2/dy^2 + d^2/dz^2)phi = 4*pi*G*rho
!  for the potential phi.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lselfgravity = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 1
! COMMUNICATED AUXILIARIES 1
!
! PENCILS PROVIDED potself; gpotself(3)
!
!***************************************************************
module Selfgravity
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'selfgravity.h'
!
!  Init Parameters
!
  real, target :: rhs_poisson_const=1.0, gravitational_const=0.0
  real, target :: tstart_selfgrav=0.0
  real, target :: tselfgrav_gentle = 0.0
  real :: kappa=0.0
!
  logical :: lselfgravity_gas=.true., lselfgravity_dust=.false.
  logical :: lselfgravity_neutrals=.false.
!
  namelist /selfgrav_init_pars/ &
      rhs_poisson_const, lselfgravity_gas, lselfgravity_dust, &
      lselfgravity_neutrals, tstart_selfgrav, gravitational_const, kappa
!
!  Run Parameters
!
  logical :: ljeans_stiffening = .false.
  integer :: nj_stiff = 8
  real :: stiff_gamma = 5. / 3.
!
  namelist /selfgrav_run_pars/ &
      rhs_poisson_const, lselfgravity_gas, lselfgravity_dust, &
      lselfgravity_neutrals, tstart_selfgrav, gravitational_const, kappa, &
      ljeans_stiffening, nj_stiff, stiff_gamma, tselfgrav_gentle
!
!  Diagnostic Indices
!
  integer :: idiag_potselfm=0, idiag_rpotselfm=0, idiag_potself2m=0, idiag_potselfmxy=0
  integer :: idiag_potselfmx=0, idiag_potselfmy=0, idiag_potselfmz=0
  integer :: idiag_gpotselfxm=0, idiag_gpotselfym=0, idiag_gpotselfzm=0
  integer :: idiag_gpotselfx2m=0, idiag_gpotselfy2m=0, idiag_gpotselfz2m=0
  integer :: idiag_gxgym=0, idiag_gxgzm=0, idiag_gygzm=0
  integer :: idiag_grgpm=0, idiag_grgzm=0, idiag_gpgzm=0
  integer :: idiag_qtoomre=0,idiag_qtoomremin=0
  integer :: idiag_jeanslength=0, idiag_ljeans2d=0
  integer :: idiag_rugpotselfm=0 ! DIAG_DOC: $\left<\rho\uv\cdot\nabla\Phi\right>$
  integer :: idiag_gpotself2m=0  ! DIAG_DOC: $\left<(\nabla\Phi)^2\right>$
!
!  Module Variables
!
  real, dimension(mz) :: rho0z = 0.0
!
  contains
!***********************************************************************
    subroutine register_selfgravity()
!
!  Initialise self gravity variables.
!
!  15-may-06/anders+jeff: adapted
!
      use FArrayManager
      use SharedVariables, only: put_shared_variable
!
!  Set indices for auxiliary variables
!
      call farray_register_auxiliary('potself',ipotself,communicated=.true.)
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  Share the variable tstart_selfgrav so that it can be used by other self-gravity modules.
!
      call put_shared_variable('tstart_selfgrav',tstart_selfgrav,caller='register_selfgravity')
!
!  Share rhs_poisson_const and gravitational_const.
!
      call put_shared_variable('rhs_poisson_const',rhs_poisson_const)
      call put_shared_variable('gravitational_const',gravitational_const)
!
!  Share tselfgrav_gentle.
!
      call put_shared_variable('tselfgrav_gentle', tselfgrav_gentle)
!
    endsubroutine register_selfgravity
!***********************************************************************
    subroutine initialize_selfgravity(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  15-may-06/anders+jeff: adapted
!
      use EquationOfState, only: get_stratz
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ierr=0
      integer :: i
!
!  Initialize gravitational potential to zero.
!
      f(:,:,:,ipotself)=0.0
!
!  If gravitational constant was set, re-define rhs_poisson_const.
!  else define the gravitational constant via rhs_poisson_const
!
      if (gravitational_const/=0.0) then
        rhs_poisson_const=4*pi*gravitational_const
      else
        gravitational_const=rhs_poisson_const/(4*pi)
      endif
!
      if (.not.lpoisson) &
        call fatal_error('initialize_selfgravity','must choose a Poisson solver in Makefile.local')
!
!  Check that density and self-potential have consistent boundary conditions.
!
      if (ldensity) then
        i = merge(irho, ilnrho, ldensity_nolog)
        if (bcx(ipotself)=='p' .and. .not.(bcx(i)=='p')) then
          if (lroot) then
            print*, 'initialize_selfgravity: potself has bcx=''p'', but the density is not'
            print*, '                        periodic! (you must set a proper boundary condition'
            print*, '                        for the potential)'
            print*, 'initialize_selfgravity: bcx=', bcx
          endif
          call fatal_error('initialize_selfgravity','')
        endif
        if (bcy(ipotself)=='p' .and. .not.(bcy(i)=='p')) then
          if (lroot) then
            print*, 'initialize_selfgravity: potself has bcy=''p'', but the density is not'
            print*, '                        periodic! (you must set a proper boundary condition'
            print*, '                        for the potential)'
            print*, 'initialize_selfgravity: bcy=', bcy
          endif
          call fatal_error('initialize_selfgravity','')
        endif
        if (bcz(ipotself)=='p' .and. .not.(bcz(i)=='p')) then
          if (lroot) then
            print*, 'initialize_selfgravity: potself has bcz=''p'', but the density is not'
            print*, '                        periodic! (you must set a proper boundary condition'
            print*, '                        for the potential)'
            print*, 'initialize_selfgravity: bcz=', bcz
          endif
          call fatal_error('initialize_selfgravity','')
        endif
      endif
!
      if (lneutraldensity) then
        if (bcx(ipotself)=='p' .and. .not.(bcx(i)=='p')) then
          if (lroot) then
            print*, 'initialize_selfgravity: potself has bcx=''p'', but the density is not'
            print*, '                        periodic! (you must set a proper boundary condition'
            print*, '                        for the potential)'
            print*, 'initialize_selfgravity: bcx=', bcx
          endif
          call fatal_error('initialize_selfgravity','')
        endif
        if (bcy(ipotself)=='p' .and. .not.(bcy(i)=='p')) then
          if (lroot) then
            print*, 'initialize_selfgravity: potself has bcy=''p'', but the density is not'
            print*, '                        periodic! (you must set a proper boundary condition'
            print*, '                        for the potential)'
            print*, 'initialize_selfgravity: bcy=', bcy
          endif
          call fatal_error('initialize_selfgravity','')
        endif
        if (bcz(ipotself)=='p' .and. .not.(bcz(i)=='p')) then
          if (lroot) then
            print*, 'initialize_selfgravity: potself has bcz=''p'', but the density is not'
            print*, '                        periodic! (you must set a proper boundary condition'
            print*, '                        for the potential)'
            print*, 'initialize_selfgravity: bcz=', bcz
          endif
          call fatal_error('initialize_selfgravity','')
        endif
      endif
!
!  Initialize the epicycle frequency for calculating Toomre Q.
!
      if (kappa==0.0) kappa=Omega
      if (lroot.and.kappa/=0.0) print*, 'initialize_selfgravity: epicycle frequency kappa = ',kappa
!
!  Get the background density stratification, if any.
!
      if (lstratz) call get_stratz(z, rho0z)
!
    endsubroutine initialize_selfgravity
!***********************************************************************
    subroutine pencil_criteria_selfgravity()
!
!  All pencils that the Selfgravity module depends on are specified here.
!
!  15-may-06/anders+jeff: adapted
!
      lpenc_requested(i_gpotself)=.true.
!
      if (ljeans_stiffening) then
        lpenc_requested(i_fpres)=.true.
        lpenc_requested(i_cs2)=.true.
        lpenc_requested(i_rho)=.true.
      endif
!
      if (idiag_potselfm/=0 .or. idiag_rpotselfm/=0 .or. idiag_potself2m/=0.0 .or. &
          idiag_potselfmx/=0 .or. idiag_potselfmy/=0 .or. idiag_potselfmz/=0) &
          lpenc_diagnos(i_potself)=.true.
!
      if (idiag_grgpm/=0 .or. idiag_grgzm/=0 .or. idiag_gpgzm/=0) then
        lpenc_diagnos(i_pomx)=.true.
        lpenc_diagnos(i_pomy)=.true.
        lpenc_diagnos(i_phix)=.true.
        lpenc_diagnos(i_phiy)=.true.
      endif
!
      if (idiag_qtoomre/=0.or.idiag_qtoomremin/=0.or.idiag_jeanslength/=0.or.idiag_ljeans2d/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_cs2)=.true.
      endif
!
      if (idiag_rugpotselfm/=0) then
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_uu)=.true.
      endif
!
      if (idiag_gpotself2m/=0) then
        lpenc_requested(i_gpotself)=.true.
      endif
!
      if (idiag_potselfmxy/=0) lpenc_diagnos2d(i_potself)=.true.
!
    endsubroutine pencil_criteria_selfgravity
!***********************************************************************
    subroutine pencil_interdep_selfgravity(lpencil_in)
!
!  Interdependency among pencils from the Selfgravity module is specified here.
!
!  15-may-06/anders+jeff: adapted
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_selfgravity
!***********************************************************************
    subroutine calc_pencils_selfgravity(f,p)
!
!  Calculate Selfgravity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  15-may-06/anders+jeff: coded
!  03-feb-11/ccyang: add Jeans stiffening
!
      use Sub, only: grad
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      logical :: first=.true.
      real, save :: gm1, c
!
      intent(inout) :: f, p
!
      if (lpencil(i_potself)) p%potself = f(l1:l2,m,n,ipotself)
      if (lpencil(i_gpotself)) then
        call grad(f,ipotself,p%gpotself)
        if (igpotselfx/=0) f(l1:l2,m,n,igpotselfx:igpotselfz)=p%gpotself
      endif
!
!  Apply Jeans stiffening to the EOS
!
      if (ljeans_stiffening) then
        if (headtt) then
          print*, 'calc_pencils_selfgravity: stiffening is applied to the EOS with '
          print*, 'calc_pencils_selfgravity: ', nj_stiff, ' points per Jeans length and adiabatic index ', stiff_gamma
        endif
        if (first) then
          gm1 = stiff_gamma - 1.
          if (dimensionality == 2) then
            c = gravitational_const * real(nj_stiff) * dxmax
          else
            c = gravitational_const * (real(nj_stiff) * dxmax)**2 / pi
          endif
          first = .false.
        endif
        p%fpres = p%fpres * spread(1. + stiff_gamma * (c * p%rho / p%cs2)**gm1, 2, 3)
      endif
!
    endsubroutine calc_pencils_selfgravity
!***********************************************************************
    subroutine calc_selfpotential(f)
!
!  Calculate the potential of the self gravity.
!
!  15-may-06/anders+jeff: coded
!
      use Particles_main, only: particles_calc_selfpotential
      use FArrayManager
      use Poisson
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension (nx,ny,nz) :: rhs_poisson
!
      integer :: k
      integer :: ipsi_real, ipsi_imag
!
      if (t>=tstart_selfgrav) then
!
!  Consider self-gravity from gas and dust density or from either one.
!
        if (ldensity.and.lselfgravity_gas) then
          if (lstratz) then
            forall(k = n1:n2) rhs_poisson(:,:,k-nghost) = rho0z(k) * (1.0 + f(l1:l2,m1:m2,k,irho))
          elseif (ldensity_nolog) then
            rhs_poisson = f(l1:l2,m1:m2,n1:n2,irho)
          else
            rhs_poisson = exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
          endif
        else
          rhs_poisson = 0.
        endif
!
!  Contribution from dust.
!
        if (ldustdensity.and.lselfgravity_dust) then
          if (ldustdensity_log) then
            if (lselfgravity_gas) then  ! No need to zero rhs
              rhs_poisson = rhs_poisson + exp(f(l1:l2,m1:m2,n1:n2,ind(1)))
            else                        ! Must zero rhs
              rhs_poisson = exp(f(l1:l2,m1:m2,n1:n2,ind(1)))
            endif
          else
            if (lselfgravity_gas) then  ! No need to zero rhs
              rhs_poisson = rhs_poisson + f(l1:l2,m1:m2,n1:n2,ind(1))
            else                        ! Must zero rhs
              rhs_poisson = f(l1:l2,m1:m2,n1:n2,ind(1))
            endif
          endif
        endif
!
!  Contribution from neutrals.
!
        if (lneutraldensity.and.lselfgravity_neutrals) then
          if (lneutraldensity_nolog) then
            if (lselfgravity_gas.or.lselfgravity_dust) then! No need to zero rhs
              rhs_poisson = rhs_poisson + f(l1:l2,m1:m2,n1:n2,irhon)
            else                        ! Must zero rhs
              rhs_poisson = exp(f(l1:l2,m1:m2,n1:n2,ilnrhon))
            endif
          else
            if (lselfgravity_gas.or.lselfgravity_dust) then! No need to zero rhs
              rhs_poisson = rhs_poisson + exp(f(l1:l2,m1:m2,n1:n2,ilnrhon))
            else                        ! Must zero rhs
              rhs_poisson = f(l1:l2,m1:m2,n1:n2,ilnrhon)
            endif
          endif
        endif
!
!  Contribution from particles is taken care of by the particle modules.
!
        if (lparticles) call particles_calc_selfpotential(f,rhs_poisson, &
            lselfgravity_gas.or.lselfgravity_dust)
!
!  Contribution from nonlinear Schroedinger equation
!
        if (lspecial) then
          ipsi_real=farray_index_by_name('psi_real')
          ipsi_imag=farray_index_by_name('psi_imag')
          if (ipsi_real>0 .and. ipsi_imag>0) then
            rhs_poisson=rhs_poisson &
              +f(l1:l2,m1:m2,n1:n2,ipsi_real)**2 &
              +f(l1:l2,m1:m2,n1:n2,ipsi_imag)**2
          endif
        endif
!
!  Send the right-hand-side of the Poisson equation to the Poisson solver and
!  receive the self-gravity potential back.
!
        call inverse_laplacian(rhs_poisson)
!
!  Put potential into f array.
!
        if (tselfgrav_gentle > 0.0 .and. t < tstart_selfgrav + tselfgrav_gentle) then
          f(l1:l2,m1:m2,n1:n2,ipotself) = 0.5 * rhs_poisson_const * &
              (1.0 - cos(pi * (t - tstart_selfgrav) / tselfgrav_gentle)) * rhs_poisson
        else
          f(l1:l2,m1:m2,n1:n2,ipotself) = rhs_poisson_const*rhs_poisson
        endif
!
      endif ! if (t>=tstart_selfgrav) then
!
    endsubroutine calc_selfpotential
!***********************************************************************
    subroutine duu_dt_selfgrav(f,df,p)
!
!  Add self gravity acceleration on gas.
!
!  15-may-06/anders+jeff: coded
!
      use Sub, only: dot_mn, dot2_mn
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  Add self-gravity acceleration on the gas and on the dust.
!
      if (t>=tstart_selfgrav) then
        if (lhydro.and.lselfgravity_gas) &
            df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - p%gpotself
        if (ldustvelocity.and.lselfgravity_dust) &
            df(l1:l2,m,n,iudx(1):iudz(1)) = df(l1:l2,m,n,iudx(1):iudz(1)) - p%gpotself
        if (lneutralvelocity.and.lselfgravity_neutrals) &
            df(l1:l2,m,n,iunx:iunz) = df(l1:l2,m,n,iunx:iunz) - p%gpotself
      endif

      call calc_diagnostics_selfgrav(p)

      call keep_compiler_quiet(f)
!
    endsubroutine duu_dt_selfgrav
!***********************************************************************
    subroutine calc_diagnostics_selfgrav(p)
!
!  Diagnostic averages.
!
      use Diagnostics
      use Sub, only: dot_mn, dot2_mn

      type (pencil_case) :: p

      real, dimension (nx) :: ugpotself, gpotself2

      if (ldiagnos) then

        call sum_mn_name(p%potself,idiag_potselfm)
        if (idiag_potself2m/=0) call sum_mn_name(p%potself**2,idiag_potself2m)
        if (idiag_rpotselfm/=0) call sum_mn_name(p%potself*p%rho,idiag_rpotselfm)
        if (idiag_rugpotselfm/=0) then
          call dot_mn(p%uu,p%gpotself,ugpotself)
          call sum_mn_name(-p%rho*ugpotself,idiag_rugpotselfm)
        endif
        if (idiag_gpotself2m/=0) then
          call dot2_mn(p%gpotself,gpotself2)
          call sum_mn_name(gpotself2,idiag_gpotself2m)
        endif
        call sum_mn_name(p%gpotself(:,1),idiag_gpotselfxm)
        call sum_mn_name(p%gpotself(:,2),idiag_gpotselfym)
        call sum_mn_name(p%gpotself(:,3),idiag_gpotselfzm)
        if (idiag_gpotselfx2m/=0) call sum_mn_name(p%gpotself(:,1)**2,idiag_gpotselfx2m)
        if (idiag_gpotselfy2m/=0) call sum_mn_name(p%gpotself(:,2)**2,idiag_gpotselfy2m)
        if (idiag_gpotselfz2m/=0) call sum_mn_name(p%gpotself(:,3)**2,idiag_gpotselfz2m)
        if (idiag_gxgym/=0) call sum_mn_name(p%gpotself(:,1)*p%gpotself(:,2),idiag_gxgym)
        if (idiag_gxgzm/=0) call sum_mn_name(p%gpotself(:,1)*p%gpotself(:,3),idiag_gxgzm)
        if (idiag_gygzm/=0) call sum_mn_name(p%gpotself(:,2)*p%gpotself(:,3),idiag_gygzm)
        if (idiag_grgpm/=0 .or. idiag_grgzm/=0 .or. idiag_gpgzm/=0) call calc_cylgrav_stresses(p)
        if (idiag_qtoomre/=0) &
             call sum_mn_name(kappa*sqrt(p%cs2)/(gravitational_const*pi*p%rho),idiag_qtoomre)
        if (idiag_qtoomremin/=0) call max_mn_name(-kappa*sqrt(p%cs2)/ &
             (gravitational_const*pi*p%rho),idiag_qtoomremin,lneg=.true.)
        if (idiag_jeanslength/=0) call max_mn_name(-sqrt(pi*p%cs2/ &
            (gravitational_const*p%rho)),idiag_jeanslength,lneg=.true.)
        if (idiag_ljeans2d/=0) call max_mn_name(-p%cs2/ &
            (gravitational_const*p%rho),idiag_ljeans2d,lneg=.true.)
      endif
!
!  1-D averages.
!
      if (l1davgfirst) then
        call yzsum_mn_name_x(p%potself,idiag_potselfmx)
        call xzsum_mn_name_y(p%potself,idiag_potselfmy)
        call xysum_mn_name_z(p%potself,idiag_potselfmz)
      endif
!
!  2-D averages.
!  Note: the name fragment "potselfm" should have contained an "r" for rho.
!
      if (l2davgfirst) then
        if (idiag_potselfmxy/=0) call zsum_mn_name_xy(p%potself*p%rho,idiag_potselfmxy)
      endif
!
    endsubroutine calc_diagnostics_selfgrav
!***********************************************************************
    subroutine calc_cylgrav_stresses(p)
!
!  Calculates cylindrical gravitational stresses in a cartesian box.
!
!  01-jul-07/wlad: coded
!
      use Diagnostics, only: sum_mn_name
!
      real, dimension(nx) :: gpotr,gpotp,gpotz
      type (pencil_case) :: p
!
      gpotr=p%gpotself(:,1)*p%pomx+p%gpotself(:,2)*p%pomy
      gpotp=p%gpotself(:,1)*p%phix+p%gpotself(:,2)*p%phiy
      gpotz=p%gpotself(:,3)
!
      if (idiag_grgpm/=0) call sum_mn_name(gpotr*gpotp,idiag_grgpm)
      if (idiag_grgzm/=0) call sum_mn_name(gpotr*gpotz,idiag_grgzm)
      if (idiag_gpgzm/=0) call sum_mn_name(gpotp*gpotz,idiag_gpgzm)
!
    endsubroutine calc_cylgrav_stresses
!***********************************************************************
    subroutine read_selfgravity_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=selfgrav_init_pars, IOSTAT=iostat)
!
    endsubroutine read_selfgravity_init_pars
!***********************************************************************
    subroutine write_selfgravity_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=selfgrav_init_pars)
!
    endsubroutine write_selfgravity_init_pars
!***********************************************************************
    subroutine read_selfgravity_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=selfgrav_run_pars, IOSTAT=iostat)
!
    endsubroutine read_selfgravity_run_pars
!***********************************************************************
    subroutine write_selfgravity_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=selfgrav_run_pars)
!
    endsubroutine write_selfgravity_run_pars
!***********************************************************************
    subroutine rprint_selfgravity(lreset,lwrite)
!
!  Reads and registers print parameters relevant for gravity advance.
!
!  16-may-06/anders+jeff: adapted
!
      use Diagnostics
      use FArrayManager, only: farray_index_append
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      integer :: iname, inamex, inamey, inamez, inamexy
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_potselfm=0; idiag_rpotselfm=0; idiag_potself2m=0; idiag_potselfmxy=0
        idiag_potselfmx=0; idiag_potselfmy=0; idiag_potselfmz=0
        idiag_gpotselfxm=0; idiag_gpotselfym=0; idiag_gpotselfzm=0
        idiag_gpotselfx2m=0; idiag_gpotselfy2m=0; idiag_gpotselfz2m=0
        idiag_gxgym=0; idiag_gxgzm=0; idiag_gygzm=0
        idiag_grgpm=0; idiag_grgzm=0; idiag_gpgzm=0
        idiag_qtoomre=0; idiag_qtoomremin=0
        idiag_jeanslength=0; idiag_ljeans2d=0
        idiag_rugpotselfm=0; idiag_gpotself2m=0
      endif
!
!  Run through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'rugpotselfm',idiag_rugpotselfm)
        call parse_name(iname,cname(iname),cform(iname),'gpotself2m',idiag_gpotself2m)
        call parse_name(iname,cname(iname),cform(iname),'potselfm', idiag_potselfm)
        call parse_name(iname,cname(iname),cform(iname),'rpotselfm',idiag_rpotselfm)
        call parse_name(iname,cname(iname),cform(iname),'potself2m',idiag_potself2m)
        call parse_name(iname,cname(iname),cform(iname),'gpotselfxm',idiag_gpotselfxm)
        call parse_name(iname,cname(iname),cform(iname),'gpotselfym',idiag_gpotselfym)
        call parse_name(iname,cname(iname),cform(iname),'gpotselfzm',idiag_gpotselfzm)
        call parse_name(iname,cname(iname),cform(iname),'gpotselfx2m',idiag_gpotselfx2m)
        call parse_name(iname,cname(iname),cform(iname),'gpotselfy2m',idiag_gpotselfy2m)
        call parse_name(iname,cname(iname),cform(iname),'gpotselfz2m',idiag_gpotselfz2m)
        call parse_name(iname,cname(iname),cform(iname),'gxgym',idiag_gxgym)
        call parse_name(iname,cname(iname),cform(iname),'gxgzm',idiag_gxgzm)
        call parse_name(iname,cname(iname),cform(iname),'gygzm',idiag_gygzm)
        call parse_name(iname,cname(iname),cform(iname),'grgpm',idiag_grgpm)
        call parse_name(iname,cname(iname),cform(iname),'grgzm',idiag_grgzm)
        call parse_name(iname,cname(iname),cform(iname),'gpgzm',idiag_gpgzm)
        call parse_name(iname,cname(iname),cform(iname),'qtoomre',idiag_qtoomre)
        call parse_name(iname,cname(iname),cform(iname),'qtoomremin',idiag_qtoomremin)
        call parse_name(iname,cname(iname),cform(iname),'jeanslength',idiag_jeanslength)
        call parse_name(iname,cname(iname),cform(iname),'ljeans2d',idiag_ljeans2d)
      enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'potselfmx',idiag_potselfmx)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'potselfmy',idiag_potselfmy)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'potselfmz',idiag_potselfmz)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'potselfmxy',idiag_potselfmxy)
      enddo
!
!  Write column where which variable is stored.
!
      if (lwr) then
        call farray_index_append('ipotself', ipotself)
      endif
!
    endsubroutine rprint_selfgravity
!***********************************************************************
endmodule Selfgravity
