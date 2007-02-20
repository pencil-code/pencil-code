! $Id: selfgravity.f90,v 1.21 2007-02-20 17:50:30 dobler Exp $

!
!  This module takes care of self gravity by solving the Poisson equation
!    (d^2/dx^2 + d^2/dy^2 + d^2/dz^2)phi = 4*pi*G*rho
!  for the potential phi.
!

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 1
! COMMUNICATED AUXILIARIES 1
!
! PENCILS PROVIDED potself,gpotself
!
!***************************************************************

module Selfgravity

  use Cdata
  use Cparam
  use Messages

  implicit none

  include 'selfgravity.h'

  real :: rhs_poisson_const=1.0
  real, target :: tstart_selfgrav=0.0

  real :: kmax=0.0
  logical :: lselfgravity_gas=.true., lselfgravity_dust=.false.
  logical :: lklimit_shear=.false., lpoisson_fftxy_discretez=.false.

  namelist /selfgrav_init_pars/ &
      rhs_poisson_const, lselfgravity_gas, lselfgravity_dust, &
      tstart_selfgrav, lpoisson_fftxy_discretez

  namelist /selfgrav_run_pars/ &
      rhs_poisson_const, lselfgravity_gas, lselfgravity_dust, &
      tstart_selfgrav, lklimit_shear, kmax, lpoisson_fftxy_discretez

  integer :: idiag_gpoten=0, idiag_gpotenmxy=0
  integer :: idiag_gpotselfxm=0, idiag_gpotselfym=0, idiag_gpotselfzm=0
  integer :: idiag_gpotselfx2m=0, idiag_gpotselfy2m=0, idiag_gpotselfz2m=0

  contains

!***********************************************************************
    subroutine register_selfgravity()
!
!  Initialise self gravity variables.
!
!  15-may-06/anders+jeff: adapted
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_selfgravity: called twice')
      first = .false.
!
!  Set indices for auxiliary variables
!
      ipotself = mvar + naux_com + 1; naux = naux + 1; naux_com = naux_com + 1
!
!  Identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: selfgravity.f90,v 1.21 2007-02-20 17:50:30 dobler Exp $")
!
!  Put variable name in array
!
      varname(ipotself) = 'potself'
!
!  Set lselfgravity.
!
      lselfgravity=.true.
!
    endsubroutine register_selfgravity
!***********************************************************************
    subroutine initialize_selfgravity()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  15-may-06/anders+jeff: adapted
!
      use SharedVariables
!
      integer :: ierr=0
!
!  Limit the wavenumber to the maximum circular region that is always available
!  in k-space. The radial wavenumber kx changes with time due to shear as
!
!    kx = kx0+qshear*Omega*t*ky
!
!  Considering the available (kx,ky) space, it turns slowly from a square to a
!  parallellogram (the hole for kx,ky<1 is ignored here):
!
!       - - - -                  - - - -
!      |       |               /       /
!      |       |    ->       /       /
!      |       |           /       /
!      |       |         /       /
!       - - - -          - - - -
!
!  To make the gravity force isotropic at small scales one can limit k to
!  the largest circular region that is present in both the square and the
!  parallellogram. The circle has radius kmax=kNy/sqrt(2). See Gammie (2001).
!
      if (lklimit_shear) kmax = kx_ny/sqrt(2.0)
!
      if (.not.lpoisson) then
        if (lroot) print*, 'initialize_selfgravity: must choose a Poisson '// &
            'solver in Makefile.local for self-gravity'
        call fatal_error('initialize_selfgravity','')
      endif
!
      if (lselfgravity_gas.and..not.(lhydro.and.ldensity)) then
        if (lroot) print*, 'initialize_selfgravity: must choose a hydro '// &
            'and a density module in Makefile.local for self-gravity'
        call fatal_error('initialize_selfgravity','')
      endif
!
      if (lselfgravity_dust.and..not.(ldustvelocity.and.ldustdensity)) then
        if (lroot) then
          print*, 'initialize_selfgravity: must choose a dust velocity '// &
              'and a dust density module in Makefile.local for '// &
              'self-gravity on the dust fluid'
        endif
        call fatal_error('initialize_selfgravity','')
      endif
!
!  Share the variable tstart_selfgrav so that it can be used by other
!  self-gravity modules.
!
      call put_shared_variable('tstart_selfgrav',tstart_selfgrav,ierr)
      if (ierr/=0) then
        if (lroot) print*, 'initialize_selfgravity: there was a problem '// &
            'when sharing tstart_selfgrav!'
        call fatal_error('initialize_selfgravity','')
      endif
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
      if (idiag_gpoten/=0 .or. idiag_gpotenmxy/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_potself)=.true.
      endif
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
      if (NO_WARN) print*, lpencil_in !(keep compiler quiet)
!
    endsubroutine pencil_interdep_selfgravity
!***********************************************************************
    subroutine calc_pencils_selfgravity(f,p)
!
!  Calculate Selfgravity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  15-may-06/anders+jeff: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(inout) :: f, p
!
      if (lpencil(i_potself)) p%potself = f(l1:l2,m,n,ipotself)
      if (lpencil(i_gpotself)) then
        call grad(f,ipotself,p%gpotself)
        if (igpotselfx/=0) f(l1:l2,m,n,igpotselfx:igpotselfz)=p%gpotself
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
      use Poisson
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (nx,ny,nz) :: rhs_poisson
!
      if (t>=tstart_selfgrav) then
!
!  Consider self-gravity from gas and dust density or from either one.
!
        if (lhydro.and.ldensity.and.lselfgravity_gas) then
          if (ldensity_nolog) then
            rhs_poisson=rhs_poisson_const*f(l1:l2,m1:m2,n1:n2,ilnrho)
          else
            rhs_poisson=rhs_poisson_const*exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
          endif
        endif
!  Dust.
        if (ldustdensity.and.ldustvelocity.and.lselfgravity_dust) then
          if (ldustdensity_log) then
            if (lselfgravity_gas) then  ! No need to zero rhs.
              rhs_poisson = rhs_poisson + &
                  rhs_poisson_const*exp(f(l1:l2,m1:m2,n1:n2,ind(1)))
            else                        ! Must zero rhs.
              rhs_poisson = rhs_poisson_const*exp(f(l1:l2,m1:m2,n1:n2,ind(1)))
            endif
          else
            if (lselfgravity_gas) then  ! No need to zero rhs.
              rhs_poisson = rhs_poisson + &
                  rhs_poisson_const*f(l1:l2,m1:m2,n1:n2,ind(1))
            else                        ! Must zero rhs.
              rhs_poisson = rhs_poisson_const*f(l1:l2,m1:m2,n1:n2,ind(1))
            endif
          endif
        endif
!
!  Contribution from particles is taken care of by the particle modules.
!
        if (lparticles) &
            call particles_calc_selfpotential(f,rhs_poisson,rhs_poisson_const, &
            lselfgravity_gas.or.lselfgravity_dust)
!
!  Send the right-hand-side of the Poisson equation to the Poisson solver and
!  receive the self-gravity potential back.
!
        if (lpoisson_fftxy_discretez) then
          call poisson_solver_fftxy_discretez(rhs_poisson)
        else
          if (kmax/=0.0) then
            call poisson_solver_fft(rhs_poisson,kmax=kmax)
          else
            call poisson_solver_fft(rhs_poisson)
          endif
        endif
!
!  Put potential into f array.
!
        f(l1:l2,m1:m2,n1:n2,ipotself) = rhs_poisson
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
      use Sub

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(out) :: df
!
!  Add self-gravity acceleration on the gas and on the dust.
!
      if (t>=tstart_selfgrav) then
        if (lhydro.and.lselfgravity_gas) &
            df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - p%gpotself
        if ( ldustvelocity.and.lselfgravity_dust) &
            df(l1:l2,m,n,iudx(1):iudz(1)) = df(l1:l2,m,n,iudx(1):iudz(1)) - p%gpotself
      endif
!
      if (ldiagnos) then
        if (idiag_gpoten/=0) call sum_mn_name(p%potself*p%rho,idiag_gpoten)
        if (idiag_gpotselfxm/=0) &
            call sum_mn_name(p%gpotself(:,1),idiag_gpotselfxm)
        if (idiag_gpotselfym/=0) &
            call sum_mn_name(p%gpotself(:,2),idiag_gpotselfym)
        if (idiag_gpotselfzm/=0) &
            call sum_mn_name(p%gpotself(:,3),idiag_gpotselfzm)
        if (idiag_gpotselfx2m/=0) &
            call sum_mn_name(p%gpotself(:,1)**2,idiag_gpotselfx2m)
        if (idiag_gpotselfy2m/=0) &
            call sum_mn_name(p%gpotself(:,2)**2,idiag_gpotselfy2m)
        if (idiag_gpotselfz2m/=0) &
            call sum_mn_name(p%gpotself(:,3)**2,idiag_gpotselfz2m)
        if (idiag_gpotenmxy/=0) call zsum_mn_name_xy(p%potself*p%rho,idiag_gpotenmxy)
      endif
!
      if (NO_WARN) print*, f, p !(keep compiler quiet)
!
    endsubroutine duu_dt_selfgrav
!***********************************************************************
    subroutine read_selfgravity_init_pars(unit,iostat)
!
!  Read self gravity init parameters
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=selfgrav_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=selfgrav_init_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_selfgravity_init_pars
!***********************************************************************
    subroutine write_selfgravity_init_pars(unit)
!
!  Write self gravity init parameters
!
      integer, intent(in) :: unit
!
      write(unit,NML=selfgrav_init_pars)
!
    endsubroutine write_selfgravity_init_pars
!***********************************************************************
    subroutine read_selfgravity_run_pars(unit,iostat)
!
!  Read self gravity run parameters
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=selfgrav_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=selfgrav_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_selfgravity_run_pars
!***********************************************************************
    subroutine write_selfgravity_run_pars(unit)
!
!  Write self gravity run parameters
!
      integer, intent(in) :: unit
!
      write(unit,NML=selfgrav_run_pars)
!
    endsubroutine write_selfgravity_run_pars
!***********************************************************************
    subroutine rprint_selfgravity(lreset,lwrite)
!
!  reads and registers print parameters relevant for gravity advance
!  dummy routine
!
!  16-may-06/anders+jeff: adapted
!
      use Sub

      logical :: lreset,lwr
      logical, optional :: lwrite
!
      integer :: iname, inamexy
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
      if (lreset) then
        idiag_gpoten=0; idiag_gpotenmxy=0
        idiag_gpotselfxm=0; idiag_gpotselfym=0; idiag_gpotselfzm=0
        idiag_gpotselfx2m=0; idiag_gpotselfy2m=0; idiag_gpotselfz2m=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'gpoten',idiag_gpoten)
        call parse_name(iname,cname(iname),cform(iname),'gpotselfxm', &
            idiag_gpotselfxm)
        call parse_name(iname,cname(iname),cform(iname),'gpotselfym', &
            idiag_gpotselfym)
        call parse_name(iname,cname(iname),cform(iname),'gpotselfzm', &
            idiag_gpotselfzm)
        call parse_name(iname,cname(iname),cform(iname),'gpotselfx2m', &
            idiag_gpotselfx2m)
        call parse_name(iname,cname(iname),cform(iname),'gpotselfy2m', &
            idiag_gpotselfy2m)
        call parse_name(iname,cname(iname),cform(iname),'gpotselfz2m', &
            idiag_gpotselfz2m)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'gpotenmxy', idiag_gpotenmxy)
      enddo
!
!  write column where which variable is stored
!
      if (lwr) then
        write(3,*) 'i_gpoten=',idiag_gpoten
        write(3,*) 'i_gpotenmxy=',idiag_gpotenmxy
        write(3,*) 'ipotself=', ipotself
      endif
!
    endsubroutine rprint_selfgravity
!***********************************************************************


endmodule Selfgravity
