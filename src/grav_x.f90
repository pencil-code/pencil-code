! $Id: grav_x.f90,v 1.2 2004-09-14 11:36:49 ajohan Exp $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Gravity

!  Gravity in the x-direction

  use cparam

  implicit none

  interface potential
    module procedure potential_global
    module procedure potential_penc
    module procedure potential_point
  endinterface
!
!  parameters used throughout entire module
!  xinfty is currently prescribed (=0)
!
  real, dimension(nx) :: gravx_profile,potx_profile
  real :: gravx=0.,xinfty=0.,xgrav=0.,dgravx=0., pot_ratio=1.

! parameters needed for compatibility
!
  real :: z1=0.,z2=1.,zref=0.,gravz=0.,zinfty,zgrav=impossible,nu_epicycle=1.
  real :: lnrho_bot,lnrho_top,ss_bot,ss_top
  real :: grav_const=1.
  real :: g0=0.,r0_pot=0.,kx_gg=1.,ky_gg=1.,kz_gg=1.
  integer :: n_pot=10
  character (len=labellen) :: grav_profile='const'
  logical :: lgravx_gas=.true.,lgravx_dust=.true.,lnumerical_equilibrium=.false.

  namelist /grav_init_pars/ &
       grav_profile,gravx,xgrav,dgravx,pot_ratio,kx_gg,ky_gg,kz_gg

!  It would be rather unusual to change the profile during the
!  run, but "adjusting" the profile slighly may be quite useful.

  namelist /grav_run_pars/ &
       grav_profile,gravx,xgrav,dgravx,pot_ratio,kx_gg,ky_gg,kz_gg, &
       lgravx_gas, lgravx_dust

  ! other variables (needs to be consistent with reset list below)
  integer :: i_curlggrms=0,i_curlggmax=0,i_divggrms=0,i_divggmax=0

  contains

!***********************************************************************
    subroutine register_gravity()
!
!  initialise gravity flags
!
!  12-jun-04/axel: adapted from grav_z
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_gravity: called twice')
      first = .false.
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: grav_x.f90,v 1.2 2004-09-14 11:36:49 ajohan Exp $")
!
      lgrav = .true.
      lgravx = .true.
      lgravz = .false.
      lgravr = .false.
!
    endsubroutine register_gravity
!***********************************************************************
    subroutine initialize_gravity()
!
!  Set up some variables for gravity; do nothing in grav_x
!
!  12-jun-04/axel: coded
!
      use CData
      use Mpicomm, only: stop_it
!
!  Different gravity profiles (gas)
!
      select case (grav_profile)

      case('const')
        if (lroot) print*,'initialize_gravity: constant gravx=',gravx
        gravx_profile=gravx
        potx_profile=-gravx*(x(l1:l2)-xinfty)
!
!  tanh profile
!  for isothermal EOS, we have 0=-cs2*dlnrho+gravx
!  pot_ratio gives the resulting ratio in the density.
!
      case('tanh')
        if (dgravx==0.) call stop_it("initialize_gravity: dgravx=0 not OK")
        if (lroot) print*,'initialize_gravity: tanh profile, gravx=',gravx
        if (lroot) print*,'initialize_gravity: xgrav,dgravx=',xgrav,dgravx
        gravx=-alog(pot_ratio)/dgravx
        gravx_profile=gravx*.5/cosh((x(l1:l2)-xgrav)/dgravx)**2
        potx_profile=-gravx*.5*(1.+tanh((x(l1:l2)-xgrav)/dgravx))*dgravx

      case('sinusoidal')
        if(headtt) print*,'initialize_gravity: sinusoidal grav, gravx=',gravx
        gravx_profile = -gravx*sin(kx_gg*x(l1:l2))

      case default
        if (lroot) print*,'initialize_gravity: grav_profile not defined'

      endselect
!
    endsubroutine initialize_gravity
!***********************************************************************
    subroutine init_gg(f,xx,yy,zz)
!
!  initialise gravity; called from start.f90
!
!  12-jun-04/axel: adapted from grav_z
! 
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
! Not doing anything (this might change if we decide to store gg)
!
      if(ip==0) print*,f,xx,yy,zz !(keep compiler quiet)
    endsubroutine init_gg
!***********************************************************************
    subroutine duu_dt_grav(f,df,uu,rho1)
!
!  add duu/dt according to gravity
!
!  12-jun-04/axel: adapted from grav_z
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: uu
      real, dimension (nx) :: rho1
      integer :: k
!
      intent(in) :: f
!
!  gravity profile in the x direction
!
      if(lhydro .and. lgravx_gas) df(l1:l2,m,n,iux) = &
          df(l1:l2,m,n,iux) + gravx_profile
      if(ldustvelocity .and. lgravx_dust) then
        do k=1,ndustspec
          df(l1:l2,m,n,iudx(k)) = df(l1:l2,m,n,iudx(k)) + gravx_profile
        enddo
      endif
!
      if(ip==0) print*,f,uu,rho1 !(keep compiler quiet)
!        
    endsubroutine duu_dt_grav
!***********************************************************************
    subroutine potential_global(xx,yy,zz,pot,pot0)
!
!  gravity potential
!
!  12-jun-04/axel: adapted from grav_z
!
      use Cdata, only: mx,my,mz
      use Mpicomm
!
      real, dimension (mx,my,mz) :: xx,yy,zz, pot
      real, optional :: pot0
!
      call stop_it("potential_global: not implemented for grav_x")
!
      if(ip==0) print*,xx(1,1,1)+yy(1,1,1)+zz(1,1,1), &
           pot(1,1,1),pot0  !(keep compiler quiet)
!
    endsubroutine potential_global
!***********************************************************************
    subroutine potential_penc(xmn,ymn,zmn,pot,pot0,grav,rmn)
!
!  calculates gravity potential and acceleration on a pencil
!
!  12-jun-04/axel: adapted from grav_z
!
      use Cdata, only: nx,lroot
!
      real, dimension (nx) :: pot,r
      real, optional :: ymn,zmn,pot0
      real, optional, dimension (nx) :: xmn,rmn
      real, optional, dimension (nx,3) :: grav
!
      real :: nu_epicycle2
      logical, save :: first=.true.
!
      intent(in) :: xmn,ymn,zmn,rmn
      intent(out) :: pot
!
!  identifier
!
      if (lroot.and.first) print*,'potential_penc: zinfty=',zinfty
!
!  the different cases are already precalculated in initialize_gravity
!
      pot=potx_profile
!
!  prevent identifier from being called more than once
!
      first=.false.
!
    endsubroutine potential_penc
!***********************************************************************
    subroutine potential_point(x,y,z,r, pot,pot0, grav)
!
!  Gravity potential in one point
!
!  12-jun-04/axel: adapted from grav_z
!
      use Mpicomm, only: stop_it
!
      real :: pot,rad
      real, optional :: x,y,z,r
      real, optional :: pot0,grav
!
      call stop_it("grav_x: potential_point not implemented")
!
      if(ip==0) print*,x,y,z,r,pot,pot0,grav     !(to keep compiler quiet)
    endsubroutine potential_point
!***********************************************************************
    subroutine rprint_gravity(lreset,lwrite)
!
!  reads and registers print parameters relevant for gravity advance
!  dummy routine
!
!  12-jun-04/axel: adapted from grav_z
!
      use Cdata
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  write column, i_XYZ, where our variable XYZ is stored
!  idl needs this even if everything is zero
!
      if (lwr) then
        write(3,*) 'i_curlggrms=',i_curlggrms
        write(3,*) 'i_curlggmax=',i_curlggmax
        write(3,*) 'i_divggrms=',i_divggrms
        write(3,*) 'i_divggmax=',i_divggmax
        write(3,*) 'igg=',igg
        write(3,*) 'igx=',igx
        write(3,*) 'igy=',igy
        write(3,*) 'igz=',igz
      endif
!
      if(ip==0) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_gravity
!***********************************************************************

endmodule Gravity
