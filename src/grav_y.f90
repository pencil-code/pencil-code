! $Id: grav_y.f90,v 1.4 2004-09-16 15:28:58 ajohan Exp $

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

  use Cdata
  use Cparam

  implicit none

  interface potential
    module procedure potential_global
    module procedure potential_penc
    module procedure potential_point
  endinterface
!
!  parameters used throughout entire module
!
  real, dimension(nx) :: gravy_pencil=0.,poty_pencil=0.
  real :: gravy=0.,yinfty=0.,ygrav=0.,dgravy=0.,pot_ratio=1.
  real :: kx_gg=1.,ky_gg=1.,kz_gg=1.
!
! parameters needed for compatibility
!
  real :: z1=0.,z2=1.,zref=0.,gravz=0.,zinfty,zgrav=impossible,nu_epicycle=1.
  real :: lnrho_bot,lnrho_top,ss_bot,ss_top
  real :: grav_const=1.
  real :: g0=0.,r0_pot=0.
  integer :: n_pot=10
  character (len=labellen) :: grav_profile='const'
  logical :: lnumerical_equilibrium=.false.

  namelist /grav_init_pars/ &
      grav_profile,gravy,ygrav,dgravy,pot_ratio,kx_gg,ky_gg,kz_gg, &
      lgravy_gas,lgravy_dust

!  It would be rather unusual to change the profile during the
!  run, but "adjusting" the profile slighly may be quite useful.

  namelist /grav_run_pars/ &
      grav_profile,gravy,ygrav,dgravy,pot_ratio,kx_gg,ky_gg,kz_gg, &
      lgravy_gas,lgravy_dust

  ! other variables (needs to be consistent with reset list below)
  integer :: i_curlggrms=0,i_curlggmax=0,i_divggrms=0,i_divggmax=0

  contains

!***********************************************************************
    subroutine register_gravity()
!
!  initialise gravity flags
!
!  16-sep-04/anders: adapted from grav_x
!
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
           "$Id: grav_y.f90,v 1.4 2004-09-16 15:28:58 ajohan Exp $")
!
      lgrav =.true.
      lgravy=.true.
      lgravy_gas =.true.
      lgravy_dust=.true.
!
    endsubroutine register_gravity
!***********************************************************************
    subroutine initialize_gravity()
!
!  Set up some variables for gravity
!
!  16-sep-04/anders: adapted from grav_x
!
    endsubroutine initialize_gravity
!***********************************************************************
    subroutine init_gg(f,xx,yy,zz)
!
!  initialise gravity; called from start.f90
!
!  16-sep-04/anders: adapted from grav_x
! 
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
! Not doing anything (this might change if we decide to store gg)
!
      if(ip==0) print*,f,xx,yy,zz !(keep compiler quiet)
!        
    endsubroutine init_gg
!***********************************************************************
    subroutine duu_dt_grav(f,df,uu,rho1)
!
!  Add duu/dt from gravity force
!
!  16-sep-04/anders: coded
!
      use Cdata
      use Mpicomm, only: stop_it
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
!  Gravity on the gas and on the dust
!
      if (lhydro .or. ldustdensity) then
!
!  Different gravity profiles
!
        select case (grav_profile)

        case('const')
          if (lroot) print*,'duu_dt_grav: constant gravy=',gravy
          gravy_pencil =  gravy
          poty_pencil  = -gravy*(y(m)-yinfty)
        case('tanh')
!
!  tanh profile
!  for isothermal EOS, we have 0=-cs2*dlnrho+gravy
!  pot_ratio gives the resulting ratio in the density.
!
          if (dgravy==0.) call stop_it("duu_dt_grav: dgravy=0 not OK")
          if (lroot) print*,'duu_dt_grav: tanh profile, gravy=',gravy
          if (lroot) print*,'duu_dt_grav: ygrav,dgravy=',ygrav,dgravy
          gravy=-alog(pot_ratio)/dgravy
          gravy_pencil= gravy*.5/cosh((y(m)-ygrav)/dgravy)**2
          poty_pencil =-gravy*.5*(1.+tanh((y(m)-ygrav)/dgravy))*dgravy

        case('sinusoidal')
          if(headtt) print*,'duu_dt_grav: sinusoidal grav, gravy=',gravy
          gravy_pencil = -gravy*sin(ky_gg*y(m))

        case('sinxsiny')
          if(headtt) print*,'duu_dt_grav: sinusoidal grav, gravy=',gravy
          gravy_pencil = -gravy*sin(kx_gg*x(l1:l2))*sin(ky_gg*y(m))

        case default
          if (lroot) print*,'duu_dt_grav: grav_profile not valid'

        endselect

      endif
!
!  Add gravity acceleration on gas and dust
!
      if(lhydro .and. lgravy_gas) & 
          df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + gravy_pencil
      if(ldustvelocity .and. lgravy_dust) then
        do k=1,ndustspec
          df(l1:l2,m,n,iudy(k)) = df(l1:l2,m,n,iudy(k)) + gravy_pencil
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
!  16-sep-04/anders: adapted from grav_x
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
!  16-sep-04/anders: adapted from grav_x
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
      pot=poty_pencil
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
!  16-sep-04/anders: adapted from grav_x
!
      use Mpicomm, only: stop_it
!
      real :: pot,rad
      real, optional :: x,y,z,r
      real, optional :: pot0,grav
!
      call stop_it("grav_y: potential_point not implemented")
!
      if(ip==0) print*,x,y,z,r,pot,pot0,grav     !(to keep compiler quiet)
    endsubroutine potential_point
!***********************************************************************
    subroutine rprint_gravity(lreset,lwrite)
!
!  reads and registers print parameters relevant for gravity advance
!  dummy routine
!
!  16-sep-04/anders: adapted from grav_x
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
