! $Id: grav_self.f90,v 1.1 2003-03-06 14:25:51 brandenb Exp $

module Gravity

!  self-gravity (solves Poisson equation)

  use Cparam

  implicit none

  interface potential
    module procedure potential_global
    module procedure potential_penc
  endinterface

  real, dimension (mx,my,mz) :: phi
  real :: grav=1.
  integer :: iterations_selfgrav=100
  logical :: ldebug_selfgrav=.false.,lself_gravity=.true.

!  NOTE: the following quantities are needed for compatibility
!  with usage of quantities from grav_z in density.f90

  real :: z1,z2,zref,gravz=-1.,zinfty,zgrav=impossible
  character (len=labellen) :: grav_profile='const'

!  The gravity potential must always be negative. However, in an plane
!  atmosphere with constant gravity, the potential goes to zero at
!  some position which is referred to as "zinfty".

  namelist /grav_init_pars/ &
       iterations_selfgrav,ldebug_selfgrav

!  It would be rather unusual to change the profile during the
!  run, but "adjusting" the profile slighly may be quite useful.

  namelist /grav_run_pars/ &
       iterations_selfgrav,ldebug_selfgrav

  contains

!***********************************************************************
    subroutine register_gravity()
!
!  initialise gravity flags
!
!  9-jan-02/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if(.not. first) call stop_it('register_grav called twice')
      first = .false.
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: grav_self.f90,v 1.1 2003-03-06 14:25:51 brandenb Exp $")
!
      lgrav = .true.
      lgravz = .false.
      lgravr = .false.
!
    endsubroutine register_gravity
!***********************************************************************
    subroutine initialize_gravity()
!
!  Set up some variables for gravity; do nothing in grav_z
!  16-jul-02/wolf: coded
!  22-nov-02/tony: renamed from setup_grav
!
    endsubroutine initialize_gravity
!***********************************************************************
    subroutine init_gg(f,xx,yy,zz)
!
!  initialise gravity; called from start.f90
!  9-jan-02/wolf: coded
!  24-nov-2002: renamed from init_grav to stay consistent
! 
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
! Not doing anything (this might change if we decide to store gg)
!
      if(ip==0) print*,f,xx,yy,zz !(keep compiler quiet)
    endsubroutine init_gg
!***********************************************************************
    subroutine duu_dt_grav(f,df)
!
!  add duu/dt according to gravity
!  (do we need f here?/AB)
!
!  9-jan-02/wolf: coded
! 28-jun-02/axel: added 'linear' gravity profile
! 28-jul-02/axel: added 'const_zero' gravity profile
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,mvar) :: grad_phi
!
      intent(in)  :: f
!
!  different gravity profiles
!
      if (headtt) print*,'duu_dt_grav='
      call grad(phi,grad_phi)
      df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz)-grad_phi
!
      if(ip==0) print*,f(1,1,1,1) !(keep compiler quiet)
    endsubroutine duu_dt_grav
!***********************************************************************
    subroutine potential_global(xx,yy,zz,pot,pot0)
!
!  gravity potential
!  16-jul-02/wolf: coded
!
      use Cdata, only: mx,my,mz
      use Mpicomm
!
      real, dimension (mx,my,mz) :: xx,yy,zz, pot
      real, optional :: pot0
!
      call stop_it("potential_global in grav_z not implemented")
!
      if(ip==0) print*,xx(1,1,1)+yy(1,1,1)+zz(1,1,1), &
           pot(1,1,1),pot0  !(keep compiler quiet)
!
    endsubroutine potential_global
!***********************************************************************
    subroutine potential_penc(xmn,ymn,zmn,pot,pot0,grav,rmn)
!
!  calculates gravity potential and gravitational acceleration
!  on a pencil.
!
!  21-jan-02/wolf: coded
!   8-jul-02/axel: activated and used for initial conditions
!
      use Cdata
      use Sub
!
      real, dimension (nx) :: xmn,pot
      real :: ymn,zmn
      real, optional :: pot0
      real, optional, dimension (nx) :: rmn
      real, optional, dimension (nx,3) :: grav
!
      intent(in) :: xmn,ymn,zmn,rmn
      intent(out) :: pot,grav
!
!  identifier
!
      if (headt) print*,'potential'
!
!  different profiles, calculate also gz=-dpot/dz
!  remember, gravz=-1 (at least negative) for z pointing upwards.
!
      pot=phi(l1:l2,m,n)
      if (present(pot0)) pot0=phi(l1,m1,n1)  !(potential at z=0)
      if (present(grav)) then
        call grad(phi,grav)
      endif
!
      if(ip==0) print*,xmn,ymn,zmn,rmn !(keep compiler quiet)
    endsubroutine potential_penc
!***********************************************************************
    subroutine self_gravity(f)
!
!  calculates gravity potential and gravitational acceleration.
!  Routine is called prior to explicit time-advance via Runge-Kutta .
!
!   1-jan-03/axel: coded
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: resid
      real :: fac,diag,om_diag,om=0.9
      integer :: iter
!
!  identifier
!
      if(lroot.and.headt) print*,'self_gravity'
!
      fac=1./dx**2
      diag=-2.*fac
      om_diag=om/diag
!
!  SOR iterations
!
      do iter=1,iterations_selfgrav
        !
        !  x-direction
        !
        resid(2:mx-1,:,:)=fac*(phi(1:mx-2,:,:)-2*phi(2:mx-1,:,:)+phi(3:mx,:,:))
        resid(1     ,:,:)=fac*(phi(  mx  ,:,:)-2*phi(1     ,:,:)+phi(2   ,:,:))
        resid(  mx  ,:,:)=fac*(phi(  mx-1,:,:)-2*phi(  mx  ,:,:)+phi(1   ,:,:))
        !
        resid=resid-grav*(exp(f(:,:,:,ilnrho))-1.)
        phi=phi-om_diag*resid
        phi=phi-sum(phi(l1:l2,m1:m2,n1:n2))/nw
        !
        if(ldebug_selfgrav) then
          print*,iter,phi(l1,m1,n1)
        endif
      enddo
!
!  debug output
!
      if(ldebug_selfgrav) write(99) iter,phi
!
    endsubroutine self_gravity
!***********************************************************************

endmodule Gravity
