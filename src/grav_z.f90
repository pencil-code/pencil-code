! $Id: grav_z.f90,v 1.23 2002-07-29 09:13:22 brandenb Exp $

module Gravity

!  Vertical gravity (for convection in a slab or a full star)
!  (The full star geometry is currently in grav_r, but in may well
!  be possible to migrate it in here.)

  use cparam

  implicit none

  interface potential
    module procedure potential_global
    module procedure potential_penc
  endinterface

!  zref is the height where rho=rho0 and cs2=cs20.
!  For a single polytrope, zinfty (calculated in the
!  density module) is the height where rho=cs2=0.

  integer :: ngrav=10
  real :: z1=0.,z2=1.,zref=0.,gravz=-1.,zinfty,zgrav=1e30
  character (len=labellen) :: grav_profile='const'

!  The gravity potential must always be negative. However, in an plane
!  atmosphere with constant gravity, the potential goes to zero at
!  some position which is referred to as "zinfty".

!  For initlnrho='piecew-poly', z1 and z2 are the interfaces between the
!  three layers of different polytropic exponent:
!
!      z
!      ^
!      |  m = m2 (top layer)
!      |
!  z2  +
!      |
!      |  m = m0 (unstable [main] layer)
!      |
!  z1  +
!      |
!      |  m = m1 (stable bottom layer)
!      |
!

!AB: Nils, could you have a look how in galactic physics (Binney & Tremaine)
!AB: the coefficient in front of .5*z^2 is called (vertical epicyclic
!AB: frequency?)
!AB: We should introduce that instead of keeping a double meaning of gravz.

  namelist /grav_init_pars/ &
       z1,z2,zref,gravz,grav_profile,zgrav

!  It would be rather unusual to change the profile during the
!  run, but "adjusting" the profile slighly may be quite useful.

  namelist /grav_run_pars/ &
       zref,gravz,grav_profile,zgrav

  contains

!***********************************************************************
    subroutine register_grav()
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
      if (.not. first) call stop_it('register_grav called twice')
      first = .false.
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: grav_z.f90,v 1.23 2002-07-29 09:13:22 brandenb Exp $")
!
      lgrav = .true.
      lgravz = .true.
      lgravr = .false.
!
    endsubroutine register_grav
!***********************************************************************
    subroutine init_grav(f,xx,yy,zz)
!
!  initialise gravity; called from start.f90
!  9-jan-02/wolf: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
! Not doing anything (this might change if we decide to store gg)
!
      if(ip==0) print*,f,xx,yy,zz !(keep compiler quiet)
    endsubroutine init_grav
!***********************************************************************
    subroutine setup_grav()
!
!  Set up some variables for gravity; do nothing in grav_z
!  16-jul-02/wolf: coded
!
    endsubroutine setup_grav
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
!
!  different gravity profiles
!
      if (grav_profile=='const') then
        if (headtt) print*,'duu_dt_grav: constant gravz=',gravz
        df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + gravz
!
!  linear gravity profile (for accretion discs)
!
      elseif (grav_profile=='const_zero') then
        if (headtt) print*,'duu_dt_grav: const_zero gravz=',gravz
        if (z(n)<=zgrav) df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+gravz
!
!  linear gravity profile (for accretion discs)
!
      elseif (grav_profile=='linear') then
        if (headtt) print*,'duu_dt_grav: linear gravz=',gravz
        df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + gravz*z(n)
      else
        if(lroot) print*,'no gravity profile given'
      endif
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
!  gravity potential
!  21-jan-02/wolf: coded
!   8-jul-02/axel: activated and used for initial conditions
!
      use Cdata, only: nx,lroot
!
      real, dimension (nx) :: xmn,pot,r
      real :: ymn,zmn
      real, optional :: pot0
      real, optional, dimension (nx) :: rmn
      real, optional, dimension (nx,3) :: grav
      logical, save :: first=.true.
!
      intent(in) :: xmn,ymn,zmn,rmn
      intent(out) :: pot,grav
!
!  identifier
!
      if (lroot.and.first) print*,'potential: zinfty=',zinfty
!
!  different profiles, calculate also gz=-dpot/dz
!  remember, gravz=-1 (at least negative) for z pointing upwards.
!
      select case(grav_profile)
        case('const')
          pot=-gravz*(zmn-zinfty)
          if (present(pot0)) pot0 = gravz*zinfty !(potential at z=0)
          if (present(grav)) then
            grav(:,1:2)=0.
            grav(:,3)=gravz
          endif
!
!  gravity is set to zero above z=zgrav
!
        case('const_zero')
          if(zmn<=zgrav) then
            pot=-gravz*(zmn-zinfty)
            if (present(grav)) then
              grav(:,1:2)=0.
              grav(:,3)=gravz
            endif
          else
            pot=-gravz*(zgrav-zinfty)
            if (present(grav)) grav=0.
          endif
          if (present(pot0)) then !(potential at z=0)
            if(0.<=zgrav) then
              pot0 = gravz*zinfty
            else
              pot0 =-gravz*(zgrav-zinfty)
            endif
          endif
!
!  gravity increases linearly with height (for accretion discs)
!
        case('linear')
          pot=-.5*gravz*(zmn**2-zinfty**2)
          if (present(pot0)) pot0 = .5*gravz*zinfty**2 !(potential at z=0)
          if (present(grav)) then
            grav=0.
            if(zmn<=zgrav) grav(:,3)=gravz*zmn
          endif
!
!  radial profile; not currently implemented
!
        case('radial')
          if (present(rmn)) then
            r=rmn
          else
            r=sqrt(xmn**2+ymn**2+zmn**2)
          endif
          !(not implemented yet; just laying out the idea)
          pot=.5*gravz*r**2/(1.+r**4)
          if (present(pot0)) pot0 = 0. ! potential at z=0
          if (present(grav)) then
            grav(:,1:3)=0.
          endif
      endselect
      first=.false.
!
    endsubroutine potential_penc
!***********************************************************************

endmodule Gravity
