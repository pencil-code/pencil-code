! $Id: grav_z.f90,v 1.17 2002-07-08 06:51:51 brandenb Exp $

module Gravity

!  Vertical gravity (for convection in a slab or a full star)
!  (The full star geometry is currently in grav_r, but in may well
!  be possible to migrate it in here.)

  implicit none

  real :: z1,z2,zref=0.,zinfty=1.5
  real :: gravz=-1.
  character (len=30) :: grav_profile='const'

!  The gravity potential must always be negative. However, in an plane
!  atmosphere with constant gravity, the potential goes to zero at
!  some position which is referred to as "zinfty".

!AB: Wolfgang, you should explain here the meaning of z1 and z2,
!AB: as well as zref.

!AB: Nils, could you have a look how in galactic physics (Binney & Tremaine)
!AB: the coefficient in front of .5*z^2 is called (vertical epicyclic frequency?)
!AB: We should introduce that instead of keeping a double meaning of gravz.

  namelist /grav_init_pars/ &
       z1,z2,zref,gravz,grav_profile,zinfty

!  It would be rather unusual to change the profile during the
!  run, but "adjusting" the profile slighly may be quite useful.

  namelist /grav_run_pars/ &
       zref,gravz,grav_profile

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
           "$Id: grav_z.f90,v 1.17 2002-07-08 06:51:51 brandenb Exp $")
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
    subroutine duu_dt_grav(f,df)
!
!  add duu/dt according to gravity
!  (do we need f here?/AB)
!
!  9-jan-02/wolf: coded
! 28-jun-02/axel: added 'linear' gravity profile
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
      elseif (grav_profile=='linear') then
        if (headtt) print*,'duu_dt_grav: linear gravz=',gravz
        df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + gravz*z(n)
      endif
!
      if(ip==0) print*,f !(keep compiler quiet)
    endsubroutine duu_dt_grav
!***********************************************************************
    subroutine potential(xmn,ymn,zmn,pot,grav,rmn)
!
!  gravity potential
!  21-jan-02/wolf: coded
!   8-jul-02/axel: activated and used for initial conditions
!
      use Cdata, only: nx,lroot
!
      real, dimension (nx) :: xmn,pot,r
      real :: ymn,zmn
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
!
      select case(grav_profile)
        case('const')
          pot=abs(gravz)*(zmn-zinfty)
          if(present(grav)) then
            grav(:,1:2)=0.
            grav(:,3)=-abs(gravz)
          endif
        case('linear')
          pot=.5*abs(gravz)*(zmn**2-zinfty**2)
          if(present(grav)) then
            grav(:,1:2)=0.
            grav(:,3)=-abs(gravz)*zmn
          endif
        case('radial')
          if (present(rmn)) then
            r=rmn
          else
            r=sqrt(xmn**2+ymn**2+zmn**2)
          endif
          !(not implemented yet; just laying out the idea)
          pot=.5*gravz*r**2/(1.+r**4)
          if(present(grav)) then
            grav(:,1:3)=0.
          endif
      endselect
      first=.false.
!
    endsubroutine potential
!***********************************************************************

endmodule Gravity
