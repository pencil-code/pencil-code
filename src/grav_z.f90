! $Id: grav_z.f90,v 1.13 2002-06-10 13:07:14 brandenb Exp $

module Gravity

!
!  Vertical gravity (for convection tests, etc.)
!

  use Cparam
  use Cdata, only: z1,z2,ztop,gravz

  implicit none

  namelist /grav_init_pars/ &
       z1,z2,ztop,gravz

  namelist /grav_run_pars/ &
       gravz

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
           "$RCSfile: grav_z.f90,v $", &
           "$Revision: 1.13 $", &
           "$Date: 2002-06-10 13:07:14 $")
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
!
      use Cdata
      use Sub
      use Slices
!
      real, dimension (mx,my,mz,mvar) :: f,df
!
if (headtt) print*,'duu_dt_grav: gravz=',gravz
      df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + gravz
!
      if(ip==0) print*,f !(keep compiler quiet)
    endsubroutine duu_dt_grav
!***********************************************************************
    subroutine potential(xmn,ymn,zmn,rmn, pot)
!
!  gravity potential
!  21-jan-02/wolf: coded
!
      use Cdata, only: nx,gravz
!
!  the following looks stupid, but otherwise it's not ok,
!  especially in 1-D
!  But for gravz, potential is never used, right??
!
      real, dimension (mx) :: xmn
      real :: ymn,zmn
      real, dimension (nx,1,1) :: rmn, pot
!
      pot = -gravz*zmn
!
      if(ip==0) print*,xmn,ymn,rmn !(keep compiler quiet)
    endsubroutine potential
!***********************************************************************

endmodule Gravity
