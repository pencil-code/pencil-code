module Gravity

!
!  Radial gravity
!

  use Cparam

  implicit none

  ! coefficients for potential
  real, dimension (5) :: cpot = (/ 5.088, -4.344, 61.36, 10.91, -13.93 /)

  contains

!***********************************************************************
    subroutine register_grav()
!
!  initialise gravity flags
!
!  10-jan-02/wolf: coded
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
!  identify version number
!
      if (lroot) call cvs_id( &
           "$RCSfile: grav_r.f90,v $", &
           "$Revision: 1.7 $", &
           "$Date: 2002-04-03 20:28:36 $")
!
      lgrav = .true.
      lgravz = .false.
      lgravr = .true.
!
    endsubroutine register_grav
!***********************************************************************
    subroutine init_grav(f,init,ampl,xx,yy,zz)
!
!  initialise gravity; called from start.f90
!  10-jan-02/wolf: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: ampl
      integer :: init
!
! Not doing anything (this might change if we decide to store gg)
!
!
    endsubroutine init_grav
!***********************************************************************
    subroutine duu_dt_grav(f,df)
!
!  add duu/dt according to gravity
!
!  10-jan-02/wolf: coded
!
      use Cdata
      use Sub
      use Global
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: er,gg
      real, dimension (nx) :: r,g_r
!
      ! Maybe we could get er explicitly, without taking the gradient?
      call grad(rr,1,gg)        ! er = grad(rr) radial unit vector
      r = rr(l1:l2,m,n)         ! There *must* be a way without global rr
      g_r = - r * poly( (/ 2*(cpot(1)*cpot(4)-cpot(2)), &
                            3*(cpot(1)*cpot(5)-cpot(3)), &
                            4*cpot(1)*cpot(3), &
                            cpot(5)*cpot(2)-cpot(3)*cpot(4), &
                            2*cpot(2)*cpot(3), &
                            cpot(3)**2  /), r) &
                 / poly( (/ 1., 0., cpot(4), cpot(5), cpot(3) /), r)**2
!      g_r = - r**2 * poly( (/ 3., 0., 1. /), r) &
!                   / poly( (/ 1., 0., 1., 1. /), r)**2
      gg = gg*spread(g_r,2,3)
      df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + gg
!
    endsubroutine duu_dt_grav
!***********************************************************************
!    subroutine potential(xmn,ymn,zmn,rmn, pot)
    subroutine potential(rr, pot)
!
!  gravity potential
!  21-jan-02/wolf: coded
!
      use Cdata, only: nx,ny,nz,gravz
      use Sub, only: poly
!
!      real, dimension (nx,1,1) :: xmn,ymn,zmn,rmn, pot
      real, dimension (mx,my,mz) :: rr,pot
!
      pot = - poly((/cpot(1), 0., cpot(2), cpot(3)/), rr) &
            / poly((/1., 0., cpot(4), cpot(5), cpot(3)/), rr)
!
    endsubroutine potential
!***********************************************************************

endmodule Gravity
