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
           "$Revision: 1.8 $", &
           "$Date: 2002-04-04 17:06:16 $")
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
!      use Global
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: evr,gg
      real, dimension (nx) :: r_mn,g_r
!
!  evr is the radial unit vector
!
      evr(:,1) = x(l1:l2)
      evr(:,2) = y(m)
      evr(:,3) = z(n)
      r_mn = sqrt(x(l1:l2)**2+y(m)**2+z(n)**2) + epsi
      evr = evr / spread(r_mn,2,3)
      g_r = - r_mn * poly( (/ 2*(cpot(1)*cpot(4)-cpot(2)), &
                              3*(cpot(1)*cpot(5)-cpot(3)), &
                              4*cpot(1)*cpot(3), &
                              cpot(5)*cpot(2)-cpot(3)*cpot(4), &
                              2*cpot(2)*cpot(3), &
                              cpot(3)**2  /), r_mn) &
                   / poly( (/ 1., 0., cpot(4), cpot(5), cpot(3) /), r_mn)**2
!      g_r = - r_mn**2 * poly( (/ 3., 0., 1. /), r_mn) &
!                   / poly( (/ 1., 0., 1., 1. /), r_mn)**2
      gg = evr*spread(g_r,2,3)
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
