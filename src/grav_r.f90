module Gravity

!
!  Radial gravity
!

  use Cparam

  implicit none

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
           "$Revision: 1.5 $", &
           "$Date: 2002-01-19 15:50:03 $")
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
!      use Mpicomm
      use Sub
      use Global
!      use Slices
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: er,gg
      real, dimension (nx) :: r,g_r
      real, dimension (5) :: c
!
      c = (/ 5.088, -4.344, 61.36, 10.91, -13.93 /) ! coefficients for pot.
!
      ! Maybe we could get er explicitly, without taking the gradient?
      call grad(rr,1,gg)        ! er = grad(rr) radial unit vector
      r = rr(l1:l2,m,n)         ! There *must* be a way without global rr
      g_r = - r * poly( (/ 2*(c(1)*c(4)-c(2)), &
                            3*(c(1)*c(5)-c(3)), &
                            4*c(1)*c(3), &
                            c(5)*c(2)-c(3)*c(4), &
                            2*c(2)*c(3), &
                            c(3)**2  /), r) &
                 / poly( (/ 1., 0., c(3), c(2) /), r)**2
!      g_r = - r**2 * poly( (/ 3., 0., 1. /), r) &
!                   / poly( (/ 1., 0., 1., 1. /), r)**2
      gg = gg*spread(g_r,2,3)*tanh(t/3.)*0.4
      df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + gg
!
    endsubroutine duu_dt_grav
!***********************************************************************

endmodule Gravity
